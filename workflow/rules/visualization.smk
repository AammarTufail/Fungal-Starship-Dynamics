# =============================================================================
# Rule file: visualization.smk
# Phase 4 — Synteny analysis and locus visualisation
#
# Two complementary visualisation approaches:
#   1. clinker  — gene-level synteny across the ToxA locus in all 3 species
#   2. pyCirclize / dna-features-viewer — nested transposon architecture diagram
# =============================================================================

rule run_clinker:
    """
    Run clinker to generate an interactive HTML synteny plot of the ToxA locus.

    clinker aligns annotated loci (GenBank format), identifies orthologous
    gene pairs, and produces a publication-quality interactive plot showing
    gene conservation and order across species.

    Input: GenBank files for the extracted ToxA locus (±50 kb) from each sample.
    Output: Interactive HTML file viewable in any web browser.

    The key biological question visualised:
    - Is the gene *order* around ToxA conserved? (synteny = shared origin)
    - Are the flanking Starship "cargo" genes the same or different?
      Different flanking genes = independent capture events

    Tool: clinker v0.0.28
    Docs: https://github.com/gamcil/clinker
    """
    input:
        gbk_files = expand(
            "{outdir}/blast/{sample}/{sample}_toxa_locus.gbk",
            outdir=config["outdir"], sample=SAMPLES
        ),
    output:
        html = "{outdir}/synteny/toxa_locus_clinker.html",
        matrix = "{outdir}/synteny/toxa_locus_clinker_matrix.csv",
    params:
        identity = config["clinker"]["min_identity"],
        outdir   = lambda wc: f"{wc.outdir}/synteny",
    log:
        "{outdir}/logs/visualization/clinker.log",
    shell:
        """
        mkdir -p {params.outdir}

        clinker {input.gbk_files} \
            --identity {params.identity} \
            --plot {output.html} \
            --output {output.matrix} \
            2>&1 | tee {log}

        echo "Clinker synteny plot: {output.html}" | tee -a {log}
        echo "Open in browser: python -m http.server 8080 (then visit localhost:8080/$(basename {output.html}))" \
            | tee -a {log}
        """


rule plot_toxa_locus:
    """
    Generate a linear locus diagram showing the nested transposon architecture:
        [Starship boundary] — [ToxTA transposon] — [ToxA gene] — [ToxTA] — [Starship]

    Uses dna-features-viewer and matplotlib for clean, publication-ready figures.
    The script reads the BLAST hit coordinates and EDTA TE annotations to
    accurately represent the nested structure.

    Output: One PNG per sample in results/figures/
    """
    input:
        locus_fasta = "{outdir}/blast/{sample}/{sample}_toxa_locus.fasta",
        te_gff3     = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.gff3",
        blast_hits  = "{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
    output:
        figure = "{outdir}/figures/{sample}_toxa_locus.png",
    log:
        "{outdir}/logs/visualization/{sample}_locus_plot.log",
    shell:
        """
        python scripts/plot_toxa_locus.py \
            --fasta {input.locus_fasta} \
            --te-gff3 {input.te_gff3} \
            --blast-hits {input.blast_hits} \
            --sample {wildcards.sample} \
            --output {output.figure} \
            2>&1 | tee {log}
        """


rule aggregate_report:
    """
    Aggregate all results into a single summary TSV report covering:
    - Assembly statistics (size, N50, contig count)
    - BUSCO completeness scores
    - Repeatome composition (% of genome per TE superfamily)
    - ToxA BLAST hit coordinates and identity
    - Y-Recombinase hit count (Starship proxy)
    - TSD sequences found

    This table is the central output of the pipeline.
    """
    input:
        assembly_stats = expand(
            "{outdir}/assembly_stats/{sample}_stats.txt",
            outdir=config["outdir"], sample=SAMPLES
        ),
        busco_summaries = expand(
            "{outdir}/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt",
            outdir=config["outdir"], sample=SAMPLES,
            lineage=config["busco"]["lineage"]
        ),
        repeat_summaries = expand(
            "{outdir}/edta/{sample}/{sample}_repeat_summary.tsv",
            outdir=config["outdir"], sample=SAMPLES
        ),
        blast_hits = expand(
            "{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
            outdir=config["outdir"], sample=SAMPLES
        ),
        yrec_hits = expand(
            "{outdir}/starship/{sample}/{sample}_yrec_hits.tblout",
            outdir=config["outdir"], sample=SAMPLES
        ),
    output:
        report = "{outdir}/report/summary_report.tsv",
    params:
        samples = SAMPLES,
    log:
        "{outdir}/logs/visualization/aggregate_report.log",
    shell:
        """
        mkdir -p $(dirname {output.report})

        python scripts/summarize_results.py \
            --assembly-stats {input.assembly_stats} \
            --busco {input.busco_summaries} \
            --repeats {input.repeat_summaries} \
            --blast {input.blast_hits} \
            --yrec {input.yrec_hits} \
            --output {output.report} \
            2>&1 | tee {log}

        echo "Summary report written to {output.report}" | tee -a {log}
        """
