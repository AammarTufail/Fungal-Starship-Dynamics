# =============================================================================
# Rule file: starship_analysis.smk
# Phase 3 — Identify Starship elements via Y-Recombinase HMM search
#           and Target Site Duplication (TSD) detection
#
# Starship elements are defined by the presence of a "captain" gene encoding
# a tyrosine-family site-specific recombinase (Tyr-Rec / DUF3435). They are
# the largest fungal transposons known (100 kb–700 kb) and are responsible
# for horizontal transfer of cargo genes including ToxA (via ToxTA).
#
# References:
#   Gluck-Thaler & Vogan 2022 (eLife): original Starship discovery
#   Vogan et al. 2023: Starship diversity and ToxA HGT
# =============================================================================

rule translate_genome:
    """
    Translate all ORFs in the genome to protein FASTA using seqkit translate.
    This is used as the target database for HMMER protein searches.

    Alternative: run a gene prediction tool (Augustus, Funannotate) first
    and use the predicted proteome. Protein prediction is more accurate but
    requires a trained species model.
    """
    input:
        fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
    output:
        proteins = "{outdir}/starship/{sample}/{sample}_sixframe.faa",
    log:
        "{outdir}/logs/starship/{sample}_translate.log",
    shell:
        """
        seqkit translate \
            --frame 6 \
            --trim \
            --min-len 50 \
            {input.fasta} \
            > {output.proteins} \
            2>{log}
        echo "Translated ORFs: $(grep -c '>' {output.proteins})" | tee -a {log}
        """


rule hmmsearch_yrec:
    """
    Search for Starship captain proteins (Y-Recombinase / DUF3435) using HMMER3.

    The DUF3435 domain (Pfam PF11976) is the defining feature of Starship
    elements. Any genome with a DUF3435-containing protein is a candidate
    for harbouring a Starship transposon.

    HMM profile source: Pfam PF11976 (download instructions in README.md)
    E-value: 1e-10 (stringent — DUF3435 is a large, well-defined domain)

    Key output columns (tblout format):
        target_name, accession, query_name, E-value, score, bias
    """
    input:
        proteins  = "{outdir}/starship/{sample}/{sample}_sixframe.faa",
        hmm_profile = config["resources"]["yrec_hmm_profile"],
    output:
        tblout   = "{outdir}/starship/{sample}/{sample}_yrec_hits.tblout",
        domtblout = "{outdir}/starship/{sample}/{sample}_yrec_hits.domtblout",
    params:
        evalue   = config["hmmer"]["yrec_evalue"],
    threads:
        config["threads"]["default"]
    log:
        "{outdir}/logs/starship/{sample}_hmmsearch.log",
    shell:
        """
        hmmsearch \
            --cpu {threads} \
            -E {params.evalue} \
            --tblout {output.tblout} \
            --domtblout {output.domtblout} \
            {input.hmm_profile} \
            {input.proteins} \
            2>&1 | tee {log}

        echo "Y-Rec hits in {wildcards.sample}: $(grep -v '^#' {output.tblout} | wc -l)" \
            | tee -a {log}
        """


rule find_tsds:
    """
    Detect Target Site Duplications (TSDs) flanking putative Starship insertions.

    Starship elements insert into specific genomic sites and leave behind short
    (4–6 bp) direct repeats at both ends — the TSDs. Identifying TSDs:
    1. Confirms the boundaries of the transposon
    2. Identifies the insertion site sequence (host genome context)
    3. Can reveal whether two insertions share the same target site preference

    Strategy:
    1. Use HMMER hits to identify the approximate location of Y-Rec ORFs
    2. Extend the search window (±500 bp from predicted Starship boundaries)
    3. Scan for perfect or near-perfect direct repeats of 4–6 bp

    This rule calls the custom Python script scripts/find_tsds.py
    """
    input:
        fasta  = "{outdir}/genomes/{sample}/{sample}.fasta",
        tblout = "{outdir}/starship/{sample}/{sample}_yrec_hits.tblout",
    output:
        tsds = "{outdir}/starship/{sample}/{sample}_tsds.tsv",
    params:
        window    = config["tsd"]["window_size"],
        min_len   = config["tsd"]["min_tsd_len"],
        max_len   = config["tsd"]["max_tsd_len"],
        max_mm    = config["tsd"]["max_mismatches"],
    log:
        "{outdir}/logs/starship/{sample}_tsds.log",
    shell:
        """
        python scripts/find_tsds.py \
            --fasta {input.fasta} \
            --hmmer-hits {input.tblout} \
            --window {params.window} \
            --min-tsd {params.min_len} \
            --max-tsd {params.max_len} \
            --max-mismatches {params.max_mm} \
            --output {output.tsds} \
            2>&1 | tee {log}
        """


rule compare_starship_structure:
    """
    Compare Starship structural properties across all samples.

    Produces a summary TSV comparing:
    - Presence/absence of Y-Rec in each sample
    - Number of Y-Rec copies (one per Starship in the genome)
    - TSD sequences found at Starship boundaries
    - Estimated Starship size (distance between TSDs)
    - Co-localisation with ToxA BLAST hits

    This table is the primary evidence for whether ToxA was captured once
    (single origin, shared Starship) or independently multiple times
    (convergent capture).
    """
    input:
        yrec_hits = expand(
            "{outdir}/starship/{sample}/{sample}_yrec_hits.tblout",
            outdir=config["outdir"], sample=SAMPLES
        ),
        tsds = expand(
            "{outdir}/starship/{sample}/{sample}_tsds.tsv",
            outdir=config["outdir"], sample=SAMPLES
        ),
        blast_hits = expand(
            "{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
            outdir=config["outdir"], sample=SAMPLES
        ),
    output:
        comparison = "{outdir}/starship/starship_structure_comparison.tsv",
    log:
        "{outdir}/logs/starship/structure_comparison.log",
    shell:
        """
        python scripts/summarize_results.py \
            --yrec-hits {input.yrec_hits} \
            --tsds {input.tsds} \
            --blast-hits {input.blast_hits} \
            --output {output.comparison} \
            2>&1 | tee {log}
        """
