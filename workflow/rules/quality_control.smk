# =============================================================================
# Rule file: quality_control.smk
# Phase 2a — Genome completeness assessment with BUSCO v5
# =============================================================================

rule run_busco:
    """
    Assess genome completeness using BUSCO v5 against the fungi_odb10 lineage.

    BUSCO searches for single-copy orthologous genes expected to be present in
    nearly all fungal genomes. A high complete BUSCOs percentage (>90%) indicates
    a high-quality, near-complete assembly. We expect all three target species
    to have >85% complete BUSCOs for the analyses downstream to be meaningful.

    Tool: BUSCO v5
    Docs: https://busco.ezlab.org/
    Lineage: fungi_odb10 (~758 orthologs, ODB release 10)

    Key output:
        short_summary.specific.<lineage>.<sample>.txt  ← the human-readable summary
        run_<lineage>/full_table.tsv                   ← full BUSCO gene table
    """
    input:
        fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
    output:
        summary = "{outdir}/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt",
        full_table = "{outdir}/busco/{sample}/run_{lineage}/full_table.tsv",
    params:
        lineage        = config["busco"]["lineage"],
        mode           = config["busco"]["mode"],
        augustus_sp    = config["busco"]["augustus_species"],
        out_path       = lambda wc: f"{wc.outdir}/busco",
        out_name       = lambda wc: wc.sample,
    threads:
        config["threads"]["heavy"]
    log:
        "{outdir}/logs/busco/{sample}_{lineage}.log",
    shell:
        """
        busco \
            --in {input.fasta} \
            --out {params.out_name} \
            --out_path {params.out_path} \
            --lineage_dataset {params.lineage} \
            --mode {params.mode} \
            --augustus_species {params.augustus_sp} \
            --cpu {threads} \
            --force \
            2>&1 | tee {log}
        """


rule busco_summary_plot:
    """
    Generate a multi-sample BUSCO summary barplot using the built-in
    `generate_plot.py` script from BUSCO.

    This produces the canonical stacked bar chart (Complete single-copy /
    Complete duplicated / Fragmented / Missing) for all samples at once,
    which is the standard figure used in genome assembly papers.

    Run this rule AFTER all per-sample BUSCO runs are complete.
    """
    input:
        summaries = expand(
            "{outdir}/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt",
            outdir=config["outdir"],
            sample=SAMPLES,
            lineage=config["busco"]["lineage"],
        ),
    output:
        plot = "{outdir}/figures/busco_summary.png",
    params:
        summary_dir = lambda wc: f"{wc.outdir}/busco",
        fig_dir     = lambda wc: f"{wc.outdir}/figures",
    log:
        "{outdir}/logs/busco/summary_plot.log",
    shell:
        """
        # Copy all summaries into one directory for the plotting script
        mkdir -p {params.summary_dir}/all_summaries

        for f in {input.summaries}; do
            cp "$f" {params.summary_dir}/all_summaries/
        done

        generate_plot.py \
            --working_directory {params.summary_dir}/all_summaries \
            2>&1 | tee {log}

        # Move the output plot to the figures directory
        mkdir -p {params.fig_dir}
        mv {params.summary_dir}/all_summaries/busco_figure.png {output.plot}
        echo "BUSCO summary plot saved to {output.plot}"
        """
