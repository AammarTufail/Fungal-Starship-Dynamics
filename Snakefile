# =============================================================================
# Snakefile — ToxTA Starship Diversity in Plant Pathogenic Fungi
# Author:  [Your Name]
# Purpose: Reproducible pipeline for ETH Zurich pre-application analysis.
#          Covers genome acquisition, BUSCO QC, EDTA repeat annotation,
#          ToxA BLAST mining, Starship Y-Rec detection, TSD identification,
#          and clinker synteny visualisation.
#
# Usage:
#   snakemake --use-conda -n               # dry run
#   snakemake --use-conda --cores 8        # local run
#   snakemake --use-conda --profile slurm/ # HPC run
# =============================================================================

configfile: "config/config.yaml"

# ─── Derive genome IDs from config ───────────────────────────────────────────
SAMPLES = list(config["genomes"].keys())
OUTDIR  = config["outdir"]
LOGDIR  = config["logdir"]

# ─── Include modular rule files ───────────────────────────────────────────────
include: "workflow/rules/data_acquisition.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/repeat_annotation.smk"
include: "workflow/rules/blast_mining.smk"
include: "workflow/rules/starship_analysis.smk"
include: "workflow/rules/visualization.smk"

# ─── Master target rule ───────────────────────────────────────────────────────
rule all:
    """
    Top-level rule. Snakemake resolves all dependencies automatically.
    Running `snakemake --cores N` will execute the full pipeline.
    """
    input:
        # Phase 1: genomes downloaded and prepared
        expand("{outdir}/genomes/{sample}/{sample}.fasta",
               outdir=OUTDIR, sample=SAMPLES),

        # Phase 2a: BUSCO quality reports
        expand("{outdir}/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt",
               outdir=OUTDIR, sample=SAMPLES,
               lineage=config["busco"]["lineage"]),

        # Phase 2b: EDTA repeat annotation
        expand("{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.gff3",
               outdir=OUTDIR, sample=SAMPLES),

        # Phase 2c: BLAST hits for ToxA
        expand("{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
               outdir=OUTDIR, sample=SAMPLES),

        # Phase 3a: Y-Recombinase HMM search
        expand("{outdir}/starship/{sample}/{sample}_yrec_hits.tblout",
               outdir=OUTDIR, sample=SAMPLES),

        # Phase 3b: TSD detection
        expand("{outdir}/starship/{sample}/{sample}_tsds.tsv",
               outdir=OUTDIR, sample=SAMPLES),

        # Phase 3c: clinker synteny plot
        "{outdir}/synteny/toxa_locus_clinker.html".format(outdir=OUTDIR),

        # Phase 4: summary report
        "{outdir}/report/summary_report.tsv".format(outdir=OUTDIR),

        # Locus visualisation
        expand("{outdir}/figures/{sample}_toxa_locus.png",
               outdir=OUTDIR, sample=SAMPLES),

        # BUSCO multi-sample summary plot
        "{outdir}/figures/busco_summary.png".format(outdir=OUTDIR),


# ─── Utility rules ────────────────────────────────────────────────────────────

rule clean:
    """
    Remove all output files. Use with caution.
    Run with: snakemake clean --cores 1
    """
    shell:
        "rm -rf {OUTDIR}"  # logs are inside OUTDIR/logs/


rule dag:
    """
    Generate a DAG of the pipeline and render it as PNG.
    Run with: snakemake dag --cores 1
    Requires: graphviz (`mamba install graphviz`)
    """
    shell:
        "snakemake --dag | dot -Tpng > pipeline_dag.png"


rule report:
    """
    Generate an HTML provenance report.
    Run with: snakemake report --cores 1
    """
    shell:
        "snakemake --report results/snakemake_report.html"
