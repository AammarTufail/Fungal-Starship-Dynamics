# =============================================================================
# Rule file: data_acquisition.smk
# Phase 1 — Download genome assemblies from NCBI FTP
# =============================================================================

# ─── Rule: Download a genome assembly from NCBI FTP ──────────────────────────
rule download_genome:
    """
    Download a single genome assembly directly from the NCBI FTP server using
    scripts/download_genomes.py.

    The script constructs the FTP URL from the accession + assembly_name
    fields in config.yaml, downloads the *_genomic.fna.gz file via wget,
    and decompresses it in-place.

    FTP URL pattern:
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/<d1>/<d2>/<d3>/
            <accession>_<assembly>/<accession>_<assembly>_genomic.fna.gz
    """
    output:
        fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
    params:
        accession     = lambda wc: config["genomes"][wc.sample]["accession"],
        assembly_name = lambda wc: config["genomes"][wc.sample]["assembly_name"],
    log:
        "{outdir}/logs/download/{sample}.log",
    shell:
        """
        python scripts/download_genomes.py \
            --accession {params.accession} \
            --assembly  {params.assembly_name} \
            --outpath   {output.fasta} \
            --log       {log}
        """


rule assembly_stats:
    """
    Generate quick assembly statistics (N50, total size, contig count) using
    assembly-stats. This gives a fast sanity check before the full BUSCO run.

    Output: a plain-text table per sample in results/assembly_stats/
    """
    input:
        fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
    output:
        stats = "{outdir}/assembly_stats/{sample}_stats.txt",
    log:
        "{outdir}/logs/assembly_stats/{sample}.log",
    shell:
        """
        assembly-stats {input.fasta} > {output.stats} 2>{log}
        cat {output.stats}
        """
