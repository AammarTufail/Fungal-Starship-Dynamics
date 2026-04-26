# =============================================================================
# Rule file: repeat_annotation.smk
# Phase 2b — Full repeatome characterisation with EDTA v2
#
# EDTA (Extensive de novo TE Annotator) is the current gold standard for
# de novo TE annotation in fungi and plants. It integrates multiple discovery
# tools (LTRharvest, LTR_Finder, TIR-Learner, HelitronScanner, RepeatModeler)
# and produces a curated, non-redundant TE library, which is then used to
# annotate the genome with RepeatMasker.
#
# IMPORTANT: EDTA is Linux-only. On macOS, use Docker (see rule edta_docker).
# On HPC clusters (Linux), use the conda rule (edta_conda).
# =============================================================================

import os

# ─── Linux / HPC rule ─────────────────────────────────────────────────────────
if not config.get("use_docker_for_edta", False):
    rule run_edta:
        """
        Run EDTA on a genome assembly (Linux / HPC only).

        Key outputs:
          *.mod.EDTA.TElib.fa    ← curated TE library (FASTA)
          *.mod.EDTA.anno/*.TE.gff3  ← genome-wide TE annotation
          *.mod.EDTA.anno/*.TE.bed   ← BED format for bedtools
          *.mod.EDTA.sum             ← summary table of TE superfamilies

        EDTA runtime: ~30 min–4 hours depending on genome size and --sensitive flag.
        Use --threads 16+ on an HPC node for reasonable runtimes.

        Tool: EDTA v2
        Docs: https://github.com/oushujun/EDTA
        """
        input:
            fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
        output:
            te_gff3 = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.gff3",
            te_lib  = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.TElib.fa",
            summary = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.sum",
        params:
            sensitive = config["edta"]["sensitive"],
            species   = config["edta"]["species"],
            anno      = config["edta"]["anno"],
            evaluate  = config["edta"]["evaluate"],
            workdir   = lambda wc: f"{wc.outdir}/edta/{wc.sample}",
            use_gpu        = config.get("gpu", {}).get("use_gpu", False),
            conda_env_lib  = config.get("gpu", {}).get("conda_env_lib", ""),
            tf_mem_growth  = config.get("gpu", {}).get("tf_memory_growth", True),
        threads:
            config["edta"]["threads"]
        log:
            "{outdir}/logs/edta/{sample}.log",
        shell:
            """
            mkdir -p {params.workdir}
            mkdir -p $(dirname {log})

            # Resolve absolute paths BEFORE changing directory — all relative
            # Snakemake paths become wrong once we cd into /tmp.
            LOGFILE=$(realpath {log})
            OUT_LIB=$(realpath -m {output.te_lib})
            OUT_SUM=$(realpath -m {output.summary})
            OUT_GFF=$(realpath -m {output.te_gff3})
            mkdir -p "$(dirname "$OUT_GFF")"

            # ── GPU / TensorFlow environment setup ────────────────────────────
            # TIR-Learner uses TensorFlow for CNN prediction. CUDA libraries live
            # inside the conda env lib dir and are not on the system LD_LIBRARY_PATH.
            # Prepend the env lib dir so TF can load libcudart / libcudnn / libcublas.
            if [ "{params.use_gpu}" = "True" ] && [ -n "{params.conda_env_lib}" ]; then
                export LD_LIBRARY_PATH="{params.conda_env_lib}:${{LD_LIBRARY_PATH:-}}"
            fi
            # Allow incremental GPU memory growth to avoid OOM on smaller GPUs.
            if [ "{params.tf_mem_growth}" = "True" ]; then
                export TF_FORCE_GPU_ALLOW_GROWTH=true
            fi
            # Silence TF info/warning spam (errors still shown).
            export TF_CPP_MIN_LOG_LEVEL=2

            # EDTA internally uses 'ln -s' to create genome.fasta.mod.
            # This fails on filesystems that do not support symlinks (e.g. NTFS,
            # exFAT external drives). Work around by running entirely in /tmp
            # and copying the final outputs back.
            EDTA_TMP=$(mktemp -d /tmp/edta_{wildcards.sample}_XXXXXX)
            trap "rm -rf \"$EDTA_TMP\"" EXIT

            cp {input.fasta} "$EDTA_TMP/{wildcards.sample}.fasta"
            cd "$EDTA_TMP"

            EDTA.pl \
                --genome {wildcards.sample}.fasta \
                --species {params.species} \
                --sensitive {params.sensitive} \
                --anno {params.anno} \
                --evaluate {params.evaluate} \
                --force 1 \
                --threads {threads} \
                2>&1 | tee "$LOGFILE"

            echo "EDTA complete for {wildcards.sample}" | tee -a "$LOGFILE"

            # Copy outputs back using absolute destination paths
            cp "$EDTA_TMP/{wildcards.sample}.fasta.mod.EDTA.TElib.fa" "$OUT_LIB"
            cp "$EDTA_TMP/{wildcards.sample}.fasta.mod.EDTA.anno/{wildcards.sample}.fasta.mod.EDTA.TEanno.sum" \
               "$OUT_SUM"
            cp "$EDTA_TMP/{wildcards.sample}.fasta.mod.EDTA.anno/{wildcards.sample}.fasta.mod.EDTA.TEanno.gff3" \
               "$OUT_GFF"
            """

# ─── macOS Docker rule ────────────────────────────────────────────────────────
else:
    rule run_edta:
        """
        Run EDTA via Docker (macOS compatible).

        Requires Docker Desktop to be running.
        Image: oushujun/edta:2.2.0

        Pull first:
            docker pull oushujun/edta:2.2.0
        """
        input:
            fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
        output:
            te_gff3 = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.gff3",
            te_lib  = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.TElib.fa",
            summary = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.sum",
        params:
            sensitive  = config["edta"]["sensitive"],
            species    = config["edta"]["species"],
            anno       = config["edta"]["anno"],
            evaluate   = config["edta"]["evaluate"],
            docker_img = config["edta_docker_image"],
            workdir    = lambda wc: f"{wc.outdir}/edta/{wc.sample}",
            abs_workdir = lambda wc, input: os.path.abspath(f"{wc.outdir}/edta/{wc.sample}"),
        threads:
            config["edta"]["threads"]
        log:
            "{outdir}/logs/edta/{sample}.log",
        shell:
            """
            mkdir -p {params.workdir}
            cp {input.fasta} {params.workdir}/{wildcards.sample}.fasta

            docker run --rm \
                -v {params.abs_workdir}:/data \
                {params.docker_img} \
                EDTA.pl \
                    --genome /data/{wildcards.sample}.fasta \
                    --species {params.species} \
                    --sensitive {params.sensitive} \
                    --anno {params.anno} \
                    --evaluate {params.evaluate} \
                    --force 1 \
                    --threads {threads} \
                2>&1 | tee {log}
            """


rule parse_edta_summary:
    """
    Parse the EDTA summary file and extract per-superfamily repeat content
    as a clean TSV table. This is used for the summary report and figures.

    Outputs a table with columns:
        Sample | TE_Superfamily | Total_bp | Genome_Fraction_pct
    """
    input:
        summary = "{outdir}/edta/{sample}/{sample}.fasta.mod.EDTA.anno/{sample}.fasta.mod.EDTA.TEanno.sum",
        fasta   = "{outdir}/genomes/{sample}/{sample}.fasta",
    output:
        tsv = "{outdir}/edta/{sample}/{sample}_repeat_summary.tsv",
    log:
        "{outdir}/logs/edta/{sample}_parse.log",
    shell:
        """
        python scripts/parse_edta_summary.py \
            --summary {input.summary} \
            --sample {wildcards.sample} \
            --output {output.tsv} \
            2>&1 | tee {log}
        """
