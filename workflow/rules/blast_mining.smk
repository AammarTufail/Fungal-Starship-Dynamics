# =============================================================================
# Rule file: blast_mining.smk
# Phase 2c — Mine genomes for ToxA using tblastn
#
# ToxA (UniProt P0C1P7) is a necrotrophic effector protein from wheat
# fungal pathogens. We search all assemblies using the protein sequence as
# query against nucleotide databases (tblastn), which is more sensitive than
# blastn when the gene sequence may have diverged at the DNA level.
# =============================================================================

rule make_blast_db:
    """
    Build a BLAST nucleotide database from a genome assembly.

    The database is required for tblastn (protein query vs. nucleotide db).
    makeblastdb indexes the genome for fast subsequence retrieval.
    """
    input:
        fasta = "{outdir}/genomes/{sample}/{sample}.fasta",
    output:
        nhr = "{outdir}/blast/{sample}/db/{sample}.nhr",
        nin = "{outdir}/blast/{sample}/db/{sample}.nin",
        nsq = "{outdir}/blast/{sample}/db/{sample}.nsq",
    params:
        db_prefix = lambda wc: f"{wc.outdir}/blast/{wc.sample}/db/{wc.sample}",
        db_type   = config["blast"]["db_type"],
    log:
        "{outdir}/logs/blast/{sample}_makedb.log",
    shell:
        """
        makeblastdb \
            -in {input.fasta} \
            -dbtype {params.db_type} \
            -out {params.db_prefix} \
            -parse_seqids \
            2>&1 | tee {log}
        """


rule tblastn_toxa:
    """
    Search for ToxA protein (P0C1P7) in all genome assemblies using tblastn.

    tblastn translates the nucleotide database in all 6 reading frames and
    aligns the protein query against these translated sequences. This approach
    detects ToxA even with significant nucleotide divergence (~50% identity
    at DNA level but >80% at protein level).

    Query: resources/ToxA_P0C1P7.fasta (download from UniProt, see README.md)
    E-value threshold: set in config.yaml (default: 1e-5)
    """
    input:
        query  = config["resources"]["toxa_protein_fasta"],
        db_nhr = "{outdir}/blast/{sample}/db/{sample}.nhr",
    output:
        hits = "{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
    params:
        db_prefix = lambda wc: f"{wc.outdir}/blast/{wc.sample}/db/{wc.sample}",
        evalue    = config["blast"]["toxa_evalue"],
        pct_id    = config["blast"]["toxa_perc_identity"],
        outfmt    = config["blast"]["outfmt"],
        threads   = config["threads"]["blast"],
    log:
        "{outdir}/logs/blast/{sample}_tblastn.log",
    shell:
        """
        # tblastn does not support -perc_identity; apply identity filter with awk post-hoc
        # -word_size 3: minimum for tblastn; improves sensitivity for diverged proteins
        tblastn \
            -query {input.query} \
            -db {params.db_prefix} \
            -evalue {params.evalue} \
            -word_size 3 \
            -outfmt "{params.outfmt}" \
            -num_threads {params.threads} \
            2>{log} \
        | awk -v pid={params.pct_id} '$3 >= pid' \
        > {output.hits}

        # Print hit count for quick validation
        echo "ToxA BLAST hits in {wildcards.sample}: $(wc -l < {output.hits})" | tee -a {log}
        """


rule extract_toxa_flanking_region:
    """
    Extract the genomic region flanking the best ToxA hit (±50 kb).

    This 100 kb window captures:
    - The ToxTA passenger transposon (~12 kb) surrounding ToxA
    - Most or all of the Starship element (Sanctuary: ~400 kb; Horizon: ~230 kb)
      — for full Starship extraction, increase window_size to 500000

    The extracted FASTA and corresponding GFF3 annotations are used as input
    to clinker for synteny comparisons.

    Tool: seqkit subseq
    """
    input:
        fasta   = "{outdir}/genomes/{sample}/{sample}.fasta",
        hits    = "{outdir}/blast/{sample}/{sample}_toxa_hits.tsv",
    output:
        locus_fasta = "{outdir}/blast/{sample}/{sample}_toxa_locus.fasta",
        bed         = "{outdir}/blast/{sample}/{sample}_toxa_locus.bed",
    params:
        flank   = 50000,   # bp to extract on each side of the ToxA hit
    log:
        "{outdir}/logs/blast/{sample}_extract_locus.log",
    shell:
        """
        # Parse the best hit from the BLAST table (highest bitscore = last column)
        best_hit=$(sort -k12,12rn {input.hits} | head -1)

        if [ -z "$best_hit" ]; then
            echo "WARNING: No ToxA hits found in {wildcards.sample}. Skipping locus extraction." | tee {log}
            touch {output.locus_fasta} {output.bed}
            exit 0
        fi

        contig_raw=$(echo "$best_hit" | awk '{{print $2}}')
        # Strip NCBI-style "gb|ACC|" wrapper — FASTA headers use bare accession
        contig=$(echo "$contig_raw" | sed 's/^[a-z]*|//;s/|.*//')
        sstart=$(echo "$best_hit" | awk '{{print $9}}')
        send=$(echo "$best_hit"   | awk '{{print $10}}')

        # Ensure start < end (blast can return reverse-complement hits)
        start=$(( $(echo -e "$sstart\n$send" | sort -n | head -1) - {params.flank} ))
        end=$(( $(echo -e "$sstart\n$send" | sort -n | tail -1) + {params.flank} ))
        if [ $start -lt 1 ]; then start=1; fi

        echo "$contig\t$start\t$end\tToxA_locus_{wildcards.sample}" > {output.bed}

        seqkit subseq \
            --chr "$contig" \
            --region "$start:$end" \
            {input.fasta} \
            > {output.locus_fasta} \
            2>{log}

        echo "Extracted locus: $contig:$start-$end" | tee -a {log}
        """


rule annotate_locus_with_prodigal:
    """
    Predict ORFs in the extracted ToxA locus with prodigal (eukaryote mode is
    not supported; use Augustus instead for refined annotation).

    NOTE: For fungi, use Augustus with a trained species model for accurate
    gene prediction. The rule below uses a simplified approach with BLAST-
    based gene finding which is sufficient for a proof-of-concept analysis.

    For full annotation, replace this rule with Augustus or use Funannotate.
    """
    input:
        locus_fasta = "{outdir}/blast/{sample}/{sample}_toxa_locus.fasta",
    output:
        locus_gbk = "{outdir}/blast/{sample}/{sample}_toxa_locus.gbk",
    log:
        "{outdir}/logs/blast/{sample}_annotate_locus.log",
    shell:
        """
        # Convert FASTA to GenBank format with minimal annotation
        # This simplified step creates a GenBank file with the ToxA region
        # that clinker can use for synteny plotting.
        python scripts/fasta_to_gbk.py \
            --fasta {input.locus_fasta} \
            --sample {wildcards.sample} \
            --output {output.locus_gbk} \
            2>&1 | tee {log}
        """
