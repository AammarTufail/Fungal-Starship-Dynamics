#!/usr/bin/env python3
"""
fasta_to_gbk.py — Convert a FASTA sequence to a minimal GenBank file for clinker.

clinker requires GenBank or FASTA+GFF3 input. This script produces a minimal
GenBank file from a FASTA sequence. If a GFF3 annotation file is provided,
features are included in the GenBank record.

Usage:
    # Minimal (no annotation):
    python fasta_to_gbk.py \\
        --fasta results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_locus.fasta \\
        --sample Pnodorum_SN15 \\
        --output results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_locus.gbk

    # With GFF3 annotation:
    python fasta_to_gbk.py \\
        --fasta results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_locus.fasta \\
        --gff3 results/edta/Pnodorum_SN15/...TE.gff3 \\
        --sample Pnodorum_SN15 \\
        --output results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_locus.gbk
"""

import argparse
import sys
from datetime import date
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def gff3_to_features(gff3_path: str, contig_id: str, seq_len: int) -> list[SeqFeature]:
    """
    Parse a GFF3 file and return BioPython SeqFeature objects for a given contig.
    Only features that fall within [0, seq_len] are included.
    """
    features = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            if chrom != contig_id:
                continue
            start_i = int(start) - 1  # GFF3 is 1-based; BioPython uses 0-based
            end_i   = int(end)
            if start_i < 0 or end_i > seq_len:
                continue
            strand_int = 1 if strand == "+" else (-1 if strand == "-" else 0)
            # Parse attributes (key=value;... format)
            qualifiers = {}
            for attr in attrs.split(";"):
                attr = attr.strip()
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    qualifiers[k] = [v]
            # Use 'gene' for gene features, 'repeat_region' for TEs
            gb_ftype = "repeat_region" if "te" in ftype.lower() else ftype
            features.append(
                SeqFeature(
                    location=FeatureLocation(start_i, end_i, strand=strand_int),
                    type=gb_ftype,
                    qualifiers=qualifiers,
                )
            )
    return features


def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to GenBank for clinker.")
    parser.add_argument("--fasta",  required=True, help="Input FASTA file")
    parser.add_argument("--gff3",   default=None,  help="Optional GFF3 annotation")
    parser.add_argument("--sample", required=True, help="Sample name (used as organism)")
    parser.add_argument("--output", required=True, help="Output GenBank (.gbk) file")
    args = parser.parse_args()

    records_in = list(SeqIO.parse(args.fasta, "fasta"))
    if not records_in:
        print(f"WARNING: Empty FASTA {args.fasta}. Writing empty GenBank.", file=sys.stderr)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text("")
        return

    records_out = []
    for rec in records_in:
        new_rec = SeqRecord(
            seq=rec.seq,
            id=rec.id,
            name=rec.id[:16],  # GenBank name field max 16 chars
            description=f"{args.sample} ToxA locus",
            annotations={
                "molecule_type": "DNA",
                "organism":      args.sample.replace("_", " "),
                "date":          date.today().strftime("%d-%b-%Y").upper(),
                "data_file_division": "PLN",
            },
        )
        # Add features from GFF3 if provided
        if args.gff3:
            feats = gff3_to_features(args.gff3, rec.id, len(rec.seq))
            new_rec.features.extend(feats)
            print(f"  Added {len(feats)} features from GFF3 to {rec.id}", file=sys.stderr)

        records_out.append(new_rec)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out_fh:
        SeqIO.write(records_out, out_fh, "genbank")

    print(f"Written {len(records_out)} records to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
