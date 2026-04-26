#!/usr/bin/env python3
"""
find_tsds.py — Target Site Duplication (TSD) finder for Starship transposons.

Starship elements generate TSDs of 4–6 bp when they insert into a new genomic
location. The TSD is a short direct repeat found at both the 5' and 3' ends of
the insertion. Identifying TSDs:
  1. Confirms the boundaries of a Starship element
  2. Identifies the preferred insertion-site sequence motif
  3. Distinguishes genuine Starship insertions from false Y-Rec positives

Strategy:
  - Parse HMMER tblout to locate Y-Recombinase-containing contigs/positions
  - Extract flanking sequence windows around candidate boundaries
  - Slide through the window looking for direct repeats of length min_tsd..max_tsd
  - Report all candidate TSDs with their positions

Usage:
    python find_tsds.py \\
        --fasta results/genomes/Pnodorum_SN15/Pnodorum_SN15.fasta \\
        --hmmer-hits results/starship/Pnodorum_SN15/Pnodorum_SN15_yrec_hits.tblout \\
        --window 500 \\
        --min-tsd 4 \\
        --max-tsd 6 \\
        --max-mismatches 0 \\
        --output results/starship/Pnodorum_SN15/Pnodorum_SN15_tsds.tsv
"""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_hmmer_tblout(tblout_path: str) -> list[dict]:
    """
    Parse HMMER3 tblout format and return a list of significant hits.

    Each hit is a dict with: target_name, e_value, score, description.
    Lines starting with '#' are comments.
    """
    hits = []
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            hits.append({
                "target_name": fields[0],
                "accession":   fields[1],
                "query_name":  fields[2],
                "e_value":     float(fields[4]),
                "score":       float(fields[5]),
                "description": " ".join(fields[18:]) if len(fields) > 18 else "",
            })
    return hits


def find_direct_repeats(seq: str, min_len: int = 4, max_len: int = 6,
                        max_mismatches: int = 0) -> list[dict]:
    """
    Find all direct repeats of length min_len..max_len within seq.

    A direct repeat is a pair of identical (or near-identical) subsequences
    separated by some distance. We scan all possible positions and lengths.

    Returns a list of dicts: {seq, start1, end1, start2, end2, length, mismatches}
    """
    results = []
    n = len(seq)

    for tsd_len in range(min_len, max_len + 1):
        for i in range(n - tsd_len * 2):
            candidate = seq[i : i + tsd_len]
            # Search for the same (or similar) sequence further downstream
            for j in range(i + tsd_len, n - tsd_len + 1):
                target = seq[j : j + tsd_len]
                mismatches = sum(a != b for a, b in zip(candidate, target))
                if mismatches <= max_mismatches:
                    results.append({
                        "tsd_sequence": candidate,
                        "start1":       i,
                        "end1":         i + tsd_len,
                        "start2":       j,
                        "end2":         j + tsd_len,
                        "length":       tsd_len,
                        "mismatches":   mismatches,
                        "gap_between":  j - (i + tsd_len),
                    })
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Detect TSDs flanking candidate Starship insertions."
    )
    parser.add_argument("--fasta",        required=True, help="Genome FASTA file")
    parser.add_argument("--hmmer-hits",   required=True, help="HMMER3 tblout file")
    parser.add_argument("--window",       type=int, default=500,
                        help="bp to extract around each Y-Rec hit boundary (default: 500)")
    parser.add_argument("--min-tsd",      type=int, default=4,
                        help="Minimum TSD length in bp (default: 4)")
    parser.add_argument("--max-tsd",      type=int, default=6,
                        help="Maximum TSD length in bp (default: 6)")
    parser.add_argument("--max-mismatches", type=int, default=0,
                        help="Maximum mismatches allowed in TSD (default: 0 = perfect)")
    parser.add_argument("--output",       required=True, help="Output TSV file")
    args = parser.parse_args()

    # ── Load genome sequences into a dict ────────────────────────────────────
    print(f"Loading genome: {args.fasta}", file=sys.stderr)
    genome: dict[str, SeqRecord] = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    print(f"  Loaded {len(genome)} sequences", file=sys.stderr)

    # ── Parse HMMER hits ──────────────────────────────────────────────────────
    print(f"Parsing HMMER hits: {args.hmmer_hits}", file=sys.stderr)
    hits = parse_hmmer_tblout(args.hmmer_hits)
    print(f"  Found {len(hits)} Y-Rec candidates", file=sys.stderr)

    if not hits:
        print("WARNING: No Y-Rec hits found. Writing empty output.", file=sys.stderr)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as out:
            out.write("contig\thit_name\te_value\tscore\twindow_start\twindow_end\t"
                      "tsd_sequence\ttsd_start1\ttsd_end1\ttsd_start2\ttsd_end2\t"
                      "tsd_length\tmismatches\tgap_between\n")
        return

    # ── For each hit, extract flanking window and search for TSDs ─────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    results_rows = []

    for hit in hits:
        contig_name = hit["target_name"]
        # Six-frame translation IDs have suffix like _frame=+1_start=1234
        # Strip suffixes to get the original contig name
        contig_base = contig_name.split()[0]
        # Find matching contig in genome (prefix match)
        matching_contig = None
        for cid in genome:
            if cid.startswith(contig_base) or contig_base.startswith(cid):
                matching_contig = cid
                break

        if matching_contig is None:
            print(f"  WARNING: Contig '{contig_base}' not found in FASTA. Skipping.",
                  file=sys.stderr)
            continue

        contig_seq = str(genome[matching_contig].seq)
        contig_len = len(contig_seq)

        # Extract a window from each end of the contig (crude boundary estimate)
        # A more precise approach would use the TE library coordinates from EDTA
        for boundary_pos in [args.window, contig_len - args.window]:
            w_start = max(0, boundary_pos - args.window)
            w_end   = min(contig_len, boundary_pos + args.window)
            window_seq = contig_seq[w_start:w_end]

            tsds = find_direct_repeats(
                window_seq,
                min_len=args.min_tsd,
                max_len=args.max_tsd,
                max_mismatches=args.max_mismatches,
            )

            for tsd in tsds:
                results_rows.append({
                    "contig":       matching_contig,
                    "hit_name":     hit["target_name"],
                    "e_value":      hit["e_value"],
                    "score":        hit["score"],
                    "window_start": w_start,
                    "window_end":   w_end,
                    **tsd,
                })

    # ── Write TSV output ──────────────────────────────────────────────────────
    columns = ["contig", "hit_name", "e_value", "score", "window_start", "window_end",
               "tsd_sequence", "start1", "end1", "start2", "end2",
               "length", "mismatches", "gap_between"]

    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in results_rows:
            out.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")

    print(f"TSD search complete. {len(results_rows)} candidate TSDs written to {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
