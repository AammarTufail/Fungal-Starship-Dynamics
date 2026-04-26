#!/usr/bin/env python3
"""
parse_edta_summary.py — Parse an EDTA summary file into a clean TSV table.

EDTA produces a summary file (<genome>.mod.EDTA.sum) that lists the
per-superfamily repeat content of the genome. This script parses that file
into a normalised TSV with columns:

    Sample | TE_Superfamily | Total_bp | Genome_Fraction_pct

This output is used by summarize_results.py to build the cross-sample
repeatome comparison table and figures.

Usage:
    python parse_edta_summary.py \\
        --summary results/edta/Pnodorum_SN15/Pnodorum_SN15.fasta.mod.EDTA.sum \\
        --sample Pnodorum_SN15 \\
        --output results/edta/Pnodorum_SN15/Pnodorum_SN15_repeat_summary.tsv
"""

import argparse
import re
import sys
from pathlib import Path


def parse_edta_sum(sum_path: str, sample: str) -> list[dict]:
    """
    Parse EDTA .sum file.

    The EDTA summary file format (v2) looks like:
    -----------------------------------------------
    Genome size (bp): 37,284,512
    ...
    Class               Coverage (bp)  Fraction (%)
    -----------------------------------------------
    LTR/Gypsy           3,245,123       8.70
    LTR/Copia           1,234,567       3.31
    TIR/Tc1-Mariner     456,789         1.22
    ...
    Total repeats       9,876,543      26.49
    -----------------------------------------------

    Returns a list of dicts with keys: Sample, TE_Superfamily, Total_bp, Genome_Fraction_pct
    """
    rows = []
    in_table = False

    with open(sum_path) as fh:
        for line in fh:
            line = line.strip()

            # Detect start of the TE table
            if re.match(r"^(Class|TE Class|Category)\s+", line, re.IGNORECASE):
                in_table = True
                continue

            # Skip separator lines
            if re.match(r"^[-=]+$", line):
                continue

            if not in_table:
                continue

            # Stop at blank lines or summary footer
            if not line:
                continue

            # Parse data rows: two formats observed in EDTA output
            # Format A: "LTR/Gypsy   3245123   8.70"
            # Format B: "LTR    Gypsy    3245123   8.70"
            m = re.match(
                r"^([\w/\-]+(?:\s+[\w/\-]+)?)\s+([\d,]+)\s+([\d.]+)",
                line,
            )
            if m:
                te_class   = m.group(1).strip().replace(" ", "/")
                total_bp   = int(m.group(2).replace(",", ""))
                fraction   = float(m.group(3))
                # Skip the "Total" summary line
                if te_class.lower().startswith("total"):
                    continue
                rows.append({
                    "Sample":             sample,
                    "TE_Superfamily":     te_class,
                    "Total_bp":          total_bp,
                    "Genome_Fraction_pct": fraction,
                })

    return rows


def main():
    parser = argparse.ArgumentParser(description="Parse EDTA summary file to TSV.")
    parser.add_argument("--summary", required=True, help="EDTA .sum file")
    parser.add_argument("--sample",  required=True, help="Sample name")
    parser.add_argument("--output",  required=True, help="Output TSV path")
    args = parser.parse_args()

    rows = parse_edta_sum(args.summary, args.sample)

    if not rows:
        print(f"WARNING: No TE rows parsed from {args.summary}. "
              "Check EDTA version and output format.", file=sys.stderr)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write("Sample\tTE_Superfamily\tTotal_bp\tGenome_Fraction_pct\n")
        for row in rows:
            out.write(
                f"{row['Sample']}\t{row['TE_Superfamily']}\t"
                f"{row['Total_bp']}\t{row['Genome_Fraction_pct']}\n"
            )

    print(f"Parsed {len(rows)} TE superfamily entries -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
