#!/usr/bin/env python3
"""
plot_toxa_locus.py — Linear locus diagram showing the nested transposon architecture.

Visualises the key nested structure:
    [Starship outer boundary]
        └── [ToxTA passenger transposon]
                └── [ToxA effector gene]

Reads:
  - Extracted locus FASTA (from blast_mining.smk extract_toxa_flanking_region rule)
  - EDTA GFF3 annotation (repeat elements in the locus)
  - BLAST hit TSV (ToxA coordinates in the locus)

Outputs a publication-quality PNG figure using dna-features-viewer.

Usage:
    python plot_toxa_locus.py \\
        --fasta results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_locus.fasta \\
        --te-gff3 results/edta/Pnodorum_SN15/Pnodorum_SN15.fasta.mod.EDTA.anno/Pnodorum_SN15.fasta.mod.TE.gff3 \\
        --blast-hits results/blast/Pnodorum_SN15/Pnodorum_SN15_toxa_hits.tsv \\
        --sample Pnodorum_SN15 \\
        --output results/figures/Pnodorum_SN15_toxa_locus.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for HPC/headless runs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow

from Bio import SeqIO


# ── Colour palette (colour-blind friendly) ───────────────────────────────────
COLOURS = {
    "ToxA":        "#E63946",   # red  — the focal gene
    "ToxTA":       "#457B9D",   # blue — passenger transposon
    "Starship":    "#2D6A4F",   # green — Starship element
    "LTR":         "#A8DADC",   # light blue — other LTR retrotransposons
    "DNA_TE":      "#F4A261",   # orange — DNA transposons
    "Unknown_TE":  "#CCC5B9",   # grey — unclassified TEs
    "Gene":        "#8338EC",   # purple — other genes
}


def parse_blast_hits(tsv_path: str) -> list[dict]:
    """Parse BLAST tabular output (outfmt 6) into a list of dicts."""
    hits = []
    # Expected columns from config outfmt:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    col_names = ["qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < len(col_names):
                continue
            hit = dict(zip(col_names, fields))
            hit["sstart"] = int(hit["sstart"])
            hit["send"]   = int(hit["send"])
            hit["pident"] = float(hit["pident"])
            hit["bitscore"] = float(hit["bitscore"])
            hits.append(hit)
    return sorted(hits, key=lambda h: h["bitscore"], reverse=True)


def parse_gff3_for_contig(gff3_path: str, contig: str,
                           locus_start: int, locus_end: int) -> list[dict]:
    """
    Extract TE annotations from a GFF3 file that fall within the locus window.
    Returns coordinates relative to locus_start (0-based for plotting).
    """
    features = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            gff_contig, source, ftype, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            if gff_contig != contig:
                continue
            if end < locus_start or start > locus_end:
                continue
            # Clip to locus window
            rel_start = max(start, locus_start) - locus_start
            rel_end   = min(end,   locus_end)   - locus_start
            # Classify TE type from the feature type or attributes
            te_class = classify_te(ftype, attrs)
            features.append({
                "start":   rel_start,
                "end":     rel_end,
                "strand":  strand,
                "te_class": te_class,
                "attrs":   attrs,
            })
    return features


def classify_te(ftype: str, attrs: str) -> str:
    """Map EDTA GFF3 feature types to display categories."""
    a = attrs.lower()
    f = ftype.lower()
    if "toxa" in a:
        return "ToxA"
    if "toxta" in a or "tc1" in a or "mariner" in a:
        return "ToxTA"
    if "starship" in a or "duf3435" in a:
        return "Starship"
    if "ltr" in f or "ltr" in a:
        return "LTR"
    if "dna" in a or "helitron" in a or "tir" in f:
        return "DNA_TE"
    return "Unknown_TE"


def draw_locus(
    locus_length: int,
    te_features: list[dict],
    toxa_start: int,
    toxa_end: int,
    toxa_strand: str,
    sample_name: str,
    output_path: str,
):
    """
    Draw a linear locus diagram with nested transposon features.

    Layout (top to bottom):
      Track 1: Starship element span
      Track 2: ToxTA passenger transposon span
      Track 3: ToxA gene arrow
      Track 4: Other TEs from EDTA
    """
    fig, ax = plt.subplots(figsize=(16, 4))
    ax.set_xlim(0, locus_length)
    ax.set_ylim(0, 5)
    ax.axis("off")

    y_base = 2.0
    track_height = 0.5
    arrow_height = 0.35

    # ── Draw all TE features ──────────────────────────────────────────────────
    for feat in te_features:
        color = COLOURS.get(feat["te_class"], COLOURS["Unknown_TE"])
        y_offset = {
            "Starship":   3.5,
            "ToxTA":      2.8,
            "ToxA":       2.1,
        }.get(feat["te_class"], y_base)

        rect = mpatches.FancyBboxPatch(
            (feat["start"], y_offset - track_height / 2),
            feat["end"] - feat["start"],
            track_height,
            boxstyle="round,pad=0.01",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.75,
        )
        ax.add_patch(rect)

    # ── Draw ToxA gene as an arrow (from BLAST coordinates) ──────────────────
    toxa_color = COLOURS["ToxA"]
    arrow_y = 2.1
    if toxa_strand == "+":
        ax.annotate(
            "", xy=(toxa_end, arrow_y), xytext=(toxa_start, arrow_y),
            arrowprops=dict(arrowstyle="->", color=toxa_color, lw=2.0),
        )
    else:
        ax.annotate(
            "", xy=(toxa_start, arrow_y), xytext=(toxa_end, arrow_y),
            arrowprops=dict(arrowstyle="->", color=toxa_color, lw=2.0),
        )
    ax.text(
        (toxa_start + toxa_end) / 2, arrow_y + 0.25,
        "ToxA", ha="center", va="bottom", fontsize=8, fontweight="bold",
        color=toxa_color,
    )

    # ── Scale bar ─────────────────────────────────────────────────────────────
    scale_len = 10_000  # 10 kb
    ax.plot([locus_length * 0.05, locus_length * 0.05 + scale_len],
            [0.4, 0.4], color="black", lw=2)
    ax.text(locus_length * 0.05 + scale_len / 2, 0.55,
            "10 kb", ha="center", va="bottom", fontsize=8)

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=COLOURS["Starship"],  label="Starship element"),
        mpatches.Patch(color=COLOURS["ToxTA"],     label="ToxTA passenger TE"),
        mpatches.Patch(color=COLOURS["ToxA"],      label="ToxA effector gene"),
        mpatches.Patch(color=COLOURS["LTR"],       label="LTR retrotransposon"),
        mpatches.Patch(color=COLOURS["DNA_TE"],    label="DNA transposon"),
        mpatches.Patch(color=COLOURS["Unknown_TE"],label="Unclassified TE"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=7, framealpha=0.9,
              ncol=3, bbox_to_anchor=(1, 1))

    # ── Title ─────────────────────────────────────────────────────────────────
    ax.set_title(
        f"ToxA locus — {sample_name.replace('_', ' ')}\n"
        f"Nested transposon architecture: ToxA ← ToxTA ← Starship",
        fontsize=10, pad=8,
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Plot ToxA locus architecture.")
    parser.add_argument("--fasta",      required=True, help="Locus FASTA (extracted)")
    parser.add_argument("--te-gff3",    required=True, help="EDTA TE GFF3 annotation")
    parser.add_argument("--blast-hits", required=True, help="BLAST tblastn hits TSV")
    parser.add_argument("--sample",     required=True, help="Sample name for title")
    parser.add_argument("--output",     required=True, help="Output PNG path")
    args = parser.parse_args()

    # ── Load locus sequence ──────────────────────────────────────────────────
    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        print(f"WARNING: Empty FASTA {args.fasta}. Cannot plot.", file=sys.stderr)
        # Write a blank placeholder figure
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, f"No ToxA locus found in {args.sample}",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=150)
        plt.close(fig)
        return

    locus_record = records[0]
    locus_length = len(locus_record.seq)
    locus_contig = locus_record.id

    # ── Parse BLAST hits for ToxA position ──────────────────────────────────
    blast_hits = parse_blast_hits(args.blast_hits)
    if not blast_hits:
        print(f"WARNING: No BLAST hits for {args.sample}.", file=sys.stderr)
        toxa_start, toxa_end, toxa_strand = 0, 0, "+"
    else:
        best = blast_hits[0]
        # The locus was extracted starting at (hit_start - flank). The ToxA
        # position within the locus is approximately at the midpoint.
        toxa_start  = max(0, best["sstart"] - 50000)   # rough offset; refine if needed
        toxa_end    = max(0, best["send"]   - 50000)
        toxa_strand = "+" if best["sstart"] < best["send"] else "-"
        if toxa_start < 0: toxa_start = 0
        if toxa_end   > locus_length: toxa_end = locus_length

    # ── Parse EDTA GFF3 for features in this locus ──────────────────────────
    # The locus BED coordinates are needed; here we approximate using locus_length
    te_features = parse_gff3_for_contig(
        args.te_gff3,
        contig=locus_contig,
        locus_start=0,
        locus_end=locus_length,
    )
    print(f"Loaded {len(te_features)} TE features for plotting.", file=sys.stderr)

    # ── Draw ─────────────────────────────────────────────────────────────────
    draw_locus(
        locus_length=locus_length,
        te_features=te_features,
        toxa_start=toxa_start,
        toxa_end=toxa_end,
        toxa_strand=toxa_strand,
        sample_name=args.sample,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
