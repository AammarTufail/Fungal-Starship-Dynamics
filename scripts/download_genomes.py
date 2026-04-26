#!/usr/bin/env python3
"""
Download a single fungal genome from the NCBI FTP, decompress, and place it
at the requested output path.

Called by Snakemake's data_acquisition rules — one invocation per sample:

    python scripts/download_genomes.py \
        --accession  GCA_000146915.2 \
        --assembly   ASM14691v2     \
        --outpath    results/genomes/Pnodorum_SN15/Pnodorum_SN15.fasta \
        --log        results/logs/download/Pnodorum_SN15.log

The NCBI FTP layout used is:
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/<d1>/<d2>/<d3>/
        <accession>_<assembly>/<accession>_<assembly>_genomic.fna.gz
"""

import argparse
import os
import re
import subprocess
import sys


# ─── FTP URL construction ─────────────────────────────────────────────────────

def build_ftp_url(accession: str, assembly_name: str) -> str:
    """Return the HTTPS FTP URL for the *_genomic.fna.gz file."""
    num_part = accession.replace("GCA_", "").split(".")[0].zfill(9)
    d1, d2, d3 = num_part[0:3], num_part[3:6], num_part[6:9]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{d1}/{d2}/{d3}"
    dirname  = f"{accession}_{assembly_name}"
    filename = f"{accession}_{assembly_name}_genomic.fna.gz"
    return f"{base}/{dirname}/{filename}"


# ─── Download + decompress ────────────────────────────────────────────────────

def download_genome(accession: str, assembly_name: str, outpath: str,
                    log_path: str) -> None:
    """Download *_genomic.fna.gz from NCBI FTP and decompress to outpath."""
    url     = build_ftp_url(accession, assembly_name)
    gz_path = outpath + ".gz"

    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    with open(log_path, "w") as log:
        log.write(f"Accession  : {accession}\n")
        log.write(f"Assembly   : {assembly_name}\n")
        log.write(f"URL        : {url}\n")
        log.write(f"Output     : {outpath}\n\n")
        log.flush()

        if os.path.isfile(outpath):
            msg = f"Output already exists, skipping download: {outpath}\n"
            log.write(msg)
            print(msg, end="")
            return

        # ── wget ────────────────────────────────────────────────────────────
        print(f"  Downloading {accession} ...")
        ret = subprocess.run(
            ["wget", "-q", "--show-progress", "-O", gz_path, url],
            timeout=600,
            stderr=subprocess.STDOUT,
            stdout=log,
        )
        if ret.returncode != 0:
            if os.path.isfile(gz_path):
                os.remove(gz_path)
            msg = f"ERROR: wget failed for {url}\n"
            log.write(msg)
            sys.exit(msg)

        # ── gunzip ──────────────────────────────────────────────────────────
        print(f"  Decompressing ...")
        subprocess.run(["gunzip", "-f", gz_path], check=True, stderr=log,
                       stdout=log)

        # gz is named outpath+".gz" → gunzip produces outpath automatically
        if not os.path.isfile(outpath):
            sys.exit(
                f"ERROR: Expected decompressed file not found: {outpath}\n"
            )

        size_mb = os.path.getsize(outpath) / 1e6
        msg = f"SUCCESS: {outpath}  ({size_mb:.1f} MB)\n"
        log.write(msg)
        print(f"  {msg}", end="")


# ─── CLI ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download one NCBI genome via FTP and decompress it."
    )
    p.add_argument("--accession",  required=True,
                   help="NCBI GenBank accession, e.g. GCA_000146915.2")
    p.add_argument("--assembly",   required=True,
                   help="NCBI assembly name, e.g. ASM14691v2")
    p.add_argument("--outpath",    required=True,
                   help="Destination FASTA path (without .gz suffix)")
    p.add_argument("--log",        required=True,
                   help="Path to write log output")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(os.path.dirname(args.log), exist_ok=True)
    download_genome(
        accession=args.accession,
        assembly_name=args.assembly,
        outpath=args.outpath,
        log_path=args.log,
    )


if __name__ == "__main__":
    main()
