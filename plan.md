# Pipeline Plan — ToxTA Starship Diversity Analysis

> **Pipeline**: Snakemake workflow for ToxA horizontal gene transfer analysis across 3 fungal genomes  
> **Samples**: *P. nodorum* SN15, *P. tritici-repentis* Pt-1C-BFP, *B. sorokiniana* ND90Pr  
> **Date**: 2026-04-23

---

## 1. What the Pipeline Does (Step-by-Step)

The pipeline has **6 phases** across **16 Snakemake rules**, processing 3 fungal genome assemblies (~35–38 Mb each).

### Phase 1 — Data Acquisition (2 rules)

| Rule | What It Does |
|---|---|
| `download_genome` | Downloads genome FASTA files from NCBI FTP (via `wget` + `gunzip`) using accession + assembly name from `config.yaml`. One download per sample. |
| `assembly_stats` | Runs `assembly-stats` on each genome to produce N50, total size, contig count, and largest scaffold. Quick sanity check. |

### Phase 2a — Quality Control (2 rules)

| Rule | What It Does |
|---|---|
| `run_busco` | Runs BUSCO v5 against the `fungi_odb10` lineage (758 orthologs) to assess genome completeness. Uses Augustus with `botrytis_cinerea` model. 16 threads per sample. |
| `busco_summary_plot` | Collects all BUSCO short summaries and generates a multi-sample stacked bar chart (Complete / Duplicated / Fragmented / Missing). |

### Phase 2b — Repeat Annotation (2 rules)

| Rule | What It Does |
|---|---|
| `run_edta` | Runs EDTA v2 (Extensive de novo TE Annotator) with `--sensitive 1` to discover and annotate all transposable elements. Integrates LTRharvest, LTR_Finder, TIR-Learner, HelitronScanner, and RepeatModeler. Produces a curated TE library, genome-wide GFF3 annotation, and summary table. **This is the slowest step by far.** |
| `parse_edta_summary` | Parses the EDTA `.sum` file into a clean TSV table with columns: Sample, TE_Superfamily, Total_bp, Genome_Fraction_pct. |

### Phase 2c — ToxA BLAST Mining (4 rules)

| Rule | What It Does |
|---|---|
| `make_blast_db` | Builds a BLAST nucleotide database (`makeblastdb`) from each genome assembly. |
| `tblastn_toxa` | Searches for the ToxA protein (UniProt P0C1P7, 118 aa) in all genomes using `tblastn` (protein vs 6-frame translated nucleotide). Post-filters hits at ≥30% identity. |
| `extract_toxa_flanking_region` | Extracts ±50 kb flanking the best ToxA BLAST hit using `seqkit subseq`. This 100 kb window captures the ToxTA transposon and surrounding Starship cargo. |
| `annotate_locus_with_prodigal` | Converts the extracted locus FASTA to GenBank format (via `fasta_to_gbk.py`) for clinker input. Despite the rule name, it uses a custom script, not Prodigal. |

### Phase 3 — Starship Detection (3 rules)

| Rule | What It Does |
|---|---|
| `translate_genome` | Six-frame translates each genome to protein FASTA using `seqkit translate` (min ORF length: 50 aa). |
| `hmmsearch_yrec` | Searches for the Starship captain protein (Y-Recombinase / DUF3435, Pfam PF13408) using HMMER3 `hmmsearch` at E-value ≤ 1e-10. |
| `find_tsds` | Detects Target Site Duplications (4–6 bp direct repeats) flanking candidate Starship insertions using a custom Python script. |

### Phase 3b — Starship Comparison (1 rule)

| Rule | What It Does |
|---|---|
| `compare_starship_structure` | Aggregates Y-Rec hits, TSDs, and BLAST results across all samples into a comparison table showing Starship presence/absence, copy number, TSD sequences, and co-localisation with ToxA. |

### Phase 4 — Visualization & Reporting (4 rules)

| Rule | What It Does |
|---|---|
| `run_clinker` | Generates an interactive HTML synteny plot comparing gene order around ToxA across all 3 species using clinker v0.0.28. |
| `plot_toxa_locus` | Creates a publication-quality PNG diagram of the nested transposon architecture (Starship → ToxTA → ToxA) per sample using matplotlib/dna-features-viewer. |
| `busco_summary_plot` | Multi-sample BUSCO bar chart (described above). |
| `aggregate_report` | Merges assembly stats, BUSCO scores, repeat composition, BLAST hits, and Y-Rec hits into a single summary TSV report. |

---

## 2. Estimated Runtime

> **System assumed**: 8-core workstation, 32 GB RAM, SSD storage, stable internet

| Phase | Rule(s) | Per Sample | × 3 Samples | Parallelisable? | Notes |
|---|---|---|---|---|---|
| **1. Download** | `download_genome` | ~40 sec | **~2 min** | ✅ Yes (3 parallel) | Depends on NCBI bandwidth |
| **1. Stats** | `assembly_stats` | ~5 sec | **< 1 min** | ✅ Yes | Trivial |
| **2a. BUSCO** | `run_busco` | ~12 min | **~36 min** | ⚠️ Sequential (16 threads each saturate CPU) | Can run 1 at a time with `heavy: 16` threads |
| **2a. BUSCO plot** | `busco_summary_plot` | — | **< 1 min** | N/A | Runs once after all BUSCO |
| **2b. EDTA** | `run_edta` | **1.5–4 hours** | **4–12 hours** | ⚠️ Sequential (8 threads each) | **DOMINATES TOTAL RUNTIME** |
| **2b. Parse EDTA** | `parse_edta_summary` | ~5 sec | **< 1 min** | ✅ Yes | Trivial |
| **2c. BLAST DB** | `make_blast_db` | ~10 sec | **< 1 min** | ✅ Yes | |
| **2c. tblastn** | `tblastn_toxa` | ~1 min | **~3 min** | ✅ Yes | 8 threads per sample |
| **2c. Extract locus** | `extract_toxa_flanking_region` | ~10 sec | **< 1 min** | ✅ Yes | |
| **2c. Annotate locus** | `annotate_locus_with_prodigal` | ~30 sec | **~1.5 min** | ✅ Yes | |
| **3. Translate** | `translate_genome` | ~20 sec | **~1 min** | ✅ Yes | |
| **3. HMM search** | `hmmsearch_yrec` | ~15 sec | **< 1 min** | ✅ Yes | |
| **3. TSD detection** | `find_tsds` | ~30 sec | **~1.5 min** | ✅ Yes | |
| **3b. Compare** | `compare_starship_structure` | — | **< 1 min** | N/A | Runs once |
| **4. Clinker** | `run_clinker` | — | **~5 min** | N/A | Runs once |
| **4. Locus plots** | `plot_toxa_locus` | ~30 sec | **~1.5 min** | ✅ Yes | |
| **4. Report** | `aggregate_report` | — | **< 1 min** | N/A | Runs once |

### Total Estimated Runtime

| Scenario | Time |
|---|---|
| **Best case** (fast CPU, `sensitive: 0`) | **~3 hours** |
| **Typical** (8 cores, `sensitive: 1`) | **~5–8 hours** |
| **Worst case** (slow disk, limited RAM, `sensitive: 1`) | **~13 hours** |

> **NOTE**: EDTA accounts for ~85–90% of total runtime. Everything else combined takes under 1 hour. Setting `edta.sensitive: 0` in `config.yaml` can cut EDTA time roughly in half.

---

## 3. Current Completion Status

| Step | Pnodorum_SN15 | Ptritici_repentis_Pt1CBF | Bsorokiniana_ND90Pr |
|---|---|---|---|
| `download_genome` | ✅ Done | ✅ Done | ✅ Done |
| `assembly_stats` | ✅ Done | ✅ Done | ✅ Done |
| `run_busco` | ✅ Done (96.7%) | ❌ **Not run** | ❌ **Not run** |
| `busco_summary_plot` | ❌ Blocked (needs all 3 BUSCO) | — | — |
| `run_edta` | ❌ **Not run** | 🔶 **Partial** (LTR step done, rest incomplete) | ❌ **Not run** |
| `parse_edta_summary` | ❌ Blocked | ❌ Blocked | ❌ Blocked |
| `make_blast_db` | ✅ Done | ✅ Done | ✅ Done |
| `tblastn_toxa` | ✅ Done | ✅ Done | ✅ Done |
| `extract_toxa_flanking_region` | ✅ Done | ✅ Done | ✅ Done |
| `annotate_locus_with_prodigal` | ✅ Done | ✅ Done | ✅ Done |
| `translate_genome` | ✅ Done | ✅ Done | ✅ Done |
| `hmmsearch_yrec` | ✅ Done | ✅ Done | ✅ Done |
| `find_tsds` | ✅ Done | ✅ Done | ✅ Done |
| `run_clinker` | ✅ Done | — | — |
| `plot_toxa_locus` | ❌ Blocked (needs EDTA GFF3) | ❌ Blocked | ❌ Blocked |
| `aggregate_report` | ❌ Blocked (needs BUSCO + EDTA) | — | — |

### Summary
- **Completed**: 13/16 rules for download/BLAST/Starship/clinker
- **Blocking bottleneck**: EDTA (0/3 complete) and BUSCO (1/3 complete)
- **Remaining wall-time**: ~5–12 hours (mostly EDTA)

---

## 4. Suggested Improvements

### 🐛 Bugs / Issues to Fix

#### 4.1 — HMM profile Pfam ID inconsistency
The `starship_analysis.smk` comments say **PF11976** (DUF3435), but the README and config correctly state **PF13408** (Y_recomb_recC). The actual HMM file (`resources/Y_recombinase.hmm`) appears to be PF13408 (correct). **Fix the misleading comments in `starship_analysis.smk` lines 48–53** to avoid confusion.

#### 4.2 — Rule name is misleading: `annotate_locus_with_prodigal`
This rule does **not** use Prodigal — it runs `fasta_to_gbk.py` which does FASTA→GenBank conversion. Rename to `convert_locus_to_gbk` or similar.

#### 4.3 — TSD finder has O(n³) complexity
`find_tsds.py` line 63–94: The `find_direct_repeats()` function uses a triple-nested loop (over TSD length × position i × position j). For a 1000 bp window with TSD lengths 4–6, this is ~1.5 billion comparisons. This is extremely slow and will produce millions of spurious hits (any random 4-mer repeats within 1 kb). Consider:
- Only looking for TSDs at the **actual element boundaries** (not all positions in the window)
- Filtering TSDs by their gap distance (should be ~element size, i.e. >10 kb)
- Using a k-mer hash approach instead of brute-force

#### 4.4 — ToxA position calculation in `plot_toxa_locus.py` is hardcoded
Lines 261–262 subtract a hardcoded `50000` (the flank size) from BLAST coordinates. But the BLAST coordinates are **relative to the full genome**, while the locus FASTA starts at `(hit_start - 50000)`. The locus BED file already has the true start coordinate — read it from the BED file instead of hardcoding.

#### 4.5 — `summarize_results.py` is called by two rules with different arguments
Both `compare_starship_structure` and `aggregate_report` call `summarize_results.py` but pass completely different argument sets (legacy vs. standard). The `compare_starship_structure` rule passes `--yrec-hits`, `--tsds`, `--blast-hits` but the script ignores `--tsds`. This rule may silently produce incomplete output.

---

### ⚡ Performance Improvements

#### 4.6 — EDTA: Consider `sensitive: 0` for initial runs
The current config has `sensitive: 1`. For a first-pass/demo run, setting `sensitive: 0` cuts EDTA time by ~50% with minimal sensitivity loss for well-assembled fungal genomes.

#### 4.7 — BUSCO: Use `threads.heavy` consistently
Config has `threads.heavy: 16` but BUSCO is set to use `threads.heavy` while EDTA uses `edta.threads: 8`. On an 8-core machine, the BUSCO `threads: 16` will oversaturate the CPU. Ensure `threads.heavy` matches your actual core count.

#### 4.8 — EDTA runs in `/tmp` — may fill up
The Linux EDTA rule (lines 63–87 of `repeat_annotation.smk`) copies the genome to `/tmp` and runs EDTA there. Fungal genomes + EDTA intermediates can use 5–10 GB per sample. On machines with a small `/tmp` partition, this will fail silently. Consider using a workspace-local temp directory.

#### 4.9 — Add `resources` directives for memory
EDTA and BUSCO can use 8–16 GB RAM. Add `resources: mem_mb=16000` to these rules to prevent out-of-memory crashes and enable proper SLURM `--mem` requests.

---

### 🔬 Scientific Improvements

#### 4.10 — Add ToxA-negative control genomes
The pipeline currently only analyses 3 ToxA-positive genomes. Adding 1–2 ToxA-negative strains (e.g., *P. tritici-repentis* race 4 DW5) would validate specificity of BLAST/HMM searches and strengthen the synteny comparison.

#### 4.11 — `fasta_to_gbk.py` produces unannotated GenBank files
Clinker's power comes from comparing **annotated** gene clusters. Currently, the GenBank files have no gene/CDS features (no GFF3 is passed to the script, despite it being supported). The clinker HTML will show bare sequence blocks with no gene arrows. **Pass the EDTA GFF3 to `fasta_to_gbk.py`** via the `--gff3` flag, or better yet, run a quick ab initio gene predictor (Augustus/Prodigal) on the locus first.

#### 4.12 — TSD detection strategy is naive
The current approach searches for short direct repeats at contig ends, which is not biologically meaningful for Starship boundary detection. A better approach would be:
1. Use the EDTA annotations + Y-Rec positions to estimate actual Starship boundaries
2. Search for TSDs only at those boundaries (±50 bp, not ±500 bp)
3. Filter TSDs by gap size (should match element size)

#### 4.13 — Missing: Multiple sequence alignment of ToxA
The pipeline mines ToxA across 3 genomes but never aligns the hits. Adding a MAFFT alignment rule would enable:
- Percent identity comparison across species
- Phylogenetic tree of ToxA copies (key figure for the paper)

---

### 🧹 Code Quality

#### 4.14 — Add a `--dry-run` validation step
Add a CI/dry-run check: `snakemake --use-conda -n --quiet` in the README or as a GitHub Action, so changes to the Snakefile are validated before commits.

#### 4.15 — Pin `clinker` version carefully
`clinker==0.0.28` is very old (2021). The latest version may have breaking API changes. Verify this version still works with the current Python 3.10 + BioPython 1.84 environment.

#### 4.16 — The `busco_summary_plot` rule copies files unnecessarily
It copies all BUSCO summaries into a `all_summaries/` subdirectory, runs the BUSCO plotting script, then moves the output. This leaves behind the temporary directory. Add a cleanup step.

#### 4.17 — Add `--forceall` or `--rerun-incomplete` guidance
The README's "Quick Start" doesn't mention how to resume after EDTA crashes (common). Add guidance for `snakemake --rerun-incomplete` or `--unlock`.

---

## 5. Recommended Run Order

If you want to re-run the pipeline from the current state:

```bash
# 1. Complete the two blocking bottlenecks (run in this order):
snakemake --use-conda --cores 8 results/busco/Ptritici_repentis_Pt1CBF/short_summary.specific.fungi_odb10.Ptritici_repentis_Pt1CBF.txt
snakemake --use-conda --cores 8 results/busco/Bsorokiniana_ND90Pr/short_summary.specific.fungi_odb10.Bsorokiniana_ND90Pr.txt

# 2. Run EDTA (the long step — consider running overnight):
snakemake --use-conda --cores 8 --forcerun run_edta

# 3. Once EDTA completes, the rest finishes in minutes:
snakemake --use-conda --cores 8
```
