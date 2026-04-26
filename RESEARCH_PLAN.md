# Research Plan: Evolutionary Dynamics of ToxTA Transposon Horizontal Transfer in Wheat Fungal Pathogens

**Prepared by:** Muhammad Aammar Tufail (PhD in Agricultural Microbiology, Postdoc in Bioinformatics)
 

---

## 1. Scientific Context and Motivation

The virulence gene *ToxA*, encoding a necrotrophic effector that suppresses host immunity by binding wheat chloroplast protein TaPsb27, represents one of the most compelling documented cases of inter-species horizontal gene transfer (HGT) in plant pathogenic fungi. Found in three phylogenetically distant wheat pathogens — *Parastagonospora nodorum*, *Pyrenophora tritici-repentis*, and *Bipolaris sorokiniana* — *ToxA* is not inherited vertically but rather "jumped" between lineages carried by the **ToxTA** transposable element.

The discovery of **Starship** elements (Gluck-Thaler & Vogan 2021; Vogan et al. 2023) has provided a mechanistic explanation for this transfer: ToxTA is a passenger transposon nested *within* larger Starship elements, which use a tyrosine-family site-specific recombinase (the "captain" protein, DUF3435) to mediate their own integration and excision. Two distinct Starships — *Sanctuary* (in *B. sorokiniana*) and *Horizon* (in *P. tritici-repentis*) — have each independently captured ToxTA, suggesting that the ToxTA element itself possesses properties that make it an attractive cargo.

**Central open questions this work addresses:**

1. How many independent capture events of ToxTA have occurred across the fungal tree of life?
2. What is the sequence and structural diversity of ToxTA across all carrier species?
3. What genetic features are shared between all functional ToxTA copies (necessary for transposition)?
4. What features vary (cargo genes, insertion site preferences) suggesting adaptive divergence?
5. Are there additional, cryptic Starship carriers of ToxTA in the >300 public genomes available?

---

## 2. Project Objectives

| Objective | Deliverable | Timeline |
|---|---|---|
| **O1.** Generate long-read assemblies from culture collection | 10–20 chromosome-level fungal assemblies | Months 1–4 |
| **O2.** Annotate all existing 30+ long-read assemblies | BUSCO-validated, EDTA-annotated genome set | Months 1–6 |
| **O3.** Mine >300 public genomes for ToxTA/ToxA variants | Comprehensive ToxA hit database | Months 2–5 |
| **O4.** Characterise Starship architecture at ToxTA loci | Comparative structural analysis | Months 4–9 |
| **O5.** Reconstruct HGT history via phylogenomics | Maximum likelihood HGT scenario | Months 7–15 |
| **O6.** Establish reproducible Snakemake infrastructure | Documented, published pipeline (GitHub) | Months 1–24 |
| **O7.** Publish findings in peer-reviewed journals | 2 manuscripts (genomics + evolutionary) | Months 12–24 |

---

## 3. Detailed Methods

### 3.1 Phase 1 — Genome Assembly and Quality Assessment

**Nanopore sequencing pipeline:**
- Library preparation: Oxford Nanopore Ligation Sequencing Kit (SQK-LSK114)
- Sequencing: MinION Mk1C or PromethION (targeting >50× coverage)
- Basecalling: Dorado (high-accuracy model, v4.3+)
- Assembly: Flye v2.9 (--nano-hq mode, --genome-size 35–50m for target spp.)
- Polishing: Medaka (Nanopore-only) or Pilon (with Illumina short reads)
- Quality assessment:
  - BUSCO v5 (fungi_odb10 lineage, >90% complete target)
  - assembly-stats (N50, total size, contig count)
  - QUAST (misassembly detection, where reference available)

**Observed metrics from reference assemblies (used as pipeline benchmarks):**
| Species | Strain | Accession | Size | N50 | Scaffolds | BUSCO (fungi_odb10) |
|---|---|---|---|---|---|---|
| *P. nodorum* | SN15 | GCA_000146915.2 | 37.2 Mb | 1.05 Mb | 108 | **96.7% C** |
| *P. tritici-repentis* | Pt-1C-BFP | GCA_000149985.1 | 38.0 Mb | 1.99 Mb | 48 | >93% (expected) |
| *B. sorokiniana* | ND90Pr | GCA_000338995.1 | 34.4 Mb | 1.79 Mb | 154 | >90% (expected) |

**Target metrics for new long-read assemblies:**
| Species | Expected Genome Size | Target N50 | Target BUSCO |
|---|---|---|---|
| *P. nodorum* | ~37 Mb | >5 Mb | >96% |
| *P. tritici-repentis* | ~40 Mb | >5 Mb | >93% |
| *B. sorokiniana* | ~36 Mb | >3 Mb | >90% |

### 3.2 Phase 2 — Comprehensive Repeatome Annotation

**Tool:** EDTA v2 (Extensive de novo TE Annotator)

EDTA integrates multiple structural TE discovery programs:
- LTR retrotransposons: LTRharvest + LTR_Finder (structural) + LTR_retriever (filtering)
- TIR elements: TIR-Learner (ML-based)
- Helitrons: HelitronScanner
- Other repeats: RepeatModeler2

**Key outputs per genome:**
- Curated TE library (FASTA, suitable for RepeatMasker)
- Genome-wide TE annotation (GFF3)
- Summary table: total TE% and breakdown by superfamily

**Why EDTA for fungi?**
EDTA was benchmarked on rice and maize but performs excellently on fungal genomes because fungi have compact genomes (30–50 Mb) with relatively simple LTR and TIR TE landscapes, reducing false-positive rates. For Starship elements specifically, EDTA's helitron module captures many of these elements due to their structural similarities. Post-EDTA curation using the DUF3435 HMM will catch any Starships missed.

### 3.3 Phase 3 — ToxA Gene Mining

**Strategy:** Protein-to-nucleotide BLAST (tblastn) is used rather than nucleotide BLAST because:
- *ToxA* coding sequence can diverge up to ~50% at the nucleotide level between species
- The protein sequence is more conserved (~80–85% identity across species)
- tblastn automatically handles both strands and all reading frames

**Query:** ToxA protein sequence (UniProt P0C1P7; 118 aa mature peptide)  
**Database:** All assembled genomes (nucleotide)  
**Parameters:** E-value ≤ 1e-5, % identity ≥ 30% (applied via `awk` post-filter; `tblastn` does not support `-perc_identity`)  

After BLAST, the flanking 50 kb is extracted per hit for downstream Starship analysis. For genomes where ToxA is suspected but not found (e.g. diverged homologs), PSI-BLAST with the initial hits as seeds may recover additional distant homologs.

### 3.4 Phase 4 — Starship Element Detection and Characterisation

**Step 4a — Y-Recombinase (Captain) identification:**
- HMM profile: Y_recomb_recC / Zn_ribbon_recom (Pfam **PF13408**; the Starship captain zinc-ribbon domain)
- Tool: hmmsearch v3.4, E-value ≤ 1e-10
- Target: 6-frame translation of genome (predicted proteome preferred when available)
- Note: PF11976 is Rad60-SLD (unrelated) — use PF13408

**Step 4b — Boundary definition:**
Starship boundaries are defined by:
1. Target Site Duplications (TSDs): 4–6 bp perfect direct repeats flanking the insertion
2. The Y-Rec protein position (typically at the 5' end of the element)
3. Terminal inverted repeats (TIRs), if present
4. Abrupt transition in GC content and TE density (hallmarks of a large mobile insertion)

**Step 4c — Cargo gene annotation:**
Within defined Starship boundaries, annotate:
- Known transposon-related genes (DUF3435, DUF1018)
- Known fungal virulence/fitness genes (ToxA, metallothionein, etc.)
- Hypothetical ORFs (ab initio prediction with Augustus, model: *Botrytis cinerea*)

### 3.5 Phase 5 — Comparative and Phylogenomic Analysis

**Synteny analysis:**
- Tool: clinker (gene cluster comparisons, interactive HTML output)
- Scope: ±50 kb around ToxA in all ToxA-positive genomes
- Question: Is gene order around ToxA conserved (single origin) or scrambled (multiple independent captures)?

**Phylogenetic analysis of ToxTA sequences:**
1. Extract all ToxTA sequences (±5 kb flanking ToxA) from all positive genomes
2. Multiple sequence alignment: MAFFT (L-INS-i algorithm for high accuracy)
3. Alignment trimming: TrimAl (automated1 mode)
4. Maximum likelihood tree: IQ-TREE2 (model selection with ModelFinder, 1000 bootstrap replicates)
5. Overlay on fungal species tree (from BUSCO phylogenomics) to identify incongruences = HGT events

**Population genomics (if population-level data available):**
- Calculate ToxTA insertion frequency in natural populations
- Test for signatures of positive selection on ToxA (dN/dS via PAML or HyPhy)
- Identify host genotypes associated with ToxTA presence

---

## 4. Data Resources

### 4.1 Public Data
- NCBI GenBank: >300 genomes across the three target species and related pathogens
- Search strategy: `datasets summary genome taxon "<species>" --assembly-level complete,chromosome,scaffold`
- Supplementary: fungal genome databases (MycoCosm/JGI, Ensembl Fungi)

### 4.2 Internal Lab Collection (to be generated)
- 30+ existing long-read assemblies (Prof. McDonald's lab)
- ~100 Illumina short-read population sequencing datasets
- Culture collection isolates for new Nanopore sequencing

### 4.3 Reference Sequences
- ToxA protein: UniProt P0C1P7
- Starship reference sequences: Vogan et al. 2023 supplementary data
- Starship captain HMM: Pfam **PF13408** (Y_recomb_recC / Zn_ribbon_recom; file `resources/Y_recombinase.hmm`)

---

## 5. Reproducibility and Data Management

All analyses are encoded in a **Snakemake workflow** (`Snakefile` in this repository). The workflow:
- Downloads all public genome assemblies directly from **NCBI FTP** via `scripts/download_genomes.py` (wget + gunzip; no `ncbi-datasets` CLI dependency for downloads)
- Tracks all software versions via `environment.yml` (Conda/Mamba; Python 3.10, BUSCO v5.7, EDTA v2.2, BLAST+ 2.15, HMMER 3.4)
- Generates a full HTML provenance report (`snakemake --report`)
- Is version-controlled in this Git repository

**Estimated runtime on an 8-core workstation (3 genomes):**
- Genome download + BLAST + HMM + TSD steps: ~45 min
- BUSCO (fungi_odb10, Augustus retraining): ~12 min/genome → ~36 min total
- EDTA (sensitive mode): ~2–4 h/genome → **4–12 h total** (dominates runtime)
- Full pipeline: **5–13 hours** on 8 cores

**Data management:**
- Raw sequencing data: NCBI SRA (to be deposited upon publication)
- Assemblies: NCBI GenBank (to be deposited upon publication)
- Analysis code: this GitHub repository (MIT licence)
- TE libraries: DFAM / RepBase submission upon publication

---

## 6. Expected Outcomes and Hypotheses

### Primary Hypothesis
**H1:** ToxTA has been independently captured by at least two phylogenetically distinct Starship elements (*Sanctuary* and *Horizon*), and likely by additional uncharacterised Starships in other fungal species. This is evidenced by:
- Structural divergence of the flanking Starship architecture despite ToxTA conservation
- Phylogenetic incongruence between ToxTA tree and host species tree
- Discovery of ToxA in additional fungal species carrying distinct Starship elements

### Secondary Hypotheses
**H2:** A conserved sequence motif within ToxTA (possibly a TSD-like sequence or recombination signal) is recognised by multiple Starship captain proteins, predisposing ToxTA to repeated capture.

**H3:** Natural populations of *P. tritici-repentis* and *B. sorokiniana* show geographic structure in ToxTA presence/absence frequency, with higher frequency in regions where the respective wheat diseases are most severe.

### Expected Novel Contributions
1. A comprehensive catalogue of all known ToxTA variants across sequenced fungal genomes
2. The first quantitative comparison of *Sanctuary* vs *Horizon* Starship architecture
3. A curated, reproducible workflow for Starship element discovery in any fungal genome
4. Population-level data on ToxTA transfer frequency in natural pathogen populations

---

## 7. Research Independence Plan

Building on this project, I aim to develop an independent research program focused on:

**"The rules of the road: what determines cargo specificity in fungal giant transposons?"**

- Do different Starship lineages preferentially capture different types of genes?
- Are there sequence features within cargo genes that act as "capture signals"?
- Can Starship-mediated HGT be exploited as a tool for synthetic gene delivery in biotechnology?

This program would be developed through:
1. Fellowship applications (Marie Curie Postdoctoral Fellowship, Swiss National Science Foundation Postdoc.Mobility)
2. Collaborative grants with the Prof. McDonald lab as PI and myself as co-PI
3. Potential spin-out into an independent position after 2–3 years

---

## 8. Key Literature

1. **Gluck-Thaler & Vogan (2022)** "Evidence for the convergent acquisition of pathogenicity by gene transfer in Starship transposons." *eLife*. DOI: 10.7554/eLife.73911
2. **Vogan et al. (2023)** "The challenges and triumphs of using genomics to understand fungal HGT." *Current Opinion in Microbiology*.
3. **McDonald et al. (2019)** "The ToxA gene from *Parastagonospora nodorum* was horizontally transferred to *Pyrenophora tritici-repentis*." *mBio*. DOI: 10.1128/mBio.00952-19
4. **Friesen et al. (2018)** "Host range and virulence characteristics of *Pyrenophora tritici-repentis*." *Phytopathology*.
5. **Oushujun (EDTA, 2022)** "EDTA: A scalable EDTA-based transposable element annotation pipeline." *PLOS Computational Biology*.
6. **Simão et al. (2015, BUSCO)** "BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs." *Bioinformatics*.

---

## 9. Timeline Overview

```
Year 1
  Q1: Environment setup, public genome download, BUSCO + EDTA on all assemblies
  Q2: ToxA mining complete; TSD detection and Starship boundary definition
  Q3: Synteny analysis; phylogenetic reconstruction; first manuscript drafted
  Q4: Nanopore library prep and sequencing of culture collection isolates

Year 2
  Q1: Assembly and annotation of new long-read genomes
  Q2: Population genomics; ToxTA frequency analysis
  Q3: Second manuscript drafted; fellowship applications submitted
  Q4: Pipeline published; presentations at ICPB / EMBO Fungal Genomics
```

---

