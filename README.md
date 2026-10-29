# MethylBench

> Reproducible analysis code for a systematic benchmark of six DNA methylation profiling technologies across diverse sequencing platforms.

---

## Overview

MethylBench systematically compares six widely used DNA methylation profiling technologies:

| Technology | Type | CpG Coverage |
|---|---|---|
| Illumina EPIC array | Array-based | ~850k CpGs |
| TWIST Methylation Panel | Targeted short-read | 2–4 million CpGs |
| Whole-Genome Enzymatic Conversion (WGEC) | Genome-wide short-read | 28–30 million CpGs |
| Reduced Representation Bisulfite Sequencing (RRBS) | Enrichment-based short-read | 1–4 million CpGs |
| Oxford Nanopore Technologies (ONT) | Long-read | 28–30 million CpGs |
| Pacific Biosciences (PacBio) | Long-read | 28–30 million CpGs |

Analyses were performed on matched blood and fibroblast samples from 5 individuals and two Genome in a Bottle (GIAB) reference samples (HG001, HG002).

---

## Repository Structure

```
MethylBench/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── scripts/
│   ├── bash/
│   │   ├── 01_modkit_pileup.sh                  # modkit: extract methylation calls from aligned ONT .bam files
│   │   ├── 02_toulligqc.sh                      # ToulligQC: QC of ONT data, intermediary report files
│   │   ├── 03_pbcpgtools.sh                     # pb-CpG-tools: extract methylation calls from PacBio alignments
│   │   ├── 04_methylseq.sh                      # nf-core/methylseq: short-read processing (RRBS, WGEC, TWIST)
│   │   └── run_methylbench.sh                   # runs the complete R analysis (scripts 08–18) in one go
│   ├── python/
│   │   └── 05_parse_toulligqc.py                # summarize ToulligQC .data files into one QC table
│   └── R/
│       ├── 06_visualize_toulligqc_summary.R     # ONT QC report visualization
│       ├── 07_generate_cpg_stats.R              # cpg_stats_methylbench.tab (input for 08)
│       ├── 08_qc_visualization.R                # Figure 3, Suppl. Figure 1, Suppl. Table S6
│       ├── 09_correlation_analysis.R            # Figure 4, Suppl. Figures 2 and 3, Suppl. Table S7
│       ├── 10_density_plots.R                   # Figure 5, Suppl. Figure 5
│       ├── 11_pca.R                             # Figure 6, Suppl. Figure 6
│       ├── 12_differential_methylation.R        # Figure 7, Suppl. Figure 7 (descriptive, no testing)
│       ├── 13_DMR_DSS_analysis.R                # Figures 8, 9A/B/D, Suppl. Figures 8 and 10, Suppl. Table S2
│       ├── 14_downsampling_sensitivity.R        # Suppl. Figure 9, Suppl. Table S5
│       ├── 15_annotation_enrichment_background.R# Figure 9C, Suppl. Figure 11, Suppl. Tables S3/S4
│       ├── 16_dmr_probe_density.R               # Suppl. Figure 12, Suppl. Table S11
│       ├── 17_imprinted_regions.R               # Suppl. Figure 4, Suppl. Table S12
│       ├── 18_coverage_by_context.R             # Suppl. Figures 1E,F and 12E, Suppl. Tables S6B and S11F
│       └── utils/
│           └── helpers.R                        # shared functions, colours, plotting helpers
│
├── data/
│   └── imprinted_gDMRs_hg38.bed                 # 67 imprinted germline DMRs (hg38), input for script 17
│
└── envs/
    ├── methylbench.yml                          # main environment (portable, pinned versions): scripts 05–18
    ├── methylbench_linux-64_explicit.txt        # same environment, bit-identical explicit export (Linux-64)
    ├── ont.yml                                  # ONT tools: modkit, ToulligQC
    └── pacbio.yml                               # PacBio tools: pb-CpG-tools
```

---

## Requirements

- [Conda](https://docs.conda.io/en/latest/) ≥ 23.x (or mamba)
- [Nextflow](https://www.nextflow.io/) and [Singularity](https://docs.sylabs.io/)/Apptainer – only for short-read preprocessing (`04_methylseq.sh`); system installation
- Internet access on the first run of scripts 13–16 and 18: annotatr downloads the hg38 CpG-island track, and `DMRcate::extractRanges()` downloads its gene annotation via ExperimentHub. Both are cached afterwards.

All conda dependencies are defined in `envs/`. The R analysis uses R 4.4.2 / Bioconductor 3.20 (DSS 2.54.0, DMRcate 3.2.1), as reported in the manuscript.

---

## Installation

```bash
git clone https://github.com/epigenetics-sb/MethylBench
cd MethylBench

# main analysis environment (portable, versions pinned)
conda env create -f envs/methylbench.yml   # environment name: methylbench

# alternatively: bit-identical environment on Linux-64
# conda create --name methylbench --file envs/methylbench_linux-64_explicit.txt

conda env create -f envs/ont.yml           # environment name: methylbench-ont
conda env create -f envs/pacbio.yml        # environment name: methylbench-pacbio
```

Nextflow and Singularity/Apptainer, needed only for `04_methylseq.sh`, are not part of the conda environments and must be available on the system (e.g. as HPC modules).

---

## Preprocessing (raw data → per-CpG methylation calls)

The bash scripts (01–04) and `06_visualize_toulligqc_summary.R` take positional arguments, which are documented in the header of each script. `05_parse_toulligqc.py` and `07_generate_cpg_stats.R` support `--help`.

```bash
conda activate methylbench-ont
bash scripts/bash/01_modkit_pileup.sh  <arguments, see script header>
bash scripts/bash/02_toulligqc.sh      <arguments, see script header>
conda deactivate

conda activate methylbench-pacbio
bash scripts/bash/03_pbcpgtools.sh     <arguments, see script header>
conda deactivate

conda activate methylbench
bash scripts/bash/04_methylseq.sh      <arguments, see script header>
python3 scripts/python/05_parse_toulligqc.py --help
Rscript scripts/R/06_visualize_toulligqc_summary.R <arguments, see script header>
Rscript scripts/R/07_generate_cpg_stats.R --help
```

The per-CpG calls of all platforms are then merged into the matrices described below.

---

## Input data

The analysis (scripts 08–16 and 18) needs **five methylation matrices** in one directory (`$DATA`) and **two statistics tables** in another (`$STATS`). File names must match exactly. Script 17 additionally needs the per-platform methylation calls (see [Inputs for script 17](#inputs-for-script-17-imprinted-gdmrs)).

### Methylation matrices (`$DATA`)

| File | Samples | Platforms | Used by |
|---|---|---|---|
| `Blood_without_EPIC.csv` | Blood1–Blood5 | ONT, WGEC, TWIST, RRBS | 09, 10, 11, 12 |
| `Fibro_without_EPIC.csv` | Fibro1–Fibro5 | ONT, WGEC, TWIST, RRBS | 09, 10, 11, 12 |
| `GIAB_without_EPIC.csv` | GIAB1, GIAB2 | ONT, WGEC, TWIST, RRBS, PacBio | 09, 10, 11 |
| `ALL_without_EPIC.csv` | all 12 samples | all sequencing platforms (union of the three files above) | 08 (`--all_path`), 13 (`--seq_path`, Tier 1), 18 (`--all_path`) |
| `ALL.csv` | all 12 samples | `ALL_without_EPIC.csv` **plus** the EPIC columns | 10 (`--epic_path`), 11, 12, 13 (`--all_path`, Tier 2), 16 (`--all_path`, optional) |

**Format** (all five files):

- Comma-separated, with header; missing values written as `NA`.
- One row per CpG. CpGs not measured by a platform or sample are `NA`. Each script selects its CpG set itself (e.g. the common five-platform set in 12, or the Tier 1/Tier 2 consensus sets in 13), so the matrices may contain the union of all CpGs.
- **Coordinate columns** `chr` and `start`:
  - `chr` with prefix (`chr1` … `chr22`, `chrX`, `chrY`).
  - `start` = 1-based position of the CpG cytosine on the forward strand, GRCh38/hg38, with both strands of a CpG merged into one row (as in `helpers.R::mergeCpGStrands()`).
- **Methylation columns** `<Method>_<Sample>`: β-value as a fraction in **[0, 1]** (not percent).
- **Coverage columns** `<Method>_cov_<Sample>`: number of reads (integer) – sequencing platforms only.
- **EPIC columns** `EPIC_<Sample>`: normalized β-value in [0, 1], no coverage column. Required in `ALL.csv` for Blood1–5 and Fibro1–5; `EPIC_GIAB1`/`EPIC_GIAB2` are used by script 10.

**Naming rules:**

- `<Method>` ∈ `ONT`, `WGEC`, `TWIST`, `RRBS`, `PacBio`, `EPIC` (case-sensitive).
- `<Sample>` ∈ `Blood1`–`Blood5`, `Fibro1`–`Fibro5`, `GIAB1` (HG001), `GIAB2` (HG002).
- **Pairing:** the paired models (`~ subject + group`) take the subject from the trailing number, so `Blood<i>` and `Fibro<i>` **must** come from the same individual.

Example header of `ALL.csv` (excerpt):

```
chr,start,ONT_Blood1,ONT_cov_Blood1,WGEC_Blood1,WGEC_cov_Blood1,...,PacBio_GIAB2,PacBio_cov_GIAB2,EPIC_Blood1,...,EPIC_GIAB2
chr1,10469,0.83,14,0.91,22,...,NA,NA,NA,...,NA
```

### Statistics tables (`$STATS`)

Both tables are tab-separated.

#### `cpg_stats_methylbench.tab` – fully reproducible

Per-sample counts of overlapping CpGs at increasing coverage thresholds, computed across all methods. It is generated from the raw per-CpG methylation files:

```bash
Rscript scripts/R/07_generate_cpg_stats.R \
  --samplesheet [samplesheet.tsv] \
  --datadir     data/ \
  --outdir      "$STATS/"
```

#### `qc_stats_methylbench.tab` – manually compiled

Per-sample × per-method QC metrics, assembled from the QC reports of each tool:

| Column | Source | Tool / File |
|---|---|---|
| `Mean_meth_general` | Global mean methylation | Bismark summary report / modkit stats |
| `Mean_meth_overlapped` | Mean methylation at overlapping CpGs | Computed from merged matrices |
| `Mean_meth_10x` | Mean methylation at ≥10× CpGs | Computed from merged matrices |
| `Passed_reads` | Fraction of reads passing QC | nf-core/methylseq MultiQC report (RRBS/WGEC/TWIST), ToulligQC report (ONT) |
| `Mean_readlength_passed` | Mean read length of passing reads | MultiQC (short-read), ToulligQC `.data` report (ONT) |
| `Mean_Cov` | Mean CpG-level coverage | Bismark coverage report / modkit / pb-CpG-tools |
| `Insert_size` | Mean insert size (TWIST only) | Picard InsertSizeMetrics via MultiQC |
| `Unique_alignments` | Number of uniquely aligned reads | Bismark alignment report / modkit stats |
| `Mean_Genome_Cov` | Mean genome-wide coverage (ONT only) | samtools coverage summary |

> **Note:** Because `qc_stats_methylbench.tab` aggregates heterogeneous per-tool reports without a common machine-readable format, it was compiled manually and is not generated by this pipeline.

### Inputs for script 17 (imprinted gDMRs)

Script 17 evaluates each platform on its own CpGs and therefore reads the per-platform methylation calls rather than the merged matrices. It is optional in `run_methylbench.sh` and runs only if `CALLS_DIR` and `SAMPLESHEET` are set; `GDMR_BED` defaults to `data/imprinted_gDMRs_hg38.bed`.

| Input | Variable / option | Format |
|---|---|---|
| Per-platform calls | `CALLS_DIR` / `--calls_dir` | `<CALLS_DIR>/<METHOD>/<Sample><suffix>`, one file per platform and sample (see below) |
| Sample sheet | `SAMPLESHEET` / `--samplesheet` | comma-separated, columns `Sample` (e.g. `Blood1`, `GIAB1`) and `Sampleset` (`Blood`, `Fibroblast`, `GIAB`) |
| Imprinted gDMRs | `GDMR_BED` / `--gdmr_bed` | BED, hg38, 0-based; columns `chr`, `start`, `end`, `name`; default in `run_methylbench.sh`: `data/imprinted_gDMRs_hg38.bed` (see below) |
| EPIC (optional) | `EPIC_TSV` / `--epic_path` | comma-separated list of tab-separated EPIC β-matrices with a `coord` column (`chr:pos`) and one column `EPIC_<Sample>` per sample |

Per-platform call files and default suffixes (changeable via `--ont_suffix`, `--bismark_suffix`, `--rrbs_suffix`, `--pacbio_suffix`):

| Platform | File | Format |
|---|---|---|
| ONT | `ONT/<Sample>_modkit_pileup.bed` | modkit pileup bedMethyl (`--cpg --combine-strands`), 0-based |
| WGEC, TWIST, RRBS (GIAB) | `<METHOD>/<Sample>.bismark.cov.gz` | Bismark coverage; per-strand or `coverage2cytosine --merge_CpG` output is detected automatically |
| RRBS (blood/fibroblast) | `RRBS/<Sample>.methylation.tab` | space-separated `chr pos coverage n_unmeth fraction`, 1-based, per strand; used if present, otherwise the Bismark file |
| PacBio | `PacBio/<Sample>.combined.bed` | pb-CpG-tools combined BED, 0-based |

All positions are converted to 1-based CpG positions (C of the + strand). With `--dss_dir`, script 17 also counts the DMCs of the primary analysis that fall into gDMRs (internal control, Suppl. Table S12E); it therefore runs after step `dss`.

#### gDMR list (`data/imprinted_gDMRs_hg38.bed`)

The 67 imprinted germline DMRs were derived from the compilation of Akbari et al. (2022, *eLife* 11:e77898, Supplementary file 1), which lists 383 intervals with parent-of-origin methylation from four genome-wide studies (Court et al. 2014; Joshi et al. 2016; Hernández Mora et al. 2018; Zink et al. 2018) in hg38 coordinates. Intervals reported by at least two of the four studies were retained (n = 68), and overlapping intervals were merged (n = 67; median width 1,278 bp). The fourth column gives the gDMR name used in Suppl. Table S12.

---

## Running the analysis

### Option A – complete run (recommended)

Run from the repository root:

```bash
conda activate methylbench

export DATA=/path/to/matrices        # the five .csv files described above
export STATS=/path/to/stats          # cpg_stats_methylbench.tab, qc_stats_methylbench.tab
export OUT=/path/to/results

bash scripts/bash/run_methylbench.sh
```

The script:
- checks that `$DATA`, `$STATS` and `$OUT` are set and that all seven input files exist;
- creates the output folders;
- runs every step in the correct order, writing one log per step to `$OUT/logs/<step>.log`;
- stops at the first failing step.

Optional variables:

| Variable | Default | Meaning |
|---|---|---|
| `N_REPS` | `5` | downsampling replicates (script 14), as in the paper |
| `CORES` | `4` | parallel workers (script 14) |
| `SKIP` | – | steps to skip, e.g. `SKIP="unpaired downsampling"` |
| `CALLS_DIR`, `SAMPLESHEET` | – | inputs for script 17; the step `imprinting` is skipped if one of them is missing |
| `GDMR_BED` | `data/imprinted_gDMRs_hg38.bed` | gDMR list for script 17 |
| `EPIC_TSV` | – | optional EPIC matrices for script 17 |

Step names: `qc`, `correlation`, `density`, `pca`, `common_set`, `dss`, `unpaired`, `downsampling`, `annotation_CpG`, `annotation_DMR`, `dmr_probe_density`, `coverage_context`, `imprinting`.

### Option B – individual steps

The same commands, step by step. Run them in this order, because later steps use earlier results:
- `dss` → `downsampling`, `annotation`, `dmr_probe_density`, `coverage_context` and `imprinting`

```bash
conda activate methylbench

export DATA=/path/to/matrices
export STATS=/path/to/stats
export OUT=/path/to/results
export N_REPS=5

mkdir -p "$OUT"/{figures,figures_unpaired,diff_meth,dmr_dss,dmr_dss_unpaired,downsampling,annotation_enrichment,dmr_density,coverage_context,imprinting}

# Figure 3, Suppl. Figure 1, Suppl. Table S6 – QC and per-CpG coverage uniformity
Rscript scripts/R/08_qc_visualization.R \
  --cpg_stats "$STATS/cpg_stats_methylbench.tab" \
  --qc_stats  "$STATS/qc_stats_methylbench.tab" \
  --all_path  "$DATA/ALL_without_EPIC.csv" \
  --outdir    "$OUT/figures/"

# Figure 4, Suppl. Figures 2 and 3, Suppl. Table S7 – cross-platform correlation
Rscript scripts/R/09_correlation_analysis.R \
  --datadir "$DATA/" \
  --outdir  "$OUT/figures/"

# Figure 5, Suppl. Figure 5 – methylation density
Rscript scripts/R/10_density_plots.R \
  --datadir   "$DATA/" \
  --outdir    "$OUT/figures/" \
  --epic_path "$DATA/ALL.csv"

# Figure 6, Suppl. Figure 6 – PCA of sample–platform profiles
Rscript scripts/R/11_pca.R \
  --all_path   "$DATA/ALL.csv" \
  --blood_path "$DATA/Blood_without_EPIC.csv" \
  --fibro_path "$DATA/Fibro_without_EPIC.csv" \
  --giab_path  "$DATA/GIAB_without_EPIC.csv" \
  --outdir     "$OUT/figures/"

# Figure 7, Suppl. Figure 7 – measurement properties on the common CpG set (no testing)
Rscript scripts/R/12_differential_methylation.R \
  --all_path   "$DATA/ALL.csv" \
  --blood_path "$DATA/Blood_without_EPIC.csv" \
  --fibro_path "$DATA/Fibro_without_EPIC.csv" \
  --outdir     "$OUT/figures/" \
  --datadir    "$OUT/diff_meth/"

# Figures 8, 9A/B/D, Suppl. Figures 8 and 10 – paired DSS + DMRcate (Tier 1 / Tier 2)
Rscript scripts/R/13_DMR_DSS_analysis.R \
  --seq_path "$DATA/ALL_without_EPIC.csv" \
  --all_path "$DATA/ALL.csv" \
  --outdir   "$OUT/figures/" \
  --datadir  "$OUT/dmr_dss/"

# Suppl. Table S2 – unpaired sensitivity run (design ~ group; tables only)
Rscript scripts/R/13_DMR_DSS_analysis.R \
  --seq_path "$DATA/ALL_without_EPIC.csv" \
  --all_path "$DATA/ALL.csv" \
  --outdir   "$OUT/figures_unpaired/" \
  --datadir  "$OUT/dmr_dss_unpaired/" \
  --unpaired

# Suppl. Figure 9, Suppl. Table S5 – depth-matched sensitivity analysis
# (reads BSseq_Tier1.rds, Tier1_DML.rds, Tier1_DMR.rds from the paired run)
Rscript scripts/R/14_downsampling_sensitivity.R \
  --dss_dir "$OUT/dmr_dss/" \
  --outdir  "$OUT/figures/" \
  --datadir "$OUT/downsampling/" \
  --n_reps  "$N_REPS" \
  --cores   4

# Figure 9C, Suppl. Figure 11, Suppl. Tables S3/S4 – annotation enrichment
Rscript scripts/R/15_annotation_enrichment_background.R \
  --dss_dir "$OUT/dmr_dss/" \
  --outdir  "$OUT/figures/" \
  --datadir "$OUT/annotation_enrichment/" \
  --level   CpG

Rscript scripts/R/15_annotation_enrichment_background.R \
  --dss_dir "$OUT/dmr_dss/" \
  --outdir  "$OUT/figures/" \
  --datadir "$OUT/annotation_enrichment/" \
  --level   DMR

# Suppl. Figure 12, Suppl. Table S11 – probe density and DMR calling
# (reads the consensus sets and DML/DMR results of the paired run)
Rscript scripts/R/16_dmr_probe_density.R \
  --dss_dir  "$OUT/dmr_dss/" \
  --all_path "$DATA/ALL.csv" \
  --outdir   "$OUT/figures/" \
  --datadir  "$OUT/dmr_density/"

# Suppl. Figures 1E,F and 12E, Suppl. Tables S6B and S11F – coverage by genomic context
# and DMC/DMR sensitivity (reads the Tier 1 results of the paired run)
Rscript scripts/R/18_coverage_by_context.R \
  --all_path "$DATA/ALL_without_EPIC.csv" \
  --dss_dir  "$OUT/dmr_dss/" \
  --outdir   "$OUT/figures/" \
  --datadir  "$OUT/coverage_context/"

# Suppl. Figure 4, Suppl. Table S12 – imprinted gDMRs (inputs: see "Inputs for script 17";
# omit --epic_path if no EPIC matrices are available)
Rscript scripts/R/17_imprinted_regions.R \
  --samplesheet "$SAMPLESHEET" \
  --calls_dir   "$CALLS_DIR" \
  --gdmr_bed    data/imprinted_gDMRs_hg38.bed \
  --epic_path   "$EPIC_TSV" \
  --dss_dir     "$OUT/dmr_dss/" \
  --outdir      "$OUT/figures/" \
  --datadir     "$OUT/imprinting/"
```

> **Note:** In multi-line commands, the backslash must be the **last** character of the line. A trailing space after `\` breaks the command.

`Rscript <script> --help` lists all options, including the statistical parameters. Defaults correspond to the manuscript: FDR < 0.05, |Δβ| ≥ 0.1, ≥ 10× in ≥ 4/5 samples per tissue, DMRcate λ = 1000 and C = 2 with ≥ 3 CpGs.

### Output overview

| Folder | Content |
|---|---|
| `$OUT/figures/` | all figure panels (PNG; Suppl. Figure 9 also as single PDF panels) |
| `$OUT/diff_meth/` | `DeltaBeta_common_set.csv`, `DeltaBeta_pairwise_correlation.csv` (Figure 7A), `Variances.csv` |
| `$OUT/dmr_dss/` | Tier 1/2 consensus sets, BSseq objects, DML/DMC/DMR tables, Jaccard and AUC tables |
| `$OUT/dmr_dss_unpaired/` | the same tables without the subject term (Suppl. Table S2) |
| `$OUT/downsampling/` | per-replicate results, `Table_S5_downsampling_*.tsv` |
| `$OUT/annotation_enrichment/` | `annotation_enrichment_{genecentric,cpgstructural}_{CpG,DMR}_{Tier1,Tier2}.tsv` (Suppl. Tables S3/S4, Additional file 1) |
| `$OUT/dmr_density/` | `TableS11A`–`TableS11E_*.tsv` (Suppl. Table S11) |
| `$OUT/coverage_context/` | `TableS6B_coverage_by_context.tsv`, `TableS11F_sensitivity_by_context.tsv`, per-sample and per-CpG tables |
| `$OUT/imprinting/` | `TableS12A`–`TableS12E_*.tsv` (Suppl. Table S12), per-gDMR summaries |
| `$OUT/figures/*.tsv/.tab` | `Coverage_Uniformity_CV_summary.tab` (Suppl. Table S6), `Table_S7_intermediate_correlation.tsv`, `Correlation_by_methylation_stratum.tsv` |

---

## Data Availability

Processed methylation matrices are available from the corresponding author upon request.

The PacBio HiFi data of the GIAB reference samples are publicly available: HG001/NA12878 from the [Platinum Pedigree dataset](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets) and HG002/NA24385 from the [PacBio website](https://www.pacb.com/connect/datasets/) (SPRQ-Nx dataset).

---

## License

This project is licensed under the MIT License – see [LICENSE](LICENSE) for details.

---

## Contact

For questions regarding the analysis code, please open a GitHub issue or contact the corresponding authors:
**Julia Schulze-Hentrich** – Department of Genetics, Saarland University
**Lukas Laufer** – Department of Genetics, Saarland University
