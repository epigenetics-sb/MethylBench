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
│   │   ├── 01_modkit_pileup.sh			# modkit, extract methylation information from aligned .bam files.
│   │   ├── 02_toulligqc.sh			# ToulligQC, perform QC analysis on alignend .bam files and create intermediary files.         
│   │   ├── 03_pbcpgtools.sh			# pb-cpg-tools, extract methylation information from PacBio alignment files.         		
│   │   └── 04_methylseq.sh			# nf-core/methylseq, run the nextflow methylseq pipeline for standard short-read data processing.         		
 ^t^b    ^t^b    ^t^t ^t^` ^t^` run_methylbench.sh              # Executes all R scripts (whole MethylBench Analysis) in one run.
│   ├── python/
│   │   ├── 05_parse_toulligqc.py	       	# Summarize over ToulligQC .data files into one QC table        
│   └── R/
│       ├── 06_visualize_toulligqc_summary.R 	# Visualization for ONT QC reports. 
│       ├── 07_generate_cpg_stats.R		# Summarize CpG information for further analysis.
│       ├── 08_qc_visualization.R 		# QC Visualization, Figure 3, Suppl. Figure 16 & Suppl. Table S6 (--all_path).
│       ├── 09_correlation_analysis.R 		# Correlation Analysis, Figure 4, Suppl. Figures 2 & 19, Suppl. Table S7.
│       ├── 10_density_plots.R 			# Methylation Density Analysis, Figure 5. 
│       ├── 11_pca.R 				# Principal Component Analysis (sample scores), Figure 6, Suppl. Figure 1.
│       ├── limma_diff_meth.R 			# Exploratory paired limma per platform (input for 12, Suppl. Figure 7).
│       ├── 12_differential_methylation.R 	# Exploratory limma & Wilcoxon analysis, Figure 7, Suppl. Figures 6-8.
│       ├── 13_DMR_DSS_analysis.R 		# Tier1/Tier2 consensus sets, paired DSS / EPIC limma (M-values), DMRcate;
│       │					# Figure 8, Figure 9A/B/D, Suppl. Figures 12-15, Suppl. Table S2 (--unpaired).
│       ├── 14_downsampling_sensitivity.R 	# Depth-matching sensitivity analysis on Tier1, Suppl. Figure 17, Suppl. Table S5.
│       ├── 15_annotation_enrichment_background.R 	# Background-corrected annotation enrichment (two hierarchies);
│       │					# Figure 9C, Suppl. Figures 9 and 18, Suppl. Tables S3/S4.
│       └── utils/
│           └── helpers.R			# Helper functionality.
│
└── envs/
    ├── ont.yml					# ONT related tools, modkit, toulligQC
    ├── pacbio.yml				# PacBio specific tool, pb-cpg-tools
    └── methylbench.yml				# Actual MethylBench Analysis
```

---

## Requirements

- [Conda](https://docs.conda.io/en/latest/) >= 23.x
- [Nextflow](https://www.nextflow.io/) >= 20.x
- [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) >= 3.x
- R >= 4.3
- Python >= 3.10

All tool-specific dependencies are managed via Conda environments defined in `envs/`.

---

## Quick Start

```bash
# Clone repository
git clone https://github.com/epigenetics-sb/MethylBench
cd MethylBench

# Create environments
cd envs/
conda env create -f ont.yml
conda env create -f pacbio.yml
conda env create -f methylbench.yml

cd ..
```
```bash
################################################################################################
# NOTE: Please check the scripts for proper argument structure and folder setup for all scripts!
################################################################################################

# Run preprocessing steps
# All bash scripts can be properly run, as the help descriptions tell you.
conda activate ont

bash scripts/bash/01_modkit_pileup.sh --help
bash scripts/bash/02_toulligqc.sh --help

conda deactivate 
conda activate pacbio

bash scripts/bash/03_pbcpgtools.sh --help

conda deactivate
conda activate methylbench

bash scripts/bash/04_methylseq.sh --help

# Run tool QC visualization
python3 scripts/python/05_parse_toulligqc.py

# Run actual analysis in R
Rscript scripts/R/06_visualize_toulligqc_summary.R
Rscript scripts/R/07_generate_cpg_stats.R
Rscript scripts/R/08_qc_visualization.R
Rscript scripts/R/09_correlation_analysis.R
Rscript scripts/R/10_density_plots.R
Rscript scripts/R/11_pca.R
Rscript scripts/R/limma_diff_meth.R          # must run before 12
Rscript scripts/R/12_differential_methylation.R
Rscript scripts/R/13_DMR_DSS_analysis.R      # paired run (default)
Rscript scripts/R/13_DMR_DSS_analysis.R --unpaired   # sensitivity run, separate --datadir
Rscript scripts/R/14_downsampling_sensitivity.R
Rscript scripts/R/15_annotation_enrichment_background.R --level CpG
Rscript scripts/R/15_annotation_enrichment_background.R --level DMR
```

---

## Data Availability

Processed methylation matrices are available from the corresponding author upon request.

GIAB reference samples (HG001/NA12878, HG002/NA24385) including PacBio methylation data are publicly available via the [PacBio website](https://www.pacb.com/connect/datasets/).

---

## Input Data: QC Statistics Tables

Two pre-computed summary tables are required as input for the downstream R analysis scripts. Both are provided in `data/stats/`:

### `cpg_stats_methylbench.tab` – fully reproducible
Per-sample counts of overlapping CpGs at increasing coverage thresholds, computed across all methods simultaneously. Generated programmatically from the raw per-CpG methylation files:

```bash
Rscript scripts/R/07_generate_cpg_stats.R \
  --samplesheet [samplesheet.tsv] \
  --datadir     data/ \
  --outdir      data/stats/
```

### `qc_stats_methylbench.tab` – manually compiled
Per-sample × per-method QC metrics. This table was assembled manually by extracting summary statistics from the QC reports of each tool:

| Column | Source | Tool / File |
|---|---|---|
| `Mean_meth_general` | Global mean methylation | Bismark summary report / modkit stats |
| `Mean_meth_overlapped` | Mean methylation at overlapping CpGs | Computed from merged matrices |
| `Mean_meth_10x` | Mean methylation at ≥10× CpGs | Computed from merged matrices |
| `Passed_reads` | Fraction of reads passing QC | nf-core/methylseq MultiQC report (RRBS/WGEC/TWIST), ToulligQC report (ONT) |
| `Mean_readlength_passed` | Mean read length of passing reads | MultiQC (short-read), ToulligQC `.data` report (ONT) |
| `Mean_Cov` | Mean CpG-level coverage | Bismark coverage report / modkit / pb-cpg-tools |
| `Insert_size` | Mean insert size (TWIST only) | Picard InsertSizeMetrics via MultiQC |
| `Unique_alignments` | Number of uniquely aligned reads | Bismark alignment report / modkit stats |
| `Mean_Genome_Cov` | Mean genome-wide coverage (ONT only) | samtools coverage summary |

> **Note:** Because `qc_stats_methylbench.tab` aggregates heterogeneous per-tool reports that do not share a common machine-readable format, it was compiled manually and is not auto-generated by this pipeline.

---

## Reproduce Individual Figures

Each R script in `scripts/R/` corresponds directly to figures and tables in the manuscript. Run `Rscript <script> --help` for all options.

```bash
# Figure 3 – QC metrics and coverage
Rscript scripts/R/08_qc_visualization.R

# Figure 4 – Cross-platform correlation
Rscript scripts/R/09_correlation_analysis.R

# Figure 5 – Methylation density distributions
Rscript scripts/R/10_density_plots.R

# Figure 6 and Suppl. Figure 1 – PCA of sample–platform profiles (scores, not loadings)
Rscript scripts/R/11_pca.R

# Figure 7, Suppl. Figures 6–8 – exploratory analysis on the common five-platform CpG set
# (paired limma on beta-values, design ~ subject + group; unpaired Wilcoxon rank-sum test).
# limma_diff_meth.R must be run first; its output directory is passed via --limma_dir.
Rscript scripts/R/limma_diff_meth.R --all_path data/matrices/ALL.csv --outdir results/limma/
Rscript scripts/R/12_differential_methylation.R --all_path data/matrices/ALL.csv \
  --limma_dir results/limma/ --outdir results/figures/ --datadir results/diff_meth/

# Figure 8, Figure 9A/B/D, Suppl. Figures 12–15 – primary analysis on the Tier1/Tier2
# consensus sets: paired Beta-Binomial model (DSS::DMLfit.multiFactor, ~ subject + group)
# for the sequencing platforms, paired limma on M-values for EPIC (Tier2), DMRcate on the
# paired per-CpG statistics. Tier 1 is built from the sequencing-only matrix,
# Tier 2 from the matrix with EPIC.
Rscript scripts/R/13_DMR_DSS_analysis.R \
  --seq_path data/matrices/ALL_without_EPIC.csv \
  --all_path data/matrices/ALL.csv \
  --outdir   results/figures/ \
  --datadir  results/dmr_dss/

# Suppl. Table S2 – sensitivity run without the subject term (design ~ group);
# writes tables only, figures are skipped. Use a separate --datadir.
Rscript scripts/R/13_DMR_DSS_analysis.R \
  --seq_path data/matrices/ALL_without_EPIC.csv \
  --all_path data/matrices/ALL.csv \
  --unpaired \
  --outdir   results/figures_unpaired/ \
  --datadir  results/dmr_dss_unpaired/

# Suppl. Figure 16, Suppl. Table S5 – depth-matching sensitivity analysis
# (needs BSseq_Tier1.rds from the paired run above)
Rscript scripts/R/14_downsampling_sensitivity.R \
  --bsseq_tier1 results/dmr_dss/BSseq_Tier1.rds

# Figure 9C, Suppl. Figures 9 and 17, Suppl. Tables S3/S4 – annotation enrichment relative to
# the tested background, with separate gene-centric and CpG-structural hierarchies.
# --level CpG: significant DMCs; --level DMR: tested CpGs located within DMRs.
Rscript scripts/R/15_annotation_enrichment_background.R --dss_dir results/dmr_dss/ \
  --outdir results/figures/ --datadir results/annotation/ --level CpG
Rscript scripts/R/15_annotation_enrichment_background.R --dss_dir results/dmr_dss/ \
  --outdir results/figures/ --datadir results/annotation/ --level DMR
```

All scripts expect preprocessed (methylation) matrices as input.

---

## Changes in the revised version

- **PCA (`11_pca.R`):** PCA is computed on sample–platform profiles (observations) × CpGs (variables) and sample scores are plotted; Method/Sample labels are parsed from the column names of the analysed matrix.
- **Paired design:** all differential methylation models include the subject term (`~ subject + group`): limma (`limma_diff_meth.R`), DSS via `DMLfit.multiFactor`/`DMLtest.multiFactor` and EPIC limma on M-values (`13_DMR_DSS_analysis.R`). DMRs are called with DMRcate from these paired statistics. `--unpaired` reproduces the analysis without the subject term.
- **Tier 1** is built from the sequencing-only matrix (`--seq_path`) and is no longer restricted to EPIC positions.
- **Annotation (`15_annotation_enrichment_background.R`):** significant sites are compared with the tested consensus background, using two independent hierarchies (gene-centric, CpG-structural) and Fisher's exact test (significant vs. non-significant tested CpGs). Replaces the former `13_annotation.R`.
- **Rank recovery (Suppl. Figure 15):** direction-aware score, `-log10(FDR) × directional agreement` with the reference platform.
- Scripts were renumbered; the former `14_DMR_DSS_analysis.R`, `16_downsampling_sensitivity.R` and `17_annotation_enrichment_background.R` are now 13, 14 and 15.

---

## License

This project is licensed under the MIT License – see [LICENSE](LICENSE) for details.

---

## Contact

For questions regarding the analysis code, please open a GitHub issue or contact the corresponding author:
**Julia Schulze-Hentrich** – Department of Genetics, Saarland University
**Lukas Laufer** - Department of Genetics, Saarland University
