# MethylBench

> Reproducible analysis code for a systematic benchmark of six DNA methylation profiling technologies across diverse sequencing platforms.

**Associated publication:**
Laufer L, Gasparoni G, Hentrich T, Sofan L, Admard J, Buena-Atienza E, Pogoda M, Ossowski S, Casadei N, Rieß O, Haack TB, Buchert R, Schulze-Hentrich J.
*MethylBench: A comprehensive benchmark of DNA methylation profiling methods across diverse sequencing platforms.*
[JOURNAL].

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
│   ├── python/
│   │   ├── 05_parse_toulligqc.py	       	# Summarize over ToulligQC .data files into one QC table        
│   └── R/
│       ├── 06_visualize_toulligqc_summary.R 	# Visualization for ONT QC reports. 
│       ├── 07_generate_cpg_stats.R		# Summarize CpG information for further analysis.
│       ├── 08_qc_visualization.R 		# QC Visualization, Figure 3. 
│       ├── 09_correlation_analysis.R 		# Correlation Analysis, Figure 4. 
│       ├── 10_density_plots.R 			# Methylation Density Analysis, Figure 5. 
│       ├── 11_pca.R 				# Principal Component Analysis, Figure 6. 
│       ├── 12_differential_methylation.R 	# Differential Methylation Analysis, Figure 7. 
│       ├── 13_annotation.R 			# Visualization for ONT QC reports.   		
│       └── helpers.R				# Helper functionality.
│
├── envs/
│   ├── environment.yml				# Basic environment, Tools and Python utility
│   ├── ont.yml					# ONT related tools, modkit, toulligQC
│   ├── pacbio.yml				# PacBio specific tool, pb-cpg-tools
│   └── r_analysis.yml				# R-related packages for R-analysis
│
└── docs/
    └── reproduction_guide.md
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

# Create and activate base environment
cd envs/
conda env create -f environment.yml
conda activate methylbench
```

---

## Data Availability

Raw sequencing data generated in this study are deposited at the **European Genome-phenome Archive (EGA)** under accession number `EGASXXXXXXX` (available upon publication).

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

Each R script in `scripts/R/` corresponds directly to a figure in the manuscript:

```bash
# Figure 3 – QC metrics and coverage
Rscript scripts/R/08_qc_visualization.R

# Figure 4 – Cross-platform correlation
Rscript scripts/R/09_correlation_analysis.R

# Figure 5 – Methylation density distributions
Rscript scripts/R/10_density_plots.R

# Figure 6 – Principal component analysis
Rscript scripts/R/11_pca.R

# Figure 7 – Differential methylation analysis
Rscript scripts/R/12_differential_methylation.R

# Figure S14 – DMC annotation
Rscript scripts/R/13_annotation.R
```

All scripts expect preprocessed (methylation) matrices as input. See `docs/reproduction_guide.md` for detailed instructions.

---

## Citation

If you use this code, please cite:

```
Laufer et al. (2025). MethylBench: A comprehensive benchmark of DNA methylation
profiling methods across diverse sequencing platforms.
[JOURNAL]. DOI: [to be added upon publication]
```

---

## License

This project is licensed under the MIT License – see [LICENSE](LICENSE) for details.

---

## Contact

For questions regarding the analysis code, please open a GitHub issue or contact the corresponding author:
**Julia Schulze-Hentrich** – Department of Genetics, Saarland University
**Lukas Laufer** - Department of Genetics, Saarland University
