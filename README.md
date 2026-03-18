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
│
├── scripts/
│   ├── bash/
│   │   ├── 01_modkit_pileup.sh			# modkit, extract methylation information from aligned .bam files.
│   │   ├── 02_toulligqc.sh			# ToulligQC, perform QC analysis on alignend .bam files and create intermediary files.         
│   │   ├── 03_pbcpgtools.sh			# pb-cpg-tools, extract methylation information from PacBio alignment files.         		
│   │   └── 04_methylseq.sh			# nf-core/methylseq, run the nextflow methylseq pipeline for standard short-read data processing.         		
│   ├── python/
│   │   ├── parse_toulligqc.py	       		# Summarize over ToulligQC .data files into one QC table        
│   └── R/
│       ├── 01_epic_preprocessing.R    		
│       ├── 02_qc_visualization.R      		
│       ├── 03_correlation_analysis.R  		
│       ├── 04_density_plots.R         		
│       ├── 05_pca.R                   		
│       ├── 06_differential_methylation.R  	
│       ├── 07_annotation.R            		
│       └── utils/
│           ├── plot_theme.R
│           └── helpers.R
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

## Reproduce Individual Figures

Each R script in `scripts/R/` corresponds directly to a figure in the manuscript:

```bash
# Figure 3 – QC metrics and coverage
Rscript scripts/R/02_qc_visualization.R

# Figure 4 – Cross-platform correlation
Rscript scripts/R/03_correlation_analysis.R

# Figure 5 – Methylation density distributions
Rscript scripts/R/04_density_plots.R

# Figure 6 – Principal component analysis
Rscript scripts/R/05_pca.R

# Figure 7 – Differential methylation analysis
Rscript scripts/R/06_differential_methylation.R

# Figure S14 – DMC annotation
Rscript scripts/R/07_annotation.R
```

All scripts expect preprocessed methylation matrices as input. See `docs/reproduction_guide.md` for detailed instructions.

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
