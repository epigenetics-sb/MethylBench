#!/usr/bin/env bash

# MethylBench - PacBio Methylation Calling
# Tool:    pb-cpg-tools (aligned_bam_to_cpg_scores)
# Input:   Aligned HiFi BAM file (PacBio, pbmm2-aligned to GRCh38)
# Output:  Per-CpG methylation scores (bed + bigwig)
# Usage:   bash 03_pbcpgtools.sh <sample_id> <bam_file> <output_dir>
# Note:    Linux x86_64 only
#          Run in methylbench-pacbio conda environment

set -euo pipefail

SAMPLE_ID="${1:?ERROR: sample_id required as \$1}"
BAM_FILE="${2:?ERROR: bam_file required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

CPGTOOLS_BIN="/path/to/cpgtools/bin/aligned_bam_to_cpg_scores"
MODEL="/path/to/cpgtools/models/pileup_calling_model.v1.tflite"
OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE_ID}.GRCh38.pbmm2"
THREADS=8

mkdir -p "${OUTPUT_DIR}"

"${CPGTOOLS_BIN}" \
    --bam "${BAM_FILE}" \
    --output-prefix "${OUTPUT_PREFIX}" \
    --model "${MODEL}" \
    --threads "${THREADS}"