#!/usr/bin/env bash

# MethylBench - ONT Methylation Calling
# Tool:    modkit pileup
# Input:   Aligned BAM file (ONT), reference genome (GRCh38)
# Output:  BED file with per-CpG 5mC and 5hmC methylation calls
# Usage:   bash 02_modkit_pileup.sh <sample_id> <bam_file> <output_dir>
# Note:    Requires modkit >= 0.6.1
#          Run in methylbench-ont conda environment

set -euo pipefail

SAMPLE_ID="${1:?ERROR: sample_id required as \$1}"
BAM_FILE="${2:?ERROR: bam_file required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

REFERENCE="/path/to/GRCh38.fa"          
THREADS=128
LOG_DIR="${OUTPUT_DIR}/logs"
OUTPUT_BED="${OUTPUT_DIR}/${SAMPLE_ID}_modkit_pileup.bed"
LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_modkit_pileup.log"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

modkit pileup \
    --reference "${REFERENCE}" \
    --modified-bases 5mC \
    --cpg \
    --combine-strands \
    --log-filepath "${LOG_FILE}" \
    --threads "${THREADS}" \
    "${BAM_FILE}" \
    "${OUTPUT_BED}"