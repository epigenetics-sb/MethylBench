#!/usr/bin/env bash

# MethylBench - ONT Quality Control
# Tool:    ToulligQC
# Input:   BAM files per ONT sample
# Output:  HTML report + data report per sample
# Usage:   bash 01_toulligqc.sh <sample_id> <bam_file> <output_dir>

set -euo pipefail

SAMPLE_ID="${1:?ERROR: sample_id required as \$1}"
BAM_FILE="${2:?ERROR: bam_file required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

THREADS=128
HTML_REPORT="${OUTPUT_DIR}/${SAMPLE_ID}_toulligqc_report.html"
DATA_REPORT="${OUTPUT_DIR}/${SAMPLE_ID}_toulligqc_report.data"

mkdir -p "${OUTPUT_DIR}"

toulligqc \
    --report-name "${SAMPLE_ID}" \
    --bam "${BAM_FILE}" \
    --html-report-path "${HTML_REPORT}" \
    --data-report-path "${DATA_REPORT}" \
    --thread "${THREADS}"