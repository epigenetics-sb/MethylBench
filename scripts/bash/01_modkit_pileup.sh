#!/usr/bin/env bash
# =============================================================================
# MethylBench - ONT methylation calling (modkit pileup)
# =============================================================================
# Per-CpG 5mC calls from aligned ONT BAM files (GRCh38), strands combined.
# 5mC and 5hmC are both requested, so that, for data with separate 5hmC calls
# (blood and fibroblasts, R9.4.1), reads called as 5hmC count as valid but not
# 5mC-modified in the 5mC rows (code m) used downstream; for 5mC-only data
# (GIAB, R10.4.1) only 5mC rows are written.
# Environment: methylbench-ont (modkit 0.6.1, as used for the manuscript)
#
# Usage:
#   bash scripts/bash/01_modkit_pileup.sh <sample_id> <bam_file> <output_dir>
# =============================================================================

set -euo pipefail

SAMPLE_ID="${1:?ERROR: sample_id required as \$1}"
BAM_FILE="${2:?ERROR: bam_file required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

REFERENCE="/path/to/GRCh38.fa"   # adjust
THREADS=128
LOG_DIR="${OUTPUT_DIR}/logs"
OUTPUT_BED="${OUTPUT_DIR}/${SAMPLE_ID}_modkit_pileup.bed"
LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_modkit_pileup.log"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

modkit pileup \
    --reference "${REFERENCE}" \
    --modified-bases 5mC 5hmC \
    --cpg \
    --combine-strands \
    --log-filepath "${LOG_FILE}" \
    --threads "${THREADS}" \
    "${BAM_FILE}" \
    "${OUTPUT_BED}"
