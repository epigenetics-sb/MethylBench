#!/usr/bin/env bash
# =============================================================================
# MethylBench - ONT methylation calling (modkit pileup)
# =============================================================================
# Per-CpG 5mC calls from aligned ONT BAM files (GRCh38), strands combined.
# Environment: methylbench-ont (modkit >= 0.6.1)
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
    --modified-bases 5mC \
    --cpg \
    --combine-strands \
    --log-filepath "${LOG_FILE}" \
    --threads "${THREADS}" \
    "${BAM_FILE}" \
    "${OUTPUT_BED}"
