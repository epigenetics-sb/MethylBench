#!/usr/bin/env bash
# =============================================================================
# MethylBench - PacBio methylation calling (pb-CpG-tools)
# =============================================================================
# Per-CpG methylation scores from pbmm2-aligned HiFi BAM files (GRCh38).
# Environment: methylbench-pacbio (Linux x86_64 only)
#
# Usage:
#   bash scripts/bash/03_pbcpgtools.sh <sample_id> <bam_file> <output_dir>
# =============================================================================

set -euo pipefail

SAMPLE_ID="${1:?ERROR: sample_id required as \$1}"
BAM_FILE="${2:?ERROR: bam_file required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

CPGTOOLS_BIN="/path/to/cpgtools/bin/aligned_bam_to_cpg_scores"   # adjust
MODEL="/path/to/cpgtools/models/pileup_calling_model.v1.tflite"    # adjust
OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE_ID}.GRCh38.pbmm2"
THREADS=8

mkdir -p "${OUTPUT_DIR}"

"${CPGTOOLS_BIN}" \
    --bam "${BAM_FILE}" \
    --output-prefix "${OUTPUT_PREFIX}" \
    --model "${MODEL}" \
    --threads "${THREADS}"
