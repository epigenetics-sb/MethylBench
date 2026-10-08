#!/usr/bin/env bash
# =============================================================================
# MethylBench - short-read preprocessing (nf-core/methylseq 1.0.0)
# =============================================================================
# WGEC and TWIST (standard mode) and RRBS (--rrbs); runs with the Singularity
# profile, requires Nextflow and Singularity/Apptainer on the system.
#
# Usage:
#   bash scripts/bash/04_methylseq.sh <wgec|twist|rrbs> <samplesheet.csv> <output_dir>
# =============================================================================

set -euo pipefail

MODE="${1:?ERROR: mode required as \$1 [wgec|twist|rrbs]}"
SAMPLESHEET="${2:?ERROR: samplesheet (CSV) required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"

GENOME="GRCh38"
PIPELINE_VERSION="1.0.0"
WORK_DIR="${OUTPUT_DIR}/work"
SINGULARITY_CACHE="${HOME}/.singularity/cache"

EXTRA_FLAGS=""
if [ "${MODE}" = "rrbs" ]
    then
        EXTRA_FLAGS="--rrbs"
fi

mkdir -p "${OUTPUT_DIR}" "${SINGULARITY_CACHE}"

export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"

nextflow run nf-core/methylseq \
    -r "${PIPELINE_VERSION}" \
    -profile singularity \
    --input "${SAMPLESHEET}" \
    --genome "${GENOME}" \
    --outdir "${OUTPUT_DIR}" \
    -work-dir "${WORK_DIR}" \
    ${EXTRA_FLAGS}
