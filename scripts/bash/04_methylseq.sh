#!/usr/bin/env bash

# MethylBench - Short-Read Methylation Preprocessing
# Pipeline: nf-core/methylseq v1.0.0
# Covers:   WGEC, TWIST (standard mode), RRBS (--rrbs flag)
# Profile:  Singularity (no local tool installation required)
# Usage:    bash 04_methylseq.sh <mode> <samplesheet> <output_dir>
#           mode: "rrbs" for adding the extra rrbs flag to the nextflow pipeline

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