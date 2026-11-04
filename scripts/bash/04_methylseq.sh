#!/usr/bin/env bash
# =============================================================================
# MethylBench - short-read preprocessing (nf-core/methylseq 3.0.0)
# =============================================================================
# WGEC and TWIST: standard mode (trimming, Bismark alignment, deduplication,
#   methylation extraction).
# RRBS: the libraries were digested with HaeIII (not MspI), so the MspI-specific
#   option --rrbs is not used; deduplication is skipped (--skip_deduplication),
#   because RRBS reads start at the restriction sites and would otherwise be
#   removed as duplicates.
# Software versions of the runs used for the manuscript: nf-core/methylseq
#   3.0.0, Bismark 0.24.2, Trim Galore 0.6.10 (Cutadapt 4.9), Nextflow 24.10.3.
# Runs with the Singularity profile; requires Nextflow and Singularity/Apptainer
# on the system.
#
# Usage:
#   bash scripts/bash/04_methylseq.sh <wgec|twist|rrbs> <samplesheet.csv> <output_dir> [genome.fa]
#
#   genome.fa  optional GRCh38 FASTA (--fasta); without it, the iGenomes
#              reference GRCh38 is used (--genome GRCh38)
# =============================================================================

set -euo pipefail

MODE="${1:?ERROR: mode required as \$1 [wgec|twist|rrbs]}"
SAMPLESHEET="${2:?ERROR: samplesheet (CSV) required as \$2}"
OUTPUT_DIR="${3:?ERROR: output_dir required as \$3}"
FASTA="${4:-}"

PIPELINE_VERSION="3.0.0"
WORK_DIR="${OUTPUT_DIR}/work"
SINGULARITY_CACHE="${HOME}/.singularity/cache"

case "${MODE}" in
    wgec|twist) EXTRA_FLAGS=() ;;
    rrbs)       EXTRA_FLAGS=(--skip_deduplication) ;;
    *) echo "ERROR: mode must be wgec, twist or rrbs (got '${MODE}')" >&2; exit 1 ;;
esac

if [ -n "${FASTA}" ]; then
    REFERENCE=(--fasta "${FASTA}")
else
    REFERENCE=(--genome GRCh38)
fi

mkdir -p "${OUTPUT_DIR}" "${SINGULARITY_CACHE}"

export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"

nextflow run nf-core/methylseq \
    -r "${PIPELINE_VERSION}" \
    -profile singularity \
    --input "${SAMPLESHEET}" \
    "${REFERENCE[@]}" \
    --outdir "${OUTPUT_DIR}" \
    -work-dir "${WORK_DIR}" \
    ${EXTRA_FLAGS[@]+"${EXTRA_FLAGS[@]}"}
