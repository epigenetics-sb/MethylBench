#!/usr/bin/env bash
# =============================================================================
# MethylBench - run the complete R analysis (scripts 08-18) in one go
# =============================================================================
# Usage (from the repository root):
#
#   conda activate methylbench
#   export DATA=/path/to/matrices      # ALL.csv, ALL_without_EPIC.csv, ...
#   export STATS=/path/to/stats        # cpg_stats_methylbench.tab, qc_stats_methylbench.tab
#   export OUT=/path/to/results
#   bash scripts/bash/run_methylbench.sh
#
# Optional environment variables:
#   N_REPS   downsampling replicates for script 14   [default: 5, as in the paper]
#   CORES    parallel workers for script 14           [default: 4]
#   SKIP     space-separated step names to skip, e.g. SKIP="downsampling unpaired"
#
# Script 17 (imprinted gDMRs) needs per-platform methylation calls instead of
# the merged matrices and runs only if CALLS_DIR and SAMPLESHEET are set:
#   CALLS_DIR    directory with per-platform calls, <CALLS_DIR>/<METHOD>/<Sample><suffix>
#   SAMPLESHEET  sample sheet (comma-separated; columns Sample, Sampleset)
#   GDMR_BED     BED file of imprinted gDMRs (hg38, 0-based; chr, start, end, name)
#                [default: data/imprinted_gDMRs_hg38.bed, shipped with the repository]
#   EPIC_TSV     optional: comma-separated EPIC beta matrices (column 'coord' = chr:pos)
#
# Every step writes a log to $OUT/logs/<step>.log. The script stops at the
# first failing step (set -e) and reports which one failed.
# =============================================================================

set -euo pipefail

N_REPS="${N_REPS:-5}"
CORES="${CORES:-4}"
SKIP="${SKIP:-}"

# ---- 0. Checks --------------------------------------------------------------
for v in DATA STATS OUT; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: environment variable \$$v is not set (see header of this script)." >&2
    exit 1
  fi
done

if [[ ! -f scripts/R/utils/helpers.R ]]; then
  echo "ERROR: run this script from the repository root (scripts/R/utils/helpers.R not found)." >&2
  exit 1
fi

required_inputs=(
  "$DATA/ALL.csv"
  "$DATA/ALL_without_EPIC.csv"
  "$DATA/Blood_without_EPIC.csv"
  "$DATA/Fibro_without_EPIC.csv"
  "$DATA/GIAB_without_EPIC.csv"
  "$STATS/cpg_stats_methylbench.tab"
  "$STATS/qc_stats_methylbench.tab"
)
missing=0
for f in "${required_inputs[@]}"; do
  if [[ ! -s "$f" ]]; then echo "ERROR: missing input file: $f" >&2; missing=1; fi
done
[[ $missing -eq 0 ]] || { echo "See README, section 'Input data', for the expected files." >&2; exit 1; }

mkdir -p "$OUT"/{figures,figures_unpaired,diff_meth,dmr_dss,dmr_dss_unpaired,downsampling,annotation_enrichment,dmr_density,imprinting,logs}

# ---- helper -----------------------------------------------------------------
run_step() {
  local name="$1"; shift
  if [[ " $SKIP " == *" $name "* ]]; then
    echo "[skip] $name"
    return 0
  fi
  echo "[$(date '+%H:%M:%S')] >>> $name"
  if ! "$@" > "$OUT/logs/$name.log" 2>&1; then
    echo "ERROR: step '$name' failed - see $OUT/logs/$name.log" >&2
    tail -n 20 "$OUT/logs/$name.log" >&2
    exit 1
  fi
}

# ---- Figure 3, Suppl. Figure 1, Suppl. Table S6 - QC and coverage uniformity -
run_step qc \
  Rscript scripts/R/08_qc_visualization.R \
    --cpg_stats "$STATS/cpg_stats_methylbench.tab" \
    --qc_stats  "$STATS/qc_stats_methylbench.tab" \
    --all_path  "$DATA/ALL_without_EPIC.csv" \
    --outdir    "$OUT/figures/"

# ---- Figure 4, Suppl. Figures 2 and 3, Suppl. Table S7 - correlation --------
run_step correlation \
  Rscript scripts/R/09_correlation_analysis.R \
    --datadir "$DATA/" \
    --outdir  "$OUT/figures/"

# ---- Figure 5, Suppl. Figure 5 - methylation density ------------------------
run_step density \
  Rscript scripts/R/10_density_plots.R \
    --datadir   "$DATA/" \
    --outdir    "$OUT/figures/" \
    --epic_path "$DATA/ALL.csv"

# ---- Figure 6, Suppl. Figure 6 - PCA ----------------------------------------
run_step pca \
  Rscript scripts/R/11_pca.R \
    --all_path   "$DATA/ALL.csv" \
    --blood_path "$DATA/Blood_without_EPIC.csv" \
    --fibro_path "$DATA/Fibro_without_EPIC.csv" \
    --giab_path  "$DATA/GIAB_without_EPIC.csv" \
    --outdir     "$OUT/figures/"

# ---- Figure 7, Suppl. Figure 7 - measurement properties (no testing) -------
run_step common_set \
  Rscript scripts/R/12_differential_methylation.R \
    --all_path   "$DATA/ALL.csv" \
    --blood_path "$DATA/Blood_without_EPIC.csv" \
    --fibro_path "$DATA/Fibro_without_EPIC.csv" \
    --outdir     "$OUT/figures/" \
    --datadir    "$OUT/diff_meth/"

# ---- Figures 8, 9A/B/D, Suppl. Figures 8 and 10 - paired DSS / DMRcate ------
run_step dss \
  Rscript scripts/R/13_DMR_DSS_analysis.R \
    --seq_path "$DATA/ALL_without_EPIC.csv" \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/figures/" \
    --datadir  "$OUT/dmr_dss/"

# ---- Suppl. Table S2 - unpaired sensitivity run (tables only) ---------------
run_step unpaired \
  Rscript scripts/R/13_DMR_DSS_analysis.R \
    --seq_path "$DATA/ALL_without_EPIC.csv" \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/figures_unpaired/" \
    --datadir  "$OUT/dmr_dss_unpaired/" \
    --unpaired

# ---- Suppl. Figure 9, Suppl. Table S5 - depth-matched sensitivity -----------
# needs BSseq_Tier1.rds, Tier1_DML.rds and Tier1_DMR.rds from step 'dss'
run_step downsampling \
  Rscript scripts/R/14_downsampling_sensitivity.R \
    --dss_dir "$OUT/dmr_dss/" \
    --outdir  "$OUT/figures/" \
    --datadir "$OUT/downsampling/" \
    --n_reps  "$N_REPS" \
    --cores   "$CORES"

# ---- Figure 9C, Suppl. Figure 11, Suppl. Tables S3/S4 - annotation ---------
# needs Tier{1,2}_consensus_CpGs.tsv and the DML/DMR tables from step 'dss'
run_step annotation_CpG \
  Rscript scripts/R/15_annotation_enrichment_background.R \
    --dss_dir "$OUT/dmr_dss/" \
    --outdir  "$OUT/figures/" \
    --datadir "$OUT/annotation_enrichment/" \
    --level   CpG

run_step annotation_DMR \
  Rscript scripts/R/15_annotation_enrichment_background.R \
    --dss_dir "$OUT/dmr_dss/" \
    --outdir  "$OUT/figures/" \
    --datadir "$OUT/annotation_enrichment/" \
    --level   DMR

# ---- Suppl. Figure 12, Suppl. Table S11 - probe density and DMR calling ----
# needs the Tier 1/Tier 2 consensus sets and DML/DMR results from step 'dss'
run_step dmr_probe_density \
  Rscript scripts/R/16_dmr_probe_density.R \
    --dss_dir  "$OUT/dmr_dss/" \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/figures/" \
    --datadir  "$OUT/dmr_density/"

# ---- Suppl. Figures 1E,F and 12E, Suppl. Tables S6B and S11F ----------------
# coverage by genomic context and DMC/DMR sensitivity (needs step 'dss')
run_step coverage_context \
  Rscript scripts/R/18_coverage_by_context.R \
    --all_path "$DATA/ALL_without_EPIC.csv" \
    --dss_dir  "$OUT/dmr_dss/" \
    --outdir   "$OUT/figures/" \
    --datadir  "$OUT/coverage_context/"

# ---- Suppl. Figure 4, Suppl. Table S12 - imprinted gDMRs (optional) ---------
# needs per-platform calls (CALLS_DIR) and SAMPLESHEET; GDMR_BED defaults to the
# gDMR list shipped with the repository. Uses the DML results from step 'dss'
# for the DMCs located in gDMRs (Suppl. Table S12E).
GDMR_BED="${GDMR_BED:-data/imprinted_gDMRs_hg38.bed}"
if [[ -n "${CALLS_DIR:-}" && -n "${SAMPLESHEET:-}" && -f "$GDMR_BED" ]]; then
  imprint_args=(
    --samplesheet "$SAMPLESHEET"
    --calls_dir   "$CALLS_DIR"
    --gdmr_bed    "$GDMR_BED"
    --dss_dir     "$OUT/dmr_dss/"
    --outdir      "$OUT/figures/"
    --datadir     "$OUT/imprinting/"
  )
  if [[ -n "${EPIC_TSV:-}" ]]; then imprint_args+=(--epic_path "$EPIC_TSV"); fi
  run_step imprinting Rscript scripts/R/17_imprinted_regions.R "${imprint_args[@]}"
else
  echo "[skip] imprinting (set CALLS_DIR and SAMPLESHEET, and GDMR_BED if not using data/imprinted_gDMRs_hg38.bed, to run script 17)"
fi

echo "[$(date '+%H:%M:%S')] MethylBench analysis finished. Results: $OUT"
