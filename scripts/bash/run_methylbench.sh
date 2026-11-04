#!/usr/bin/env bash
# =============================================================================
# MethylBench - run the complete R analysis (scripts 07-18) in one go
# =============================================================================
# Usage (from the repository root):
#
#   conda activate methylbench
#   export DATA=/path/to/matrices      # ALL.csv, ALL_without_EPIC.csv, ... (written by 07b)
#   export STATS=/path/to/stats        # cpg_stats_methylbench.tab (written by 07),
#                                      # qc_stats_methylbench.tab (compiled manually)
#   export OUT=/path/to/results
#   export CALLS_SHEET=calls.csv       # per-platform calls: columns Sample, Method, File
#   export SAMPLESHEET=samplesheet.csv # columns Sample, Sampleset
#   export EPIC_TSV=EPIC_blood_fibro.csv,EPIC_GIAB.csv   # EPIC beta matrices (0-based)
#   bash scripts/bash/run_methylbench.sh
#
# With CALLS_SHEET, SAMPLESHEET and EPIC_TSV set, the run starts from the
# per-platform calls (WGEC, TWIST: Bismark coverage2cytosine --merge_CpG;
# RRBS: methylation table (.methylation.tab); ONT: modkit pileup; PacBio:
# pb-CpG-tools; formats: README, "Inputs for script 17") and first builds the matrices
# (step 'matrices', 07b) and cpg_stats_methylbench.tab (step 'cpg_stats', 07);
# script 17 (step 'imprinting') reads the same calls. Without them, the
# matrices and statistics already in $DATA and $STATS are used and script 17
# is skipped.
#
# Optional environment variables:
#   N_REPS   downsampling replicates for script 14   [default: 5, as in the paper]
#   CORES    parallel workers for script 14           [default: 4]
#   SKIP     space-separated step names to skip, e.g. SKIP="downsampling unpaired"
#   GDMR_BED BED file of imprinted gDMRs (hg38, 0-based; chr, start, end, name)
#            [default: data/imprinted_gDMRs_hg38.bed, shipped with the repository]
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

FROM_CALLS=0
if [[ -n "${CALLS_SHEET:-}" || -n "${SAMPLESHEET:-}" || -n "${EPIC_TSV:-}" ]]; then
  for v in CALLS_SHEET SAMPLESHEET EPIC_TSV; do
    if [[ -z "${!v:-}" ]]; then
      echo "ERROR: to start from the per-platform calls, set CALLS_SHEET, SAMPLESHEET and EPIC_TSV (\$$v is missing)." >&2
      exit 1
    fi
  done
  for f in "$CALLS_SHEET" "$SAMPLESHEET" ${EPIC_TSV//,/ }; do
    [[ -s "$f" ]] || { echo "ERROR: missing input file: $f" >&2; exit 1; }
  done
  FROM_CALLS=1
fi

mkdir -p "$DATA" "$STATS" "$OUT"/{figures,figures_unpaired,figures_unsmoothed,diff_meth,dmr_dss,dmr_dss_unpaired,dmr_dss_unsmoothed,downsampling,annotation_enrichment,dmr_density,coverage_context,imprinting,logs}

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

# ---- Merged matrices and cpg_stats_methylbench.tab from the per-platform calls
if [[ $FROM_CALLS -eq 1 ]]; then
  run_step matrices \
    Rscript scripts/R/07b_build_matrices.R \
      --samplesheet "$SAMPLESHEET" \
      --calls_sheet "$CALLS_SHEET" \
      --epic_path   "$EPIC_TSV" \
      --outdir      "$DATA/"

  run_step cpg_stats \
    Rscript scripts/R/07_generate_cpg_stats.R \
      --samplesheet "$SAMPLESHEET" \
      --calls_sheet "$CALLS_SHEET" \
      --outdir      "$STATS/"
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

# ---- Figure 3, Suppl. Figure 1, Suppl. Table S4 - QC and coverage uniformity -
run_step qc \
  Rscript scripts/R/08_qc_visualization.R \
    --cpg_stats "$STATS/cpg_stats_methylbench.tab" \
    --qc_stats  "$STATS/qc_stats_methylbench.tab" \
    --all_path  "$DATA/ALL_without_EPIC.csv" \
    --outdir    "$OUT/figures/"

# ---- Figure 4, Suppl. Figures 2 and 3, Suppl. Table S5 - correlation --------
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

# ---- Suppl. Table S7 - unpaired sensitivity run (tables only) ---------------
run_step unpaired \
  Rscript scripts/R/13_DMR_DSS_analysis.R \
    --seq_path "$DATA/ALL_without_EPIC.csv" \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/figures_unpaired/" \
    --datadir  "$OUT/dmr_dss_unpaired/" \
    --unpaired

# ---- Suppl. Table S7 - sensitivity run without DSS smoothing (tables only) --
run_step unsmoothed \
  Rscript scripts/R/13_DMR_DSS_analysis.R \
    --seq_path "$DATA/ALL_without_EPIC.csv" \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/figures_unsmoothed/" \
    --datadir  "$OUT/dmr_dss_unsmoothed/" \
    --no_smoothing

run_step sensitivity_summary \
  python3 scripts/python/compare_dss_runs.py \
    --run primary="$OUT/dmr_dss" --run unsmoothed="$OUT/dmr_dss_unsmoothed" \
    --run unpaired="$OUT/dmr_dss_unpaired" -o "$OUT/dss_sensitivity_summary.tsv"

# ---- Suppl. Figure 9, Suppl. Table S8 - depth-matched sensitivity -----------
# needs BSseq_Tier1.rds, Tier1_DML.rds and Tier1_DMR.rds from step 'dss'
run_step downsampling \
  Rscript scripts/R/14_downsampling_sensitivity.R \
    --dss_dir "$OUT/dmr_dss/" \
    --outdir  "$OUT/figures/" \
    --datadir "$OUT/downsampling/" \
    --n_reps  "$N_REPS" \
    --cores   "$CORES"

# ---- Figure 9C, Suppl. Figure 11, Suppl. Tables S9/S10 - annotation ---------
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

# ---- Suppl. Figures 1E,F and 12E, Suppl. Tables S4B and S11F ----------------
# coverage by genomic context and DMC/DMR sensitivity (needs step 'dss')
run_step coverage_context \
  Rscript scripts/R/18_coverage_by_context.R \
    --all_path "$DATA/ALL_without_EPIC.csv" \
    --dss_dir  "$OUT/dmr_dss/" \
    --outdir   "$OUT/figures/" \
    --datadir  "$OUT/coverage_context/"

# ---- Suppl. Figure 4, Suppl. Table S6 - imprinted gDMRs --------------------
# reads the per-platform calls (CALLS_SHEET) instead of the matrices; uses the
# DML results from step 'dss' for the DMCs located in gDMRs (Suppl. Table S6E).
GDMR_BED="${GDMR_BED:-data/imprinted_gDMRs_hg38.bed}"
if [[ $FROM_CALLS -eq 1 && -f "$GDMR_BED" ]]; then
  run_step imprinting \
    Rscript scripts/R/17_imprinted_regions.R \
      --samplesheet "$SAMPLESHEET" \
      --calls_sheet "$CALLS_SHEET" \
      --epic_path   "$EPIC_TSV" \
      --gdmr_bed    "$GDMR_BED" \
      --dss_dir     "$OUT/dmr_dss/" \
      --outdir      "$OUT/figures/" \
      --datadir     "$OUT/imprinting/"
else
  echo "[skip] imprinting (needs CALLS_SHEET, SAMPLESHEET, EPIC_TSV and GDMR_BED)"
fi

echo "[$(date '+%H:%M:%S')] MethylBench analysis finished. Results: $OUT"
