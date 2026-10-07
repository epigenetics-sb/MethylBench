#!/usr/bin/env bash
# =============================================================================
# MethylBench - run the complete R analysis (scripts 08-15) in one go
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

mkdir -p "$OUT"/{figures,figures_unpaired,diff_meth,dmr_dss,dmr_dss_unpaired,downsampling,annotation_enrichment,logs}

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

# ---- Figure 3, Suppl. Figure 16, Suppl. Table S6 - QC and coverage uniformity
run_step qc \
  Rscript scripts/R/08_qc_visualization.R \
    --cpg_stats "$STATS/cpg_stats_methylbench.tab" \
    --qc_stats  "$STATS/qc_stats_methylbench.tab" \
    --all_path  "$DATA/ALL_without_EPIC.csv" \
    --outdir    "$OUT/figures/"

# ---- Figure 4, Suppl. Figures 2 and 19, Suppl. Table S7 - correlation -------
run_step correlation \
  Rscript scripts/R/09_correlation_analysis.R \
    --datadir "$DATA/" \
    --outdir  "$OUT/figures/"

# ---- Figure 5, Suppl. Figures 3 and 4 - methylation density -----------------
run_step density \
  Rscript scripts/R/10_density_plots.R \
    --datadir   "$DATA/" \
    --outdir    "$OUT/figures/" \
    --epic_path "$DATA/ALL.csv"

# ---- Figure 6, Suppl. Figure 1 - PCA ----------------------------------------
run_step pca \
  Rscript scripts/R/11_pca.R \
    --all_path   "$DATA/ALL.csv" \
    --blood_path "$DATA/Blood_without_EPIC.csv" \
    --fibro_path "$DATA/Fibro_without_EPIC.csv" \
    --giab_path  "$DATA/GIAB_without_EPIC.csv" \
    --outdir     "$OUT/figures/"

# ---- Exploratory paired limma (input for script 12) -------------------------
run_step limma \
  Rscript scripts/R/limma_diff_meth.R \
    --all_path "$DATA/ALL.csv" \
    --outdir   "$OUT/diff_meth/"

# ---- Figure 7, Suppl. Figures 6-8 - exploratory limma / Wilcoxon ------------
run_step diff_meth \
  Rscript scripts/R/12_differential_methylation.R \
    --all_path   "$DATA/ALL.csv" \
    --blood_path "$DATA/Blood_without_EPIC.csv" \
    --fibro_path "$DATA/Fibro_without_EPIC.csv" \
    --limma_dir  "$OUT/diff_meth/" \
    --outdir     "$OUT/figures/" \
    --datadir    "$OUT/diff_meth/"

# ---- Figures 8, 9A/B/D, Suppl. Figures 12-15 - paired DSS / DMRcate ---------
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

# ---- Suppl. Figure 17, Suppl. Table S5 - depth-matched sensitivity ----------
# needs BSseq_Tier1.rds, Tier1_DML.rds and Tier1_DMR.rds from step 'dss'
run_step downsampling \
  Rscript scripts/R/14_downsampling_sensitivity.R \
    --dss_dir "$OUT/dmr_dss/" \
    --outdir  "$OUT/figures/" \
    --datadir "$OUT/downsampling/" \
    --n_reps  "$N_REPS" \
    --cores   "$CORES"

# ---- Figure 9C, Suppl. Figures 9 and 18, Suppl. Tables S3/S4 - annotation ---
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

echo "[$(date '+%H:%M:%S')] MethylBench analysis finished. Results: $OUT"
