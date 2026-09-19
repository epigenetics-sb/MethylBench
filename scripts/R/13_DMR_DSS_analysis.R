#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Differential Methylation Analysis, Tier 1/Tier 2 (Section 2.5.2
# and 2.5.3): Figure 8, Figure 9, Supplementary Figure 15
# =============================================================================
# Description:
#   Single source of truth for the primary (coverage-matched consensus set)
#   differential methylation analysis:
#     - Builds the Tier 1 (sequencing methods only) and Tier 2 (additionally
#       intersected with EPIC) consensus CpG sets and the corresponding
#       BSseq objects.
#     - Subject-paired single-CpG (DML) testing for the four sequencing
#       methods (DSS Beta-Binomial via DMLfit.multiFactor/DMLtest.multiFactor)
#       and for EPIC (limma on M-values), both with a `~ subject + group`
#       design.
#     - Region-level DMRs via DMRcate, fed directly from the paired per-CpG
#       test statistics above.
#     - Figure 8 (CpG-level concordance): pairwise Delta-beta scatter, UpSet,
#       exclusivity, hyper/hypo barplot, Jaccard heatmap, coverage
#       distribution on the consensus set.
#     - Figure 9 (region-level concordance): DMR counts, width, pairwise
#       Jaccard, CpG-structural context of DMRs.
#     - Supplementary Figure 15: threshold-free rank-recovery (ROC/AUC of
#       each sequencing platform's continuous DMC score against every other
#       platform's called DMCs as reference).
#
# CHANGELOG (post-review consolidation, Reviewer 1 / Major Comments 2 & 3):
#   This script previously only covered the DSS side of the analysis and
#   fed an internal, unpaired DSS::callDMR() side-computation (its own
#   dmr_DSS_counts/width/Jaccard figures). A separate, hardcoded-path,
#   not-yet-reviewed script turned out to be the actual source of Figure 8,
#   Figure 9 (via DMRcate) and the primary-stage EPIC-on-M-values limma
#   analysis, independently rebuilding the SAME Tier 1/Tier 2 consensus and
#   BSseq objects a second time. Both have now been merged into this single
#   script, chronologically ordered to match how Section 2.5 of the
#   manuscript reads, so there is exactly one place that defines the Tier
#   1/Tier 2 consensus, one paired DML/DMC computation per platform, and
#   one script producing everything derived from it (Figures 8 and 9,
#   Supplementary Figure 15). 15_DMR_DMRcate_analysis.R independently
#   recomputed the same consensus/BSseq construction and an UNPAIRED
#   DSS::DMLtest() + DMRcate a third time (writing to the same
#   BSseq_Tier{1,2}.rds / Tier{1,2}_DML_significant_<method>.tsv filenames
#   as this script, which would silently overwrite one another if both
#   were ever run into the same --datadir) and has been removed.
#
#   1) FIXED missing subject-level pairing:
#      - Sequencing methods: single-CpG testing now uses
#        DSS::DMLfit.multiFactor()/DMLtest.multiFactor() with a
#        `~ subject + group` design (run_dss_tissue_paired()), subject IDs
#        parsed from the BSseq sample names (get_subject_id_dss()), not
#        assumed from column position -- same principle as the PCA
#        label-desync fix in 11_pca.R and the limma fix in
#        limma_diff_meth.R.
#      - EPIC (Tier 2 only, no read counts available for a Beta-Binomial
#        model): differential methylation is assessed with limma on
#        M-values (beta2m(), clipped near 0/1), design `~ subject + group`,
#        exactly mirroring limma_diff_meth.R's fix. This is the primary-
#        analysis EPIC step described in the Methods ("EPIC was analyzed
#        using limma on M-values"), as distinct from the exploratory,
#        beta-scale, all-platform comparison in limma_diff_meth.R /
#        12_differential_methylation.R (Section 2.5.1).
#   2) REMOVED the old unpaired DSS::DMLtest()+callDMR() side-analysis and
#      its dmr_DSS_counts/width/Jaccard outputs. DSS::callDMR() only
#      accepts the two-group DMLtest() object (it needs the smoothed
#      mu1/mu2/diff/areastat fields that only that function produces) and
#      cannot be given the paired DMLtest.multiFactor() result at all, so
#      it could only ever be run unpaired -- and its outputs did not
#      correspond to any figure or supplementary figure in the manuscript.
#      Region-level DMRs are now called exclusively via DMRcate (see
#      below), which has no such limitation: it consumes a generic per-CpG
#      stat/p-value/effect-size table regardless of which model produced
#      it, so feeding it the paired test directly makes Figure 9 paired
#      with no workaround needed. Do not describe region-level DMR calling
#      in the manuscript as unpaired -- only DSS's own two-group DMR
#      caller (no longer used here) had that limitation.
#   3) FIXED WGBS/WGEC naming throughout the merged-in sections (the
#      donor script pre-dated the repository-wide rename); a compatibility
#      shim renames legacy "WGBS_*" input columns to "WGEC_*" if present.
#   4) REMOVED the CpG-/genic-feature annotation composition that the donor
#      script also computed (background-free, partial duplicate of
#      15_annotation_enrichment_background.R); 15 is the single,
#      background-corrected source for that analysis. ADDED a CpG-
#      structural-only context annotation of the called DMRs (Figure 9C),
#      which was not implemented anywhere else in the repository.
#   5) REMOVED a broken Tier-2 "Sequencing_withEPIC" Delta-beta scatter
#      call from the donor script (its DML list never carried an "EPIC"
#      entry, so it silently plotted only the four sequencing methods
#      under a misleading label; it also did not correspond to any panel
#      in the manuscript -- Figure 8A is Tier 1 only, six platform pairs;
#      EPIC vs. TWIST is Figure 7B, from 12_differential_methylation.R).
#   6) Hardcoded absolute paths in the donor script replaced with the
#      --outdir/--datadir arguments already used throughout this script.
#
#   ACTION ITEM (not addressed here): 14_downsampling_sensitivity.R keeps
#   its own local, unpaired run_dss_tissue()/call_dss_dmr() copies as its
#   "before" baseline (it does not import from this script). For full
#   consistency with the now fully paired Tier 1 analysis here, it should
#   be updated to depth-match against the paired model instead. Left as-is
#   for now since it directly answers a separate reviewer comment on its
#   own terms.
#
#   RENUMBERING: this script was previously 14_DMR_DSS_analysis.R. The
#   exploratory-stage, single-hierarchy, no-background 13_annotation.R
#   (Reviewer 1's Major Comment 3) has been retired -- its per-method BED
#   output was consumed only by 17_annotation_enrichment_background.R,
#   which has been rewritten to read Tier 1/Tier 2 DMC/consensus files
#   directly from this script's --datadir instead (fixing an additional,
#   previously unnoticed scope mismatch: it had compared exploratory-stage
#   significant DMCs against a Tier 1/2 background). With 13_annotation.R
#   gone, every subsequent script shifts down by one number: this script
#   14 -> 13, 16_downsampling_sensitivity.R -> 14, and
#   17_annotation_enrichment_background.R -> 15.
#
# Input:
#   --all_path      Path to the combined methylation matrix
#                   (EPIC + sequencing methods).
#
#   Column naming convention:
#     EPIC_Blood1..n, EPIC_Fibro1..n
#     ONT_Blood1..n,  ONT_Fibro1..n
#     TWIST_Blood1..n, TWIST_Fibro1..n
#     WGEC_Blood1..n, WGEC_Fibro1..n
#     RRBS_Blood1..n, RRBS_Fibro1..n
#     ONT_cov_Blood1..n, ...
#     WGEC_cov_Blood1..n, ...
#
# Output:
#   --outdir        Directory for figures (Figure 8, Figure 9, Supplementary
#                   Figure 15 panels).
#   --datadir       Directory for DML/DMC/DMR tables and intermediate RDS
#                   files (BSseq_Tier{1,2}.rds, Tier{1,2}_consensus_CpGs.tsv,
#                   Tier{1,2}_DML_<method>.tsv [full paired test],
#                   Tier{1,2}_DML_significant_<method>.tsv [called DMCs],
#                   Tier2_EPIC_* [paired EPIC limma-on-M-values results]).
#
# Usage:
#   Rscript scripts/R/13_DMR_DSS_analysis.R \
#     --all_path    data/matrices/ALL_with_EPIC.csv \
#     --outdir      results/figures/ \
#     --datadir     results/dmr_dss/
#
# Author: MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(bsseq)
  library(BiocGenerics)
  library(GenomicRanges)
  library(DMRcate)
  library(limma)
  library(annotatr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(circlize)
  library(grid)
  library(scales)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path",
    type = "character",
    help = "Path to combined methylation matrix [required]",
    metavar = "FILE"
  ),
  make_option("--outdir",
    type    = "character",
    default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"
  ),
  make_option("--datadir",
    type    = "character",
    default = "results/dmr_dss/",
    help    = "Output directory for DML/DMC/DMR tables [default: results/dmr_dss/]",
    metavar = "DIR"
  ),
  make_option("--min_cov",
    type    = "double",
    default = 10,
    help    = "Minimum coverage per CpG [default: 10]",
    metavar = "FLOAT"
  ),
  make_option("--min_samples",
    type    = "integer",
    default = 4,
    help    = "Minimum number of samples passing coverage [default: 4]",
    metavar = "INT"
  ),
  make_option("--fdr_cutoff",
    type    = "double",
    default = 0.05,
    help    = "FDR threshold for significant DMCs [default: 0.05]",
    metavar = "FLOAT"
  ),
  make_option("--delta_cutoff",
    type    = "double",
    default = 0.1,
    help    = "Delta-beta threshold for significant DMCs [default: 0.1]",
    metavar = "FLOAT"
  ),
  make_option("--lambda",
    type    = "double",
    default = 1000,
    help    = "DMRcate smoothing bandwidth lambda, bp [default: 1000]",
    metavar = "FLOAT"
  ),
  make_option("--C",
    type    = "double",
    default = 2,
    help    = "DMRcate scaling parameter C [default: 2]",
    metavar = "FLOAT"
  ),
  make_option("--min_cpgs",
    type    = "integer",
    default = 3,
    help    = "Minimum number of CpGs per DMR [default: 3]",
    metavar = "INT"
  ),
  make_option("--genome",
    type    = "character",
    default = "hg38",
    help    = "Genome build [default: hg38]",
    metavar = "STRING"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$all_path)) stop("ERROR: --all_path is required")
if (!file.exists(opt$all_path)) stop(paste("File not found:", opt$all_path))

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

MIN_COV       <- opt$min_cov
MIN_SAMPLES   <- opt$min_samples
FDR_CUTOFF    <- opt$fdr_cutoff
DELTA_CUTOFF  <- opt$delta_cutoff
LAMBDA        <- opt$lambda
C_PARAM       <- opt$C
MIN_CPGS      <- opt$min_cpgs
GENOME        <- opt$genome

METHOD_PREFIX <- c(
  ONT   = "ONT",
  TWIST = "TWIST",
  WGEC  = "WGEC",
  RRBS  = "RRBS"
)

METHOD_COLORS <- c(
  ONT   = "#E69F00",
  TWIST = "#009E73",
  WGEC  = "#5654E9",
  RRBS  = "#0072B2"
)

METHOD_COLORS_EPIC <- c(METHOD_COLORS, EPIC = "#CC79A7")

TISSUES <- c("Blood", "Fibro")

cat("[1/8] Loading methylation matrix...\n")
combined_df <- fread(
  opt$all_path,
  header = TRUE,
  sep = ",",
  na.strings = "NA"
)

# --- WGBS -> WGEC compatibility shim ----------------------------------------
# Guards against an input matrix that still uses the legacy "WGBS_*"
# column naming from before the repository-wide WGEC rename.
wgbs_cols <- grep("^WGBS(_cov)?_(Blood|Fibro)[0-9]+$", colnames(combined_df), value = TRUE)
if (length(wgbs_cols) > 0) {
  cat(sprintf("Renaming %d legacy WGBS_* column(s) to WGEC_*...\n", length(wgbs_cols)))
  setnames(combined_df, wgbs_cols, sub("^WGBS", "WGEC", wgbs_cols))
}

required_coord_cols <- c("chr", "start")
missing_coord_cols <- setdiff(required_coord_cols, colnames(combined_df))
if (length(missing_coord_cols) > 0) {
  stop(
    "Required coordinate columns not found: ",
    paste(missing_coord_cols, collapse = ", ")
  )
}

combined_df[, cpg_id := paste0(chr, ":", start)]

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

get_meth_cols <- function(df, method, tissue) {
  prefix <- METHOD_PREFIX[[method]]
  grep(paste0("^", prefix, "_", tissue, "[0-9]+$"), colnames(df), value = TRUE)
}

get_cov_cols <- function(df, method, tissue) {
  prefix <- METHOD_PREFIX[[method]]
  grep(paste0("^", prefix, "_cov_", tissue, "[0-9]+$"), colnames(df), value = TRUE)
}

get_passing_cpgs_seq <- function(df, method, tissue,
                                 min_cov = MIN_COV, min_samples = MIN_SAMPLES) {
  cov_cols <- get_cov_cols(df, method, tissue)
  if (length(cov_cols) == 0) {
    stop(sprintf("[%s | %s] No coverage columns found.", method, tissue))
  }
  cov_mat <- as.matrix(df[, ..cov_cols])
  storage.mode(cov_mat) <- "numeric"
  pass <- rowSums(!is.na(cov_mat) & cov_mat >= min_cov) >= min_samples
  df$cpg_id[pass]
}

get_passing_cpgs_epic <- function(df, tissue, min_samples = MIN_SAMPLES) {
  beta_cols <- grep(paste0("^EPIC_", tissue, "[0-9]+$"), colnames(df), value = TRUE)
  if (length(beta_cols) == 0) {
    stop(sprintf("[EPIC | %s] No EPIC columns found.", tissue))
  }
  beta_mat <- as.matrix(df[, ..beta_cols])
  pass <- rowSums(!is.na(beta_mat)) >= min_samples
  df$cpg_id[pass]
}

make_bsseq_from_combined <- function(df, method, tissue) {
  beta_cols <- get_meth_cols(df, method, tissue)
  cov_cols  <- get_cov_cols(df, method, tissue)

  if (length(beta_cols) == 0) stop(sprintf("[%s | %s] No methylation columns found.", method, tissue))
  if (length(cov_cols) == 0) stop(sprintf("[%s | %s] No coverage columns found.", method, tissue))
  if (length(beta_cols) != length(cov_cols)) {
    stop(sprintf("[%s | %s] Number of methylation columns (%d) does not match number of coverage columns (%d).",
                 method, tissue, length(beta_cols), length(cov_cols)))
  }

  cov_mat <- as.matrix(df[, ..cov_cols]); storage.mode(cov_mat) <- "numeric"
  beta_mat <- as.matrix(df[, ..beta_cols]); storage.mode(beta_mat) <- "numeric"
  meth_mat <- round(beta_mat * cov_mat)

  na_mask <- is.na(cov_mat) | cov_mat <= 0
  cov_mat[na_mask] <- 0
  meth_mat[is.na(meth_mat) | na_mask] <- 0

  sample_names <- sub(paste0("^", METHOD_PREFIX[[method]], "_"), paste0(method, "_"), beta_cols)

  BSseq(
    chr         = as.character(df$chr),
    pos         = as.integer(df$start),
    M           = unname(meth_mat),
    Cov         = unname(cov_mat),
    sampleNames = sample_names
  )
}

#' Parse subject/individual ID directly from a sample/column name (e.g.
#' "ONT_Blood3" -> "3", "EPIC_Fibro3" -> "3"). Blood/Fibro samples sharing
#' a trailing number are assumed to come from the same individual (matched
#' sampling; see Methods). Deriving this from the actual names, rather than
#' assuming blood/fibro samples line up positionally, means pairing cannot
#' silently break if sample order changes upstream -- same principle as the
#' PCA label-desync fix in 11_pca.R and get_subject_id() in
#' limma_diff_meth.R.
get_subject_id <- function(sample_names) {
  m  <- regmatches(sample_names, regexec("(?:Blood|Fibro)([0-9]+)$", sample_names))
  ok <- lengths(m) == 2
  if (!all(ok)) {
    stop("get_subject_id(): could not parse subject ID from: ",
         paste(sample_names[!ok], collapse = ", "))
  }
  vapply(m, `[[`, character(1), 2)
}

#' Per-CpG mean-methylation difference (Blood - Fibro), computed directly
#' from the raw M/Cov counts of a combined BSseq object. This is the
#' effect-size companion to the paired DMLtest.multiFactor() p-values,
#' since that function tests significance but does not itself report a
#' delta-beta-like effect size.
compute_delta_beta_bsseq <- function(bs_combined, blood_samples, fibro_samples) {
  M    <- getCoverage(bs_combined, type = "M")
  Cov  <- getCoverage(bs_combined, type = "Cov")
  beta <- M / Cov
  beta[Cov == 0] <- NA
  rowMeans(beta[, blood_samples, drop = FALSE], na.rm = TRUE) -
    rowMeans(beta[, fibro_samples, drop = FALSE], na.rm = TRUE)
}

#' Subject-paired single-CpG DML test (design ~ subject + group).
#' DSS::DMLtest() only supports a plain two-group comparison with no
#' covariate argument; DSS::DMLfit.multiFactor()/DMLtest.multiFactor()
#' support an arbitrary design formula, so subject pairing is added here,
#' mirroring the `~ subject + group` fix in limma_diff_meth.R.
#'
#' NOTE: the returned test data frame (chr, pos, stat, pvals, fdrs) has no
#' smoothed mu1/mu2/diff/areastat fields, so it cannot be passed to
#' DSS::callDMR() -- see the CHANGELOG. Region-level DMRs are instead
#' called from these per-CpG statistics via DMRcate (dss_to_cpgannotated()/
#' run_dmrcate() below), which has no such restriction.
run_dss_tissue_paired <- function(bs_blood, bs_fibro) {
  blood_samples <- sampleNames(bs_blood)
  fibro_samples <- sampleNames(bs_fibro)
  bs_combined   <- BiocGenerics::combine(bs_blood, bs_fibro)

  subject_ids <- get_subject_id(sampleNames(bs_combined))

  design <- data.frame(
    subject = factor(subject_ids, levels = sort(unique(subject_ids))),
    group   = factor(
      c(rep("Blood", length(blood_samples)), rep("Fibro", length(fibro_samples))),
      levels = c("Fibro", "Blood")
    )
  )

  fit  <- DMLfit.multiFactor(bs_combined, design = design, formula = ~ subject + group)
  test <- DMLtest.multiFactor(fit, coef = "groupBlood")

  delta_beta <- compute_delta_beta_bsseq(bs_combined, blood_samples, fibro_samples)
  stopifnot(nrow(test) == length(delta_beta))
  test$delta_beta <- delta_beta

  test
}

#' Build a DMRcate CpGannotated object from a generic per-CpG test-statistic
#' table. DMRcate does not care which model produced stat/rawpval/diff --
#' unlike DSS::callDMR(), it has no requirement that the input come from an
#' unpaired two-group DMLtest(). This is what makes it possible to feed it
#' the paired sequencing test (chr/pos/stat/pvals/delta_beta/fdrs) and the
#' paired EPIC limma-on-M-values test on equal footing.
make_cpgannotated <- function(chr, pos, stat, rawpval, diff, fdr, ids = NULL) {
  gr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start = pos, width = 1),
    stat     = stat,
    rawpval  = rawpval,
    diff     = diff,
    ind.fdr  = fdr,
    is.sig   = !is.na(fdr) & fdr < FDR_CUTOFF & abs(diff) >= DELTA_CUTOFF
  )
  if (is.null(ids)) ids <- paste0(chr, ":", pos)
  names(gr) <- ids
  new("CpGannotated", ranges = gr)
}

dss_test_to_cpgannotated <- function(test_df) {
  pval_col <- intersect(c("pvals", "pval"), colnames(test_df))[1]
  make_cpgannotated(
    chr = test_df$chr, pos = test_df$pos, stat = test_df$stat,
    rawpval = test_df[[pval_col]], diff = test_df$delta_beta, fdr = test_df$fdrs
  )
}

run_dmrcate <- function(cpg_ann, method_label, lambda = LAMBDA, C = C_PARAM, min_cpgs = MIN_CPGS) {
  cat(sprintf("  [%s] DMRcate...\n", method_label))
  dmr <- tryCatch(
    dmrcate(cpg_ann, lambda = lambda, C = C, min.cpgs = min_cpgs),
    error = function(e) {
      cat(sprintf("    [%s] ERROR: %s\n", method_label, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dmr)) return(NULL)
  dmr_gr <- extractRanges(dmr, genome = GENOME)
  cat(sprintf("    [%s] %d DMRs\n", method_label, length(dmr_gr)))
  dmr_gr
}

jaccard_generic <- function(a, b) {
  if (is.null(a) || is.null(b) || length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_dmr <- function(gr1, gr2) {
  if (is.null(gr1) || is.null(gr2) || length(gr1) == 0 || length(gr2) == 0) return(NA_real_)
  gr1 <- GenomicRanges::reduce(gr1); gr2 <- GenomicRanges::reduce(gr2)
  inter <- GenomicRanges::intersect(gr1, gr2)
  union <- GenomicRanges::reduce(c(gr1, gr2))
  if (length(union) == 0) return(NA_real_)
  sum(width(inter)) / sum(width(union))
}

build_jaccard_mat <- function(item_list, jaccard_fun) {
  keys <- names(item_list)[vapply(item_list, function(x) !is.null(x) && length(x) > 0, logical(1))]
  if (length(keys) == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
  mat <- outer(keys, keys, Vectorize(function(i, j) jaccard_fun(item_list[[i]], item_list[[j]])))
  rownames(mat) <- colnames(mat) <- keys
  mat
}

write_jaccard <- function(mat, path) {
  if (length(mat) == 0) return(invisible(NULL))
  fwrite(as.data.table(mat, keep.rownames = "Method"), path, sep = "\t")
}

# -------------------------------------------------------------------------
# 2. Determine consensus CpG sets
# -------------------------------------------------------------------------

cat("[2/8] Determining consensus CpG sets...\n")

passing_cpgs <- list()

for (method in names(METHOD_PREFIX)) {
  blood_pass <- get_passing_cpgs_seq(combined_df, method, "Blood")
  fibro_pass <- get_passing_cpgs_seq(combined_df, method, "Fibro")
  both_pass  <- intersect(blood_pass, fibro_pass)
  cat(sprintf("  [%s] Blood: %d | Fibro: %d | Blood ∩ Fibro: %d\n",
              method, length(blood_pass), length(fibro_pass), length(both_pass)))
  passing_cpgs[[method]] <- both_pass
}

epic_blood_pass <- get_passing_cpgs_epic(combined_df, "Blood")
epic_fibro_pass <- get_passing_cpgs_epic(combined_df, "Fibro")
epic_both_pass  <- intersect(epic_blood_pass, epic_fibro_pass)
cat(sprintf("  [EPIC] Blood: %d | Fibro: %d | Blood ∩ Fibro: %d\n",
            length(epic_blood_pass), length(epic_fibro_pass), length(epic_both_pass)))
passing_cpgs[["EPIC"]] <- epic_both_pass

consensus_tier1 <- Reduce(intersect, passing_cpgs[names(METHOD_PREFIX)])
consensus_tier2 <- intersect(consensus_tier1, passing_cpgs[["EPIC"]])

cat(sprintf("  Tier 1 consensus: %d CpGs\n", length(consensus_tier1)))
cat(sprintf("  Tier 2 consensus: %d CpGs\n", length(consensus_tier2)))

combined_tier1 <- combined_df[combined_df$cpg_id %in% consensus_tier1]
combined_tier2 <- combined_df[combined_df$cpg_id %in% consensus_tier2]

stopifnot(
  nrow(combined_tier1) == length(consensus_tier1),
  nrow(combined_tier2) == length(consensus_tier2)
)

fwrite(data.table(cpg_id = consensus_tier1), file.path(opt$datadir, "Tier1_consensus_CpGs.tsv"), sep = "\t")
fwrite(data.table(cpg_id = consensus_tier2), file.path(opt$datadir, "Tier2_consensus_CpGs.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# 3. Build BSseq objects
# -------------------------------------------------------------------------

cat("[3/8] Building BSseq objects...\n")

bsseq_t1 <- list()
bsseq_t2 <- list()

for (method in names(METHOD_PREFIX)) {
  for (tissue in TISSUES) {
    key <- paste0(method, "_", tissue)
    bsseq_t1[[key]] <- make_bsseq_from_combined(combined_tier1, method, tissue)
    bsseq_t2[[key]] <- make_bsseq_from_combined(combined_tier2, method, tissue)
    cat(sprintf("  [%s] Tier1: %d CpGs | Tier2: %d CpGs | samples: %d\n",
                key, nrow(bsseq_t1[[key]]), nrow(bsseq_t2[[key]]), ncol(bsseq_t1[[key]])))
  }
}

saveRDS(bsseq_t1, file.path(opt$datadir, "BSseq_Tier1.rds"))
saveRDS(bsseq_t2, file.path(opt$datadir, "BSseq_Tier2.rds"))

# -------------------------------------------------------------------------
# 4. Paired DML/DMC + DMRcate per sequencing method (Tier 1 and Tier 2)
# -------------------------------------------------------------------------

cat("[4/8] Running paired DML testing and DMRcate per sequencing method...\n")

run_tier <- function(bsseq_list, tier_label) {

  dml_list <- list()  # full paired test, all consensus CpGs
  dmc_list <- list()  # significant subset (chr/pos/diff/fdr)
  dmr_list <- list()  # DMRcate regions

  for (method in names(METHOD_PREFIX)) {
    key_blood <- paste0(method, "_Blood")
    key_fibro <- paste0(method, "_Fibro")

    cat(sprintf("  [%s | %s] DMLtest.multiFactor (paired: ~subject+group)...\n", tier_label, method))

    test_df <- run_dss_tissue_paired(bsseq_list[[key_blood]], bsseq_list[[key_fibro]])
    dml_list[[method]] <- test_df

    dmc <- test_df[!is.na(test_df$fdrs) & test_df$fdrs < FDR_CUTOFF &
                     abs(test_df$delta_beta) > DELTA_CUTOFF, , drop = FALSE]
    dmc_list[[method]] <- data.frame(
      chr = dmc$chr, pos = dmc$pos, diff = dmc$delta_beta, fdr = dmc$fdrs
    )

    cat(sprintf("    [%s | %s] %d significant DMCs\n", tier_label, method, nrow(dmc)))

    fwrite(as.data.table(test_df),
           file.path(opt$datadir, paste0(tier_label, "_DML_", method, ".tsv")), sep = "\t")
    fwrite(as.data.table(dmc_list[[method]]),
           file.path(opt$datadir, paste0(tier_label, "_DML_significant_", method, ".tsv")), sep = "\t")

    cat(sprintf("  [%s | %s] DMRcate...\n", tier_label, method))
    dmr <- run_dmrcate(dss_test_to_cpgannotated(test_df), method)
    dmr_list[[method]] <- dmr
    if (!is.null(dmr) && length(dmr) > 0) {
      fwrite(as.data.table(as.data.frame(dmr)),
             file.path(opt$datadir, paste0(tier_label, "_DMR_", method, ".tsv")), sep = "\t")
    }
  }

  saveRDS(dml_list, file.path(opt$datadir, paste0(tier_label, "_DML.rds")))
  saveRDS(dmc_list, file.path(opt$datadir, paste0(tier_label, "_DMC.rds")))
  saveRDS(dmr_list, file.path(opt$datadir, paste0(tier_label, "_DMR.rds")))

  list(dml = dml_list, dmc = dmc_list, dmr = dmr_list)
}

seq_t1 <- run_tier(bsseq_t1, "Tier1")
seq_t2 <- run_tier(bsseq_t2, "Tier2")

# -------------------------------------------------------------------------
# 5. EPIC: paired limma on M-values (Tier 2 only)
# -------------------------------------------------------------------------
# EPIC provides no methylated/unmethylated read counts, so it cannot enter
# the DSS Beta-Binomial pipeline above; per the Methods, its primary-stage
# differential methylation is assessed with limma on M-values instead.
# Design is `~ subject + group`, mirroring the sequencing methods' pairing
# and limma_diff_meth.R's exploratory-stage fix.

cat("[5/8] Running paired EPIC limma-on-M-values (Tier 2)...\n")

epic_blood_cols <- grep("^EPIC_Blood[0-9]+$", colnames(combined_tier2), value = TRUE)
epic_fibro_cols <- grep("^EPIC_Fibro[0-9]+$", colnames(combined_tier2), value = TRUE)

epic_beta <- as.matrix(combined_tier2[, c(epic_blood_cols, epic_fibro_cols), with = FALSE])
rownames(epic_beta) <- combined_tier2$cpg_id
epic_beta <- epic_beta[complete.cases(epic_beta), , drop = FALSE]

beta_to_m <- function(beta, offset = 0.001) {
  beta <- pmin(pmax(beta, offset), 1 - offset)
  log2(beta / (1 - beta))
}
epic_m <- beta_to_m(epic_beta)

epic_subject <- factor(get_subject_id(c(epic_blood_cols, epic_fibro_cols)))
epic_group   <- factor(
  c(rep("Blood", length(epic_blood_cols)), rep("Fibro", length(epic_fibro_cols))),
  levels = c("Fibro", "Blood")
)
epic_design <- model.matrix(~ epic_subject + epic_group)
epic_fit    <- eBayes(lmFit(epic_m, epic_design))

epic_delta_beta <- rowMeans(epic_beta[, epic_blood_cols, drop = FALSE]) -
  rowMeans(epic_beta[, epic_fibro_cols, drop = FALSE])

epic_test <- topTable(epic_fit, coef = "epic_groupBlood", number = Inf,
                       adjust.method = "BH", sort.by = "none")
epic_test$chr        <- sub(":.*", "", rownames(epic_test))
epic_test$pos        <- as.integer(sub(".*:", "", rownames(epic_test)))
epic_test$delta_beta <- epic_delta_beta[rownames(epic_test)]
# Harmonize column naming with the DSS paired test output (chr/pos/stat/
# pvals/fdrs/delta_beta) so all five Tier-2 methods share one schema in
# Tier2_DML.rds, rather than mixing limma's P.Value/adj.P.Val naming in
# for just one of the five list entries.
epic_test$stat  <- epic_test$t
epic_test$pvals <- epic_test$P.Value
epic_test$fdrs  <- epic_test$adj.P.Val

epic_dmc <- epic_test[!is.na(epic_test$adj.P.Val) & epic_test$adj.P.Val < FDR_CUTOFF &
                        abs(epic_test$delta_beta) >= DELTA_CUTOFF, , drop = FALSE]

cat(sprintf("  [EPIC] %d significant DMCs\n", nrow(epic_dmc)))

fwrite(as.data.table(epic_test), file.path(opt$datadir, "Tier2_DML_EPIC.tsv"), sep = "\t")
fwrite(data.table(chr = epic_dmc$chr, pos = epic_dmc$pos,
                   diff = epic_dmc$delta_beta, fdr = epic_dmc$adj.P.Val),
       file.path(opt$datadir, "Tier2_DML_significant_EPIC.tsv"), sep = "\t")

cat("  [EPIC] DMRcate...\n")
epic_cpg_ann <- make_cpgannotated(
  chr = epic_test$chr, pos = epic_test$pos, stat = epic_test$t,
  rawpval = epic_test$P.Value, diff = epic_test$delta_beta, fdr = epic_test$adj.P.Val,
  ids = rownames(epic_test)
)
epic_dmr <- run_dmrcate(epic_cpg_ann, "EPIC")
if (!is.null(epic_dmr) && length(epic_dmr) > 0) {
  fwrite(as.data.table(as.data.frame(epic_dmr)),
         file.path(opt$datadir, "Tier2_DMR_EPIC.tsv"), sep = "\t")
}

# Fold EPIC into the Tier 2 lists so downstream Fig. 8/9 code can treat all
# five Tier-2 platforms uniformly.
seq_t2$dml[["EPIC"]] <- epic_test
seq_t2$dmc[["EPIC"]] <- data.frame(chr = epic_dmc$chr, pos = epic_dmc$pos,
                                    diff = epic_dmc$delta_beta, fdr = epic_dmc$adj.P.Val)
seq_t2$dmr[["EPIC"]] <- epic_dmr

# -------------------------------------------------------------------------
# 6. Figure 8 -- CpG-level (DMC) concordance
# -------------------------------------------------------------------------

cat("[6/8] Generating Figure 8 (CpG-level DMC concordance)...\n")

make_cpg_ids <- function(dmc_list) lapply(dmc_list, function(df) paste0(df$chr, ":", df$pos))

dmc_cpg_ids_t1 <- make_cpg_ids(seq_t1$dmc)
dmc_cpg_ids_t2 <- make_cpg_ids(seq_t2$dmc)

# 6a) Hyper/hypo count barplot ------------------------------------------

plot_dmc_counts <- function(dmc_list, tier_label, method_col) {
  count_df <- bind_rows(lapply(names(dmc_list), function(m) {
    df <- dmc_list[[m]]
    if (!m %in% names(method_col)) return(NULL)
    data.frame(method = m, n_hyper = sum(df$diff > 0, na.rm = TRUE),
               n_hypo = sum(df$diff < 0, na.rm = TRUE))
  })) |>
    tidyr::pivot_longer(c(n_hyper, n_hypo), names_to = "direction", values_to = "count") |>
    dplyr::mutate(
      direction = factor(direction, levels = c("n_hyper", "n_hypo"),
                          labels = c("Hypermethylated (Blood > Fibro)", "Hypomethylated (Blood < Fibro)")),
      method = factor(method, levels = names(method_col))
    )

  ggplot(count_df, aes(x = method, y = count, fill = direction)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("#d73027", "#4575b4")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste("DMC Counts –", tier_label), x = NULL, y = "#DMCs", fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 22),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22),
          text = element_text(size = 22),
          axis.text.x = element_text(color = method_col[levels(count_df$method)])) +
    guides(fill = guide_legend(nrow = 2))
}

ggsave(file.path(opt$outdir, "dmc_counts_Tier1.png"),
       plot_dmc_counts(seq_t1$dmc, "Tier1", METHOD_COLORS), width = 9, height = 7, dpi = 300)
ggsave(file.path(opt$outdir, "dmc_counts_Tier2.png"),
       plot_dmc_counts(seq_t2$dmc, "Tier2", METHOD_COLORS_EPIC), width = 9, height = 7, dpi = 300)

# 6b) Exclusivity -----------------------------------------------------------

plot_dmc_exclusivity <- function(dmc_cpg_ids, tier_label, method_col) {
  keys <- names(dmc_cpg_ids)
  all_cpgs <- unique(unlist(dmc_cpg_ids))
  membership <- sapply(dmc_cpg_ids, function(ids) all_cpgs %in% ids)

  excl_df <- bind_rows(lapply(keys, function(m) {
    my_idx <- which(membership[, m])
    n <- length(my_idx)
    if (n == 0) return(data.frame(method = m, Consensus = 0, Partial = 0, Exclusive = 0))
    rs <- rowSums(membership[my_idx, , drop = FALSE])
    data.frame(method = m,
               Consensus = sum(rs == length(keys)) / n,
               Partial   = sum(rs %in% 2:(length(keys) - 1)) / n,
               Exclusive = sum(rs == 1) / n)
  })) |>
    tidyr::pivot_longer(-method, names_to = "category", values_to = "fraction") |>
    dplyr::mutate(method = factor(method, levels = names(method_col)),
                  category = factor(category, levels = c("Consensus", "Partial", "Exclusive")))

  ggplot(excl_df, aes(x = method, y = fraction, fill = category)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(Consensus = "#1a9641", Partial = "#fdae61", Exclusive = "#d7191c")) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste("DMC Exclusivity –", tier_label), x = NULL, y = "Fraction of DMCs [%]", fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 22),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22),
          text = element_text(size = 22),
          axis.text.x = element_text(color = method_col[levels(excl_df$method)]))
}

ggsave(file.path(opt$outdir, "dmc_exclusivity_Tier1.png"),
       plot_dmc_exclusivity(dmc_cpg_ids_t1, "Tier1", METHOD_COLORS), width = 9, height = 7, dpi = 300)
ggsave(file.path(opt$outdir, "dmc_exclusivity_Tier2.png"),
       plot_dmc_exclusivity(dmc_cpg_ids_t2, "Tier2", METHOD_COLORS_EPIC), width = 9, height = 7, dpi = 300)

# 6c) DMC-level Jaccard heatmap ----------------------------------------------

plot_jaccard_heatmap <- function(jmat, tier_label, method_col, level_label) {
  keys <- intersect(rownames(jmat), names(method_col))
  jmat <- jmat[keys, keys, drop = FALSE]
  col_fun <- colorRamp2(c(0, 0.5, 1), c("#f7fbff", "#6baed6", "#08306b"))
  label_colors <- method_col[keys]

  top_ann  <- HeatmapAnnotation(Method = keys, col = list(Method = setNames(method_col[keys], keys)),
                                 show_annotation_name = FALSE, show_legend = FALSE)
  left_ann <- rowAnnotation(Method = keys, col = list(Method = setNames(method_col[keys], keys)),
                             show_annotation_name = FALSE, show_legend = FALSE)

  ht <- Heatmap(
    jmat, name = "Jaccard\nIndex", col = col_fun, na_col = "grey80",
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- jmat[i, j]
      if (!is.na(val)) grid.text(sprintf("%.3f", val), x, y,
                                  gp = gpar(fontsize = 11, col = ifelse(val > 0.5, "white", "black")))
    },
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 12, fontface = "bold", col = label_colors),
    column_names_gp = gpar(fontsize = 12, fontface = "bold", col = label_colors),
    column_title = paste0("Pairwise ", level_label, " Jaccard Index – ", tier_label),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    left_annotation = left_ann, top_annotation = top_ann
  )

  png(file.path(opt$outdir, paste0("jaccard_", tolower(level_label), "_", tier_label, ".png")),
      units = "in", width = 9, height = 7, res = 500)
  draw(ht)
  dev.off()
}

jaccard_dmc_t1 <- build_jaccard_mat(dmc_cpg_ids_t1, jaccard_generic)
jaccard_dmc_t2 <- build_jaccard_mat(dmc_cpg_ids_t2, jaccard_generic)
write_jaccard(jaccard_dmc_t1, file.path(opt$datadir, "DMC_Jaccard_Tier1.tsv"))
write_jaccard(jaccard_dmc_t2, file.path(opt$datadir, "DMC_Jaccard_Tier2.tsv"))
plot_jaccard_heatmap(jaccard_dmc_t1, "Tier1", METHOD_COLORS, "DMC")
plot_jaccard_heatmap(jaccard_dmc_t2, "Tier2", METHOD_COLORS_EPIC, "DMC")

# 6d) UpSet -------------------------------------------------------------

plot_dmc_upset <- function(dmc_cpg_ids, tier_label, method_col) {
  all_cpgs <- unique(unlist(dmc_cpg_ids))
  upset_data <- as.data.frame(sapply(dmc_cpg_ids, function(ids) all_cpgs %in% ids))
  keys <- names(method_col)[names(method_col) %in% colnames(upset_data)]

  assay_type <- c(ONT = "LongRead", TWIST = "ShortRead", WGEC = "ShortRead",
                   RRBS = "ShortRead", EPIC = "Array")
  stripes <- data.frame(set = keys, Assay = assay_type[keys])

  p <- ComplexUpset::upset(
    upset_data, intersect = keys,
    queries = lapply(keys, function(k) upset_query(set = k, fill = method_col[[k]])),
    set_sizes = (
      upset_set_size(geom = geom_bar(width = 0.8)) + ylab("DMCs") +
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        theme(text = element_text(size = 18))
    ),
    base_annotations = list(
      "Overlapping DMCs" = (
        ComplexUpset::intersection_size(width = 0.8, counts = FALSE, mapping = aes(fill = "bar")) +
          scale_fill_manual(values = c(bar = "#555555"), guide = "none") +
          scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
          theme(plot.background = element_rect(fill = "grey92", color = NA), text = element_text(size = 14))
      )
    ),
    stripes = upset_stripes(
      mapping = aes(color = Assay),
      colors = c(LongRead = "#FFE6CC", ShortRead = "#E1D5E7", Array = "#D5E8D4"),
      data = stripes
    ),
    name = "Method specific overlap", min_size = 10
  ) +
    theme(text = element_text(size = 14), axis.text = element_text(size = 12),
          axis.title = element_text(size = 13), strip.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 13)) +
    guides(fill = guide_legend(title = "Method")) +
    patchwork::plot_annotation(title = paste("Significant DMC overlap –", tier_label))

  ggsave(file.path(opt$outdir, paste0("upset_dmc_", tier_label, ".png")), p, width = 14, height = 7, dpi = 300)
}

plot_dmc_upset(dmc_cpg_ids_t1, "Tier1", METHOD_COLORS)
plot_dmc_upset(dmc_cpg_ids_t2, "Tier2", METHOD_COLORS_EPIC)

# 6e) Pairwise Delta-beta scatter (Tier 1 only -- matches Fig. 8A, six
#     sequencing-platform pairs; EPIC vs. TWIST is Fig. 7B, from
#     12_differential_methylation.R, not reproduced here) ------------------

plot_delta_beta_scatter <- function(dml_list, tier_label, method_col) {
  methods <- names(dml_list)
  pairs   <- combn(methods, 2, simplify = FALSE)

  plot_list <- lapply(pairs, function(pair) {
    m1 <- pair[1]; m2 <- pair[2]
    df1 <- data.frame(cpg_id = paste0(dml_list[[m1]]$chr, ":", dml_list[[m1]]$pos), diff1 = dml_list[[m1]]$delta_beta)
    df2 <- data.frame(cpg_id = paste0(dml_list[[m2]]$chr, ":", dml_list[[m2]]$pos), diff2 = dml_list[[m2]]$delta_beta)
    merged <- inner_join(df1, df2, by = "cpg_id")
    r_val <- round(cor(merged$diff1, merged$diff2, method = "pearson", use = "complete.obs"), 3)

    ggplot(merged, aes(x = diff1, y = diff2)) +
      geom_hex(bins = 80) +
      scale_fill_viridis_c(option = "magma", trans = "log10", name = "CpGs (log10)") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
      geom_smooth(method = "lm", se = FALSE, color = "white", linewidth = 0.6) +
      annotate("text", x = -0.8, y = 0.9, label = paste0("r = ", r_val),
               hjust = 0, size = 4.5, fontface = "bold", color = "white") +
      labs(title = paste(m1, "vs.", m2), x = paste0("Δβ ", m1), y = paste0("Δβ ", m2)) +
      coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            axis.text = element_text(size = 12), axis.title = element_text(size = 12), text = element_text(size = 12))
  })

  ncols <- ceiling(sqrt(length(plot_list)))
  nrows <- ceiling(length(plot_list) / ncols)

  p_combined <- patchwork::wrap_plots(plot_list, ncol = ncols) +
    patchwork::plot_annotation(
      title = paste("Pairwise Δβ Concordance –", tier_label),
      subtitle = "All consensus CpGs | red dashed = identity line | r = Pearson",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
                    plot.subtitle = element_text(hjust = 0.5, size = 16))
    )

  ggsave(file.path(opt$outdir, paste0("scatter_delta_beta_", tier_label, ".png")),
         p_combined, width = ncols * 4, height = nrows * 4, dpi = 300)
}

plot_delta_beta_scatter(seq_t1$dml, "Tier1", METHOD_COLORS)

# 6f) Coverage distribution on the consensus set -----------------------------

plot_coverage_distribution <- function(combined_consensus, method_col, min_cov = MIN_COV) {
  cov_df <- bind_rows(lapply(names(METHOD_PREFIX), function(m) {
    cov_cols <- get_cov_cols(combined_consensus, m, "Blood")
    cov_cols <- c(cov_cols, get_cov_cols(combined_consensus, m, "Fibro"))
    cov_mat  <- as.matrix(combined_consensus[, ..cov_cols])
    data.frame(method = m, sample = rep(cov_cols, each = nrow(cov_mat)), coverage = as.vector(cov_mat))
  })) |>
    dplyr::filter(!is.na(coverage)) |>
    dplyr::mutate(method = factor(method, levels = names(METHOD_PREFIX)),
                  tissue = ifelse(grepl("Blood", sample), "Blood", "Fibro"))

  p <- ggplot(cov_df, aes(x = method, y = coverage, fill = method)) +
    geom_violin(alpha = 0.8, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.9) +
    geom_hline(yintercept = min_cov, linetype = "dashed", color = "red", linewidth = 0.7) +
    annotate("text", x = 0.5, y = min_cov + 2, label = paste0("min cov = ", min_cov, "×"),
             hjust = 0, color = "red", size = 3.5) +
    scale_fill_manual(values = method_col[names(METHOD_PREFIX)]) +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~ tissue) +
    labs(title = "Coverage Distribution on Consensus CpG Set (Tier 1)", x = NULL, y = "Coverage (×, log10)") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          strip.background = element_rect(fill = "grey90", color = NA),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22), text = element_text(size = 22),
          axis.text.x = element_text(size = 20, colour = method_col[names(METHOD_PREFIX)]))

  ggsave(file.path(opt$outdir, "coverage_distribution_consensus_Tier1.png"), p, width = 9, height = 7, dpi = 300)
}

plot_coverage_distribution(combined_tier1, METHOD_COLORS)

# -------------------------------------------------------------------------
# 7. Figure 9 -- Region-level (DMR) concordance
# -------------------------------------------------------------------------

cat("[7/8] Generating Figure 9 (DMR-level concordance)...\n")

make_dmr_summary <- function(dmr_list, tier_label) {
  bind_rows(lapply(names(dmr_list), function(method) {
    dmr <- dmr_list[[method]]
    if (is.null(dmr) || length(dmr) == 0) {
      return(data.frame(Tier = tier_label, Method = method, DMRs = 0L,
                         MedianWidth = NA_real_, MedianCpGs = NA_real_))
    }
    data.frame(Tier = tier_label, Method = method, DMRs = length(dmr),
               MedianWidth = median(width(dmr), na.rm = TRUE),
               MedianCpGs  = median(dmr$no.cpgs, na.rm = TRUE))
  }))
}

dmr_summary <- bind_rows(make_dmr_summary(seq_t1$dmr, "Tier1"), make_dmr_summary(seq_t2$dmr, "Tier2"))
fwrite(as.data.table(dmr_summary), file.path(opt$datadir, "DMRcate_DMR_summary.tsv"), sep = "\t")

plot_dmr_counts <- function(dmr_list, tier_label, method_col) {
  count_df <- bind_rows(lapply(names(dmr_list), function(m) {
    dmr <- dmr_list[[m]]
    if (is.null(dmr) || length(dmr) == 0 || !m %in% names(method_col)) return(NULL)
    data.frame(method = m, n_hyper = sum(dmr$meandiff > 0, na.rm = TRUE),
               n_hypo = sum(dmr$meandiff < 0, na.rm = TRUE))
  })) |>
    tidyr::pivot_longer(c(n_hyper, n_hypo), names_to = "direction", values_to = "count") |>
    dplyr::mutate(
      direction = factor(direction, levels = c("n_hyper", "n_hypo"),
                          labels = c("Hypermethylated (Blood > Fibro)", "Hypomethylated (Blood < Fibro)")),
      method = factor(method, levels = names(method_col))
    )

  ggplot(count_df, aes(x = method, y = count, fill = direction)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("#d73027", "#4575b4")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste("DMR Counts per Method –", tier_label), x = NULL, y = "#DMRs", fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 22),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22), text = element_text(size = 22),
          axis.text.x = element_text(colour = method_col[levels(count_df$method)])) +
    guides(fill = guide_legend(nrow = 2))
}

ggsave(file.path(opt$outdir, "dmr_counts_Tier1.png"),
       plot_dmr_counts(seq_t1$dmr, "Tier1", METHOD_COLORS), width = 9, height = 7, dpi = 300)
ggsave(file.path(opt$outdir, "dmr_counts_Tier2.png"),
       plot_dmr_counts(seq_t2$dmr, "Tier2", METHOD_COLORS_EPIC), width = 9, height = 7, dpi = 300)

plot_dmr_width <- function(dmr_list, tier_label, method_col) {
  width_df <- bind_rows(lapply(names(dmr_list), function(m) {
    dmr <- dmr_list[[m]]
    if (is.null(dmr) || length(dmr) == 0 || !m %in% names(method_col)) return(NULL)
    data.frame(method = m, width_bp = width(dmr))
  })) |>
    dplyr::mutate(method = factor(method, levels = names(method_col)))

  ggplot(width_df, aes(x = method, y = width_bp, fill = method)) +
    geom_violin(alpha = 0.8, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = method_col) +
    scale_y_log10(labels = scales::comma) +
    labs(title = paste("DMR Width Distribution –", tier_label), x = NULL, y = "DMR width (bp, log10)") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 22),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22), text = element_text(size = 22),
          axis.text.x = element_text(colour = method_col[levels(width_df$method)]))
}

ggsave(file.path(opt$outdir, "dmr_width_Tier1.png"),
       plot_dmr_width(seq_t1$dmr, "Tier1", METHOD_COLORS), width = 9, height = 7, dpi = 300)
ggsave(file.path(opt$outdir, "dmr_width_Tier2.png"),
       plot_dmr_width(seq_t2$dmr, "Tier2", METHOD_COLORS_EPIC), width = 9, height = 7, dpi = 300)

jaccard_dmr_t1 <- build_jaccard_mat(seq_t1$dmr, jaccard_dmr)
jaccard_dmr_t2 <- build_jaccard_mat(seq_t2$dmr, jaccard_dmr)
write_jaccard(jaccard_dmr_t1, file.path(opt$datadir, "DMRcate_DMR_Jaccard_Tier1.tsv"))
write_jaccard(jaccard_dmr_t2, file.path(opt$datadir, "DMRcate_DMR_Jaccard_Tier2.tsv"))
plot_jaccard_heatmap(jaccard_dmr_t1, "Tier1", METHOD_COLORS, "DMR")
plot_jaccard_heatmap(jaccard_dmr_t2, "Tier2", METHOD_COLORS_EPIC, "DMR")

# Figure 9C: CpG-structural context of the called DMRs. This panel did not
# exist anywhere in the repository before this consolidation. It uses a
# single, CpG-structural-only hierarchy (Island > Shore > Shelf > Open Sea)
# -- deliberately NOT combined with gene-centric categories in one shared
# priority list, which is exactly the artifact behind Reviewer 1's Major
# Comment 3 on the now-retired, exploratory-stage 13_annotation.R (that
# script's number has been reassigned to this one; its own single-hierarchy,
# no-background analysis has been fully superseded by
# 15_annotation_enrichment_background.R and removed). No background
# comparison is made here either (the manuscript text for Fig. 9C only
# describes composition, not an enrichment test); for a background-
# corrected enrichment analysis at the DMR level, see
# 15_annotation_enrichment_background.R with --level DMR.
CPG_PRIORITY <- c(
  paste0(GENOME, "_cpg_islands"), paste0(GENOME, "_cpg_shores"),
  paste0(GENOME, "_cpg_shelves"), paste0(GENOME, "_cpg_inter")
)
cpg_db <- build_annotations(genome = GENOME, annotations = CPG_PRIORITY)

annotate_dmr_cpg_context <- function(dmr_list) {
  bind_rows(lapply(names(dmr_list), function(m) {
    dmr <- dmr_list[[m]]
    if (is.null(dmr) || length(dmr) == 0) return(NULL)
    ann <- annotate_regions(regions = dmr, annotations = cpg_db, ignore.strand = TRUE, quiet = TRUE)
    df <- as.data.frame(ann)
    df$priority <- match(df$annot.type, CPG_PRIORITY)
    df |>
      dplyr::group_by(seqnames, start, end) |>
      dplyr::slice_min(priority, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(method = m,
                    context = dplyr::case_when(
                      grepl("islands", annot.type) ~ "CpG Island",
                      grepl("shores",  annot.type) ~ "CpG Shore",
                      grepl("shelves", annot.type) ~ "CpG Shelf",
                      grepl("inter",   annot.type) ~ "Open Sea",
                      TRUE ~ "Other"
                    ))
  }))
}

plot_dmr_cpg_context <- function(dmr_list, tier_label, method_col) {
  ctx_df <- annotate_dmr_cpg_context(dmr_list)
  if (is.null(ctx_df) || nrow(ctx_df) == 0) return(invisible(NULL))

  plot_df <- ctx_df |>
    dplyr::filter(method %in% names(method_col)) |>
    dplyr::count(method, context) |>
    dplyr::group_by(method) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      context = factor(context, levels = c("CpG Island", "CpG Shore", "CpG Shelf", "Open Sea")),
      method  = factor(method, levels = names(method_col))
    )

  p <- ggplot(plot_df, aes(x = method, y = frac, fill = context)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("CpG Island" = "#2166ac", "CpG Shore" = "#74add1",
                                  "CpG Shelf" = "#abd9e9", "Open Sea" = "#e0f3f8")) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste("Genomic CpG Context of DMRs –", tier_label),
         x = NULL, y = "Fraction of DMRs", fill = "CpG Context") +
    theme_bw() +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22), text = element_text(size = 22),
          axis.text.x = element_text(color = method_col[levels(plot_df$method)]))

  ggsave(file.path(opt$outdir, paste0("dmr_cpg_context_", tier_label, ".png")), p, width = 9, height = 7, dpi = 300)
}

plot_dmr_cpg_context(seq_t1$dmr, "Tier1", METHOD_COLORS)
plot_dmr_cpg_context(seq_t2$dmr, "Tier2", METHOD_COLORS_EPIC)

# -------------------------------------------------------------------------
# 8. Supplementary Figure 15 -- threshold-free rank-recovery (ROC/AUC)
# -------------------------------------------------------------------------

cat("[8/8] Generating Supplementary Figure 15 (rank-recovery ROC/AUC)...\n")

plot_roc_vs_reference <- function(dml_list, reference, tier_label, method_col) {
  library(pROC)

  ref_df <- dml_list[[reference]]
  ref_df$cpg_id <- paste0(ref_df$chr, ":", ref_df$pos)
  ref_df$is_dmc <- !is.na(ref_df$fdrs) & ref_df$fdrs < FDR_CUTOFF & abs(ref_df$delta_beta) >= DELTA_CUTOFF

  methods_plot <- setdiff(names(dml_list), reference)

  roc_list <- bind_rows(lapply(methods_plot, function(m) {
    df <- dml_list[[m]]
    df$cpg_id <- paste0(df$chr, ":", df$pos)
    merged <- inner_join(ref_df[, c("cpg_id", "is_dmc")], df[, c("cpg_id", "fdrs", "delta_beta")], by = "cpg_id") |>
      dplyr::filter(!is.na(fdrs))
    merged$score <- -log10(merged$fdrs + 1e-300) * sign(merged$delta_beta)
    roc_obj <- roc(merged$is_dmc, merged$score, quiet = TRUE)
    data.frame(method = m, fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities,
               auc = as.numeric(auc(roc_obj)))
  }))

  auc_labels <- roc_list |> dplyr::distinct(method, auc) |>
    dplyr::mutate(label = paste0(method, " (AUC=", round(auc, 3), ")"))
  label_map <- setNames(auc_labels$label, auc_labels$method)

  p <- ggplot(roc_list, aes(x = fpr, y = tpr, color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = method_col[methods_plot], labels = label_map[methods_plot]) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste("DMC Recovery vs.", reference, "–", tier_label),
         subtitle = paste("Reference:", reference, "| Score = −log10(FDR) × sign(Δβ)"),
         x = "False Positive Rate (1 − Specificity)", y = "True Positive Rate (Sensitivity)", color = NULL) +
    theme_bw() +
    theme(legend.position = c(0.7, 0.3), plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          plot.subtitle = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size = 22), axis.title = element_text(size = 22), text = element_text(size = 22))

  ggsave(file.path(opt$outdir, paste0("roc_vs_", reference, "_", tier_label, ".png")),
         p, width = 9, height = 7, dpi = 300)
}

for (ref in names(METHOD_PREFIX)) {
  plot_roc_vs_reference(seq_t1$dml, reference = ref, tier_label = "Tier1", method_col = METHOD_COLORS)
}

# -------------------------------------------------------------------------
# Save session information
# -------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(opt$datadir, "sessionInfo.txt"))

cat("\nDone.\n")
cat("  DML/DMC/DMR tables : ", normalizePath(opt$datadir), "\n", sep = "")
cat("  Figures            : ", normalizePath(opt$outdir), "\n", sep = "")
