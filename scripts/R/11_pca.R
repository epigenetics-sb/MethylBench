!/usr/bin/env Rscript
# =============================================================================
# MethylBench – PCA Analysis
# =============================================================================
# Description:
#   Principal Component Analysis of CpG methylation profiles across methods
#   and samples. PCA is computed on the transposed matrix (samples/methods
#   as observations, CpGs as variables) and PCA SCORES (pca$x) are plotted
#   to visualize sample-level separation by method.
#
#   Six PCA plots are generated:
#     1. Blood  – with EPIC    (EPIC sites only, all methods)
#     2. Fibro  – with EPIC    (EPIC sites only, all methods)
#     3. Blood  – without EPIC (EPIC sites only, sequencing methods)
#     4. Fibro  – without EPIC (EPIC sites only, sequencing methods)
#     5. Blood  – without EPIC (all overlapping CpGs at 10x, sequencing methods)
#     6. Fibro  – without EPIC (all overlapping CpGs at 10x, sequencing methods)
#
# CHANGELOG (post-review fix, Reviewer 1 / Major comment 1):
#   1) FIXED label/column desync:
#      Method/Sample labels used to be built by hand as a SEPARATE vector
#      (e.g. rep(SEQ_METHODS, each = 5): method-major order), while the
#      matrix columns were produced by meth_cols() in SAMPLE-MAJOR order
#      (methods cycle fastest: ONT_Blood1, TWIST_Blood1, WGBS_Blood1,
#      RRBS_Blood1, ONT_Blood2, ...). The two orderings did not match, so
#      many points in Figure 6 carried the wrong Method color / Sample
#      shape. Labels are now parsed directly out of colnames(mat) via
#      parse_meth_labels(), so they cannot desync from the matrix again —
#      there is only one source of truth (the column names actually used
#      to build `mat`), not two independently-maintained orderings.
#   2) FIXED loadings-vs-scores confusion:
#      PCA used to be run on the untransposed matrix (CpGs as rows, i.e.
#      as "observations"; sample-method columns as "variables"), and
#      pca$rotation (variable LOADINGS) was plotted as if it showed
#      sample/method separation. PCA is now run on the TRANSPOSED matrix
#      (samples/methods as rows/observations, CpGs as columns/variables,
#      so each CpG is centered across samples — the standard convention
#      for sample-level PCA), and pca$x (sample SCORES) is plotted, which
#      is the conventional way to visualize separation between samples.
#
# Input:
#   --all_path    Path to ALL.csv (merged matrix with EPIC columns)
#                 Expected column naming after buildMergedMatrix(include_epic=TRUE):
#                   EPIC_Blood1..5, EPIC_Fibro1..5, EPIC_GIAB1..2
#                   ONT_Blood1..5, WGBS_Blood1..5, TWIST_Blood1..5, RRBS_Blood1..5
#                   ONT_Fibro1..5, WGBS_Fibro1..5, TWIST_Fibro1..5, RRBS_Fibro1..5
#                   ONT_cov_*, WGBS_cov_*, TWIST_cov_*, RRBS_cov_*
#   --blood_path  Path to Blood_without_EPIC.csv
#   --fibro_path  Path to Fibro_without_EPIC.csv
#   --outdir      Output directory for figures
#
# Usage:
#   Rscript 05_pca.R \
#     --all_path   data/matrices/ALL_with_EPIC.csv \
#     --blood_path data/matrices/Blood_without_EPIC.csv \
#     --fibro_path data/matrices/Fibro_without_EPIC.csv \
#     --outdir     results/figures/
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path",
    type    = "character",
    help    = "Path to ALL_with_EPIC.csv (EPIC + sequencing matrix) [required]",
    metavar = "FILE"
  ),
  make_option("--blood_path",
    type    = "character",
    help    = "Path to Blood_without_EPIC.csv [required]",
    metavar = "FILE"
  ),
  make_option("--fibro_path",
    type    = "character",
    help    = "Path to Fibro_without_EPIC.csv [required]",
    metavar = "FILE"
  ),
  make_option("--outdir",
    type    = "character",
    default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$all_path))   stop("ERROR: --all_path is required")
if (is.null(opt$blood_path)) stop("ERROR: --blood_path is required")
if (is.null(opt$fibro_path)) stop("ERROR: --fibro_path is required")
for (p in c(opt$all_path, opt$blood_path, opt$fibro_path)) {
  if (!file.exists(p)) stop(paste("File not found:", p))
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

col.vec    <- get_colors()
BLOOD_IDS  <- paste0("Blood", 1:5)
FIBRO_IDS  <- paste0("Fibro",  1:5)
SEQ_METHODS <- c("ONT", "TWIST", "WGEC", "RRBS")  # sequencing-only, no EPIC

# Column name helpers
meth_cols <- function(methods, samples, prefix_map = c("WGEC" = "WGBS")) {
  # Returns methylation column names in method x sample order
  # prefix_map: remap method label to actual column prefix (WGEC -> WGBS)
  sapply(samples, function(smp) {
    sapply(methods, function(m) {
      pfx <- if (m %in% names(prefix_map)) prefix_map[[m]] else m
      paste0(pfx, "_", smp)
    })
  }) |> as.vector()
}

cov_cols <- function(methods, samples, prefix_map = c("WGEC" = "WGBS")) {
  sapply(samples, function(smp) {
    sapply(methods, function(m) {
      pfx <- if (m %in% names(prefix_map)) prefix_map[[m]] else m
      paste0(pfx, "_cov_", smp)
    })
  }) |> as.vector()
}

#' Parse Method and Sample labels directly from a methylation matrix's
#' own column names (e.g. "ONT_Blood1", "WGBS_Fibro3", "EPIC_Blood2").
#'
#' This replaces the previous approach of maintaining a second,
#' independently-ordered label vector (rep(...,each=...) / rep(...,times=...))
#' alongside the column-building logic in meth_cols()/cov_cols(). Deriving
#' the labels from the actual colnames() of the matrix that goes into
#' prcomp() makes it structurally impossible for labels and columns to
#' fall out of sync, regardless of what order meth_cols() happens to
#' produce columns in.
parse_meth_labels <- function(col_names, prefix_map = c("WGBS" = "WGEC")) {
  m  <- regmatches(
    col_names,
    regexec("^([A-Za-z0-9]+)_((?:Blood|Fibro)[0-9]+)$", col_names)
  )
  ok <- lengths(m) == 3
  if (!all(ok)) {
    stop(
      "parse_meth_labels(): could not parse method/sample from column(s): ",
      paste(col_names[!ok], collapse = ", ")
    )
  }
  methods <- vapply(m, `[[`, character(1), 2)
  samples <- vapply(m, `[[`, character(1), 3)
  methods <- unname(ifelse(methods %in% names(prefix_map),
                            prefix_map[methods], methods))
  list(methods = methods, samples = samples)
}

runPCA <- function(mat, label = "") {
  # mat: CpGs (rows) x sample/method columns (colnames like "ONT_Blood1")
  labels <- parse_meth_labels(colnames(mat))

  # Transpose so that samples/methods are the OBSERVATIONS (rows) and CpGs
  # are the VARIABLES (columns). prcomp() then centers each CpG across
  # samples (the standard convention for sample-level PCA), and pca$x
  # (scores) has exactly one row per input column, in the original column
  # order -- so `labels` is guaranteed to line up with it.
  # scale. = FALSE: beta-values already live on a common [0,1] scale;
  # scaling to unit variance would inflate the influence of low-variance
  # (e.g. near-0 or near-1) CpGs relative to biologically variable ones.
  pca      <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pca_mat  <- as.data.frame(pca$x)
  var_expl <- summary(pca)$importance[2, ] * 100  # proportion of variance

  cat(sprintf("\n--- PCA: %s ---\n", label))
  print(summary(pca))

  pca_mat$Method <- labels$methods
  pca_mat$Sample <- labels$samples

  return(list(pca_mat = pca_mat, var_expl = var_expl))
}

plot_pca <- function(pca_mat, var_expl, filename, height = 12, width = 15) {

  pc1_label <- sprintf("PC1 (%.2f%%)", var_expl["PC1"])
  pc2_label <- sprintf("PC2 (%.2f%%)", var_expl["PC2"])

  p <- ggplot(pca_mat, aes(x = PC1, y = PC2, color = Method, shape = Sample)) +
    geom_point(size = 9) +
    scale_color_manual(values = col.vec) +
    labs(x = pc1_label, y = pc2_label) +
    theme_bw() +
    theme(
      axis.text  = element_text(size = 26),
      axis.title = element_text(size = 26),
      text       = element_text(size = 26)
    )

  ggsave(p,
    filename = file.path(opt$outdir, filename),
    height = height, width = width, dpi = 300
  )
  cat(sprintf("  Saved: %s\n", filename))
}

cat("[1/3] Loading matrices...\n")

all   <- fread(opt$all_path,   header = TRUE, sep = ",", na.strings = "NA")
blood <- fread(opt$blood_path, header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(opt$fibro_path, header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  ALL   : %d CpGs x %d columns\n", nrow(all),   ncol(all)))
cat(sprintf("  Blood : %d CpGs x %d columns\n", nrow(blood), ncol(blood)))
cat(sprintf("  Fibro : %d CpGs x %d columns\n", nrow(fibro), ncol(fibro)))

cat("[2/3] Running PCAs on ALL matrix (with EPIC)...\n")

# Blood with EPIC -------------------------------------------------------
bl_epic_cols <- c(
  paste0("EPIC_",  BLOOD_IDS),
  meth_cols(SEQ_METHODS, BLOOD_IDS)
)
bl_epic_cols <- intersect(bl_epic_cols, colnames(all))

bl_epic      <- all[, ..bl_epic_cols]
bl_epic      <- bl_epic[complete.cases(bl_epic)]

res <- runPCA(bl_epic, label = "Blood with EPIC")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_EPIC.png")

# Fibro with EPIC -------------------------------------------------------
fb_epic_cols <- c(
  paste0("EPIC_",  FIBRO_IDS),
  meth_cols(SEQ_METHODS, FIBRO_IDS)
)
fb_epic_cols <- intersect(fb_epic_cols, colnames(all))

fb_epic      <- all[, ..fb_epic_cols]
fb_epic      <- fb_epic[complete.cases(fb_epic)]

res <- runPCA(fb_epic, label = "Fibro with EPIC")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_EPIC.png")

cat("[2/3] Running PCAs on ALL matrix (sequencing methods on EPIC sites)...\n")

# Blood without EPIC, on EPIC sites ------------------------------------
bl_seq_cols  <- meth_cols(SEQ_METHODS, BLOOD_IDS)
bl_seq_cols  <- intersect(bl_seq_cols, colnames(all))

# Use EPIC completeness to restrict to EPIC-covered CpGs
epic_mask_bl <- complete.cases(all[, ..bl_epic_cols])
bl_seq       <- all[epic_mask_bl, ..bl_seq_cols]
bl_seq       <- bl_seq[complete.cases(bl_seq)]

res <- runPCA(bl_seq, label = "Blood without EPIC (on EPIC sites)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_noEPIC_onEPICsites.png")

# Fibro without EPIC, on EPIC sites ------------------------------------
fb_seq_cols  <- meth_cols(SEQ_METHODS, FIBRO_IDS)
fb_seq_cols  <- intersect(fb_seq_cols, colnames(all))

epic_mask_fb <- complete.cases(all[, ..fb_epic_cols])
fb_seq       <- all[epic_mask_fb, ..fb_seq_cols]
fb_seq       <- fb_seq[complete.cases(fb_seq)]

res <- runPCA(fb_seq, label = "Fibro without EPIC (on EPIC sites)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_noEPIC_onEPICsites.png")

cat("[3/3] Running PCAs on without-EPIC matrices (10x filter, all CpGs)...\n")

# Blood without EPIC, 10x ----------------------------------------------
bl_cov_cols <- cov_cols(SEQ_METHODS, BLOOD_IDS)
bl_cov_cols <- intersect(bl_cov_cols, colnames(blood))

bl_10x      <- extractCovDf(blood, threshold = 10, cov_cols = bl_cov_cols)
bl_meth_cols <- meth_cols(SEQ_METHODS, BLOOD_IDS)
bl_meth_cols <- intersect(bl_meth_cols, colnames(bl_10x))
bl_10x       <- bl_10x[, ..bl_meth_cols]
bl_10x       <- bl_10x[complete.cases(bl_10x)]

res <- runPCA(bl_10x, label = "Blood without EPIC (10x, all CpGs)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_noEPIC.png")

# Fibro without EPIC, 10x ---------------------------------------------
fb_cov_cols  <- cov_cols(SEQ_METHODS, FIBRO_IDS)
fb_cov_cols  <- intersect(fb_cov_cols, colnames(fibro))

fb_10x       <- extractCovDf(fibro, threshold = 10, cov_cols = fb_cov_cols)
fb_meth_cols <- meth_cols(SEQ_METHODS, FIBRO_IDS)
fb_meth_cols <- intersect(fb_meth_cols, colnames(fb_10x))
fb_10x       <- fb_10x[, ..fb_meth_cols]
fb_10x       <- fb_10x[complete.cases(fb_10x)]

res <- runPCA(fb_10x, label = "Fibro without EPIC (10x, all CpGs)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_noEPIC.png")

cat(sprintf("\nDone. 6 figures written to: %s\n", opt$outdir))
