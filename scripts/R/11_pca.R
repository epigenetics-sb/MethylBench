#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - principal component analysis of sample-platform profiles
# =============================================================================
# PCA per tissue with sample-platform profiles as observations and CpGs as
# variables (centered, not scaled); sample scores are plotted.
# Figure 6 (A/B: with EPIC; C/D: sequencing, >= 10x; E/F: sequencing on the
# EPIC CpG set) and Suppl. Figure 1 (GIAB, unfiltered and >= 10x).
#
# Usage:
#   Rscript scripts/R/11_pca.R --all_path ALL.csv \
#     --blood_path Blood_without_EPIC.csv --fibro_path Fibro_without_EPIC.csv \
#     --giab_path GIAB_without_EPIC.csv --outdir results/figures/
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
  make_option("--giab_path",
    type    = "character",
    help    = "Path to GIAB_without_EPIC.csv [required]",
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
if (is.null(opt$giab_path))  stop("ERROR: --giab_path is required")
for (p in c(opt$all_path, opt$blood_path, opt$fibro_path, opt$giab_path)) {
  if (!file.exists(p)) stop(paste("File not found:", p))
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

col.vec    <- get_colors()
BLOOD_IDS  <- paste0("Blood", 1:5)
FIBRO_IDS  <- paste0("Fibro",  1:5)
GIAB_IDS   <- paste0("GIAB",   1:2)
SEQ_METHODS      <- c("ONT", "TWIST", "WGEC", "RRBS")
SEQ_METHODS_GIAB <- c("ONT", "PacBio", "TWIST", "WGEC", "RRBS")

meth_cols <- function(methods, samples) {
  sapply(samples, function(smp) {
    sapply(methods, function(m) paste0(m, "_", smp))
  }) |> as.vector()
}

cov_cols <- function(methods, samples) {
  sapply(samples, function(smp) {
    sapply(methods, function(m) paste0(m, "_cov_", smp))
  }) |> as.vector()
}

# Method and sample labels are parsed from the column names of the analysed
# matrix, so labels cannot get out of order with the columns.
parse_meth_labels <- function(col_names) {
  m  <- regmatches(
    col_names,
    regexec("^([A-Za-z0-9]+)_((?:Blood|Fibro|GIAB)[0-9]+)$", col_names)
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
  list(methods = methods, samples = samples)
}

runPCA <- function(mat, label = "") {
  labels <- parse_meth_labels(colnames(mat))

  # Transposed: profiles = observations, CpGs = variables.
  pca      <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  pca_mat  <- as.data.frame(pca$x)
  var_expl <- summary(pca)$importance[2, ] * 100

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

# ---- 1. Load matrices -------------------------------------------------------

cat("[1/4] Loading matrices...\n")

all   <- fread(opt$all_path,   header = TRUE, sep = ",", na.strings = "NA")
blood <- fread(opt$blood_path, header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(opt$fibro_path, header = TRUE, sep = ",", na.strings = "NA")
giab  <- fread(opt$giab_path,  header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  ALL   : %d CpGs x %d columns\n", nrow(all),   ncol(all)))
cat(sprintf("  Blood : %d CpGs x %d columns\n", nrow(blood), ncol(blood)))
cat(sprintf("  Fibro : %d CpGs x %d columns\n", nrow(fibro), ncol(fibro)))
cat(sprintf("  GIAB  : %d CpGs x %d columns\n", nrow(giab),  ncol(giab)))

# ---- 2. Figure 6A/B: EPIC + sequencing, CpGs complete in all profiles -------

cat("[2/4] Running PCAs on ALL matrix (with EPIC)...\n")

bl_epic_cols <- c(
  paste0("EPIC_",  BLOOD_IDS),
  meth_cols(SEQ_METHODS, BLOOD_IDS)
)
bl_epic_cols <- intersect(bl_epic_cols, colnames(all))

bl_epic      <- all[, ..bl_epic_cols]
bl_epic      <- bl_epic[complete.cases(bl_epic)]

res <- runPCA(bl_epic, label = "Blood with EPIC")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_EPIC.png")

fb_epic_cols <- c(
  paste0("EPIC_",  FIBRO_IDS),
  meth_cols(SEQ_METHODS, FIBRO_IDS)
)
fb_epic_cols <- intersect(fb_epic_cols, colnames(all))

fb_epic      <- all[, ..fb_epic_cols]
fb_epic      <- fb_epic[complete.cases(fb_epic)]

res <- runPCA(fb_epic, label = "Fibro with EPIC")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_EPIC.png")

# ---- Figure 6E/F: sequencing only, same CpG set as 6A/B ---------------------

cat("[2/4] Running PCAs on ALL matrix (sequencing methods on EPIC sites)...\n")

bl_seq_cols  <- meth_cols(SEQ_METHODS, BLOOD_IDS)
bl_seq_cols  <- intersect(bl_seq_cols, colnames(all))

epic_mask_bl <- complete.cases(all[, ..bl_epic_cols])
bl_seq       <- all[epic_mask_bl, ..bl_seq_cols]
bl_seq       <- bl_seq[complete.cases(bl_seq)]

res <- runPCA(bl_seq, label = "Blood without EPIC (on EPIC sites)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_noEPIC_onEPICsites.png")

fb_seq_cols  <- meth_cols(SEQ_METHODS, FIBRO_IDS)
fb_seq_cols  <- intersect(fb_seq_cols, colnames(all))

epic_mask_fb <- complete.cases(all[, ..fb_epic_cols])
fb_seq       <- all[epic_mask_fb, ..fb_seq_cols]
fb_seq       <- fb_seq[complete.cases(fb_seq)]

res <- runPCA(fb_seq, label = "Fibro without EPIC (on EPIC sites)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_noEPIC_onEPICsites.png")

# ---- 3. Figure 6C/D: sequencing only, >= 10x in every profile ---------------

cat("[3/4] Running PCAs on without-EPIC matrices (10x filter, all CpGs)...\n")

bl_cov_cols <- cov_cols(SEQ_METHODS, BLOOD_IDS)
bl_cov_cols <- intersect(bl_cov_cols, colnames(blood))

bl_10x      <- extractCovDf(blood, threshold = 10, cov_cols = bl_cov_cols)
bl_meth_cols <- meth_cols(SEQ_METHODS, BLOOD_IDS)
bl_meth_cols <- intersect(bl_meth_cols, colnames(bl_10x))
bl_10x       <- bl_10x[, ..bl_meth_cols]
bl_10x       <- bl_10x[complete.cases(bl_10x)]

res <- runPCA(bl_10x, label = "Blood without EPIC (10x, all CpGs)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Blood_noEPIC.png")

fb_cov_cols  <- cov_cols(SEQ_METHODS, FIBRO_IDS)
fb_cov_cols  <- intersect(fb_cov_cols, colnames(fibro))

fb_10x       <- extractCovDf(fibro, threshold = 10, cov_cols = fb_cov_cols)
fb_meth_cols <- meth_cols(SEQ_METHODS, FIBRO_IDS)
fb_meth_cols <- intersect(fb_meth_cols, colnames(fb_10x))
fb_10x       <- fb_10x[, ..fb_meth_cols]
fb_10x       <- fb_10x[complete.cases(fb_10x)]

res <- runPCA(fb_10x, label = "Fibro without EPIC (10x, all CpGs)")
plot_pca(res$pca_mat, res$var_expl, "PCA_Fibro_noEPIC.png")

# ---- 4. Suppl. Figure 1: GIAB, unfiltered and >= 10x ------------------------

cat("[4/4] Running GIAB PCA (Supplementary Figure 1, five sequencing methods incl. PacBio)...\n")

giab_meth_cols <- meth_cols(SEQ_METHODS_GIAB, GIAB_IDS)
giab_meth_cols <- intersect(giab_meth_cols, colnames(giab))

giab_unfiltered <- giab[, ..giab_meth_cols]
giab_unfiltered <- giab_unfiltered[complete.cases(giab_unfiltered)]

res <- runPCA(giab_unfiltered, label = "GIAB (unfiltered)")
plot_pca(res$pca_mat, res$var_expl, "PCA_GIAB_unfiltered.png")

giab_cov_cols <- cov_cols(SEQ_METHODS_GIAB, GIAB_IDS)
giab_cov_cols <- intersect(giab_cov_cols, colnames(giab))

giab_10x        <- extractCovDf(giab, threshold = 10, cov_cols = giab_cov_cols)
giab_10x_meth   <- intersect(giab_meth_cols, colnames(giab_10x))
giab_10x        <- giab_10x[, ..giab_10x_meth]
giab_10x        <- giab_10x[complete.cases(giab_10x)]

res <- runPCA(giab_10x, label = "GIAB (10x coverage filter)")
plot_pca(res$pca_mat, res$var_expl, "PCA_GIAB_10x.png")

cat(sprintf("\nDone. 8 figures written to: %s\n", opt$outdir))
