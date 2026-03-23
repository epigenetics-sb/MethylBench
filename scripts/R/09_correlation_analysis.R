#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Correlation Analysis
# =============================================================================
# Description:
#   Computes pairwise Pearson correlations between all methylation platforms
#   across increasing coverage thresholds. Generates per-sampleset correlation
#   line plots (Blood, Fibroblast, GIAB) and an extended high-coverage
#   ONT vs TWIST comparison for GIAB.
#
# Input:
#   --datadir   Directory containing merged methylation matrices
#               (Blood_without_EPIC.csv, Fibro_without_EPIC.csv,
#                GIAB_without_EPIC.csv) as produced by buildMergedMatrix()
#   --outdir    Output directory for figures
#
# Usage:
#   Rscript 03_correlation_analysis.R \
#     --datadir data/matrices/ \
#     --outdir  results/figures/
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
  make_option("--datadir",
    type    = "character",
    help    = "Directory containing merged methylation matrices [required]",
    metavar = "DIR"
  ),
  make_option("--outdir",
    type    = "character",
    default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$datadir)) stop("ERROR: --datadir is required")
if (!dir.exists(opt$datadir)) stop(paste("Directory not found:", opt$datadir))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

METHODS_NO_PACBIO <- c("ONT" = "ONT", "WGEC" = "WGBS",
                        "RRBS" = "RRBS", "TWIST" = "TWIST")
METHODS_PACBIO    <- c("ONT" = "ONT", "WGEC" = "WGBS", "RRBS" = "RRBS",
                        "TWIST" = "TWIST", "PacBio" = "PacBio")
COVERAGES         <- c(0, 5, 10, 15)
COVERAGES_HIGH    <- c(0, 5, 10, 15, 20, 25, 30, 35, 40)
COV_MAX_PLOT      <- 15

# Shared plot theme
theme_corr <- function() {
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text    = element_text(size = 22),
    axis.title.x = element_text(size = 22, vjust = -0.5),
    text         = element_text(size = 22),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22),
    legend.position = "bottom",
    legend.title    = element_blank()
  )
}

cat("[1/5] Loading merged methylation matrices...\n")

blood <- fread(file.path(opt$datadir, "Blood_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(file.path(opt$datadir, "Fibro_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")
giab  <- fread(file.path(opt$datadir, "GIAB_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  Blood : %d CpGs\n", nrow(blood)))
cat(sprintf("  Fibro : %d CpGs\n", nrow(fibro)))
cat(sprintf("  GIAB  : %d CpGs\n", nrow(giab)))

cat("[2/5] Computing correlations...\n")

corr_blood <- computeCorrAcrossCoverages(
  data      = blood,
  samples   = paste0("Blood", 1:5),
  methods   = METHODS_NO_PACBIO,
  coverages = COVERAGES
)

corr_fibro <- computeCorrAcrossCoverages(
  data      = fibro,
  samples   = paste0("Fibro", 1:5),
  methods   = METHODS_NO_PACBIO,
  coverages = COVERAGES
)

corr_giab <- computeCorrAcrossCoverages(
  data      = giab,
  samples   = c("GIAB1", "GIAB2"),
  methods   = METHODS_PACBIO,
  coverages = COVERAGES
)

corr_highcov_blood <- computeCorrAcrossCoverages(
  data      = blood,
  samples   = "Blood3",
  methods   = METHODS_NO_PACBIO,
  coverages = COVERAGES
)
corr_highcov_fibro <- computeCorrAcrossCoverages(
  data      = fibro,
  samples   = "Fibro4",
  methods   = METHODS_NO_PACBIO,
  coverages = COVERAGES
)
corr_highcov_giab <- computeCorrAcrossCoverages(
  data      = giab,
  samples   = "GIAB2",
  methods   = METHODS_PACBIO,
  coverages = COVERAGES
)
corr_highcov <- rbind(corr_highcov_blood, corr_highcov_fibro, corr_highcov_giab)

cat("[3/5] Computing ONT vs TWIST high-coverage correlation...\n")

ont_twist_rows <- lapply(COVERAGES_HIGH, function(cov_thr) {
  lapply(c("GIAB1", "GIAB2"), function(smp) {
    if (cov_thr == 0) {
      filtered <- giab
    } else {
      cov_cols <- c(
        paste0("ONT_cov_",   smp),
        paste0("TWIST_cov_", smp)
      )
      cov_cols <- cov_cols[cov_cols %in% colnames(giab)]
      filtered <- extractCovDf(giab, threshold = cov_thr, cov_cols = cov_cols)
    }
    corr_val <- cor(
      filtered[[paste0("ONT_",   smp)]],
      filtered[[paste0("TWIST_", smp)]],
      use = "pairwise.complete.obs"
    )
    data.frame(
      Sample      = smp,
      Correlation = corr_val,
      Cov         = cov_thr,
      stringsAsFactors = FALSE
    )
  })
})
corr_ont_twist <- rbindlist(unlist(ont_twist_rows, recursive = FALSE))
corr_ont_twist[, Coverage2 := ifelse(Cov == 0, "None", paste0(Cov, "x"))]
corr_ont_twist[, Coverage2 := factor(
  Coverage2,
  levels = ifelse(COVERAGES_HIGH == 0, "None", paste0(COVERAGES_HIGH, "x"))
)]

cat("[4/5] Generating figures...\n")

plot_corr <- function(df, ncols, title) {
  ggplot(
    df[df$Coverage <= COV_MAX_PLOT, ],
    aes(x = Coverage2, y = Correlation, color = Comparison)
  ) +
    geom_point(size = 4) +
    geom_line(aes(group = Comparison)) +
    scale_color_manual(values = C25) +
    facet_wrap(~sample, ncol = ncols) +
    labs(
      title = title,
      x     = "Coverage Filter",
      y     = "Pearson Correlation coefficient"
    ) +
    guides(color = guide_legend(nrow = 2)) +
    scale_x_discrete(
      labels = ifelse(COVERAGES == 0, "None", paste0(COVERAGES, "x"))
    ) +
    theme_corr()
}

# ---- 5.1 Blood correlations -------------------------------------------------

p_blood <- plot_corr(
  corr_blood, ncols = 5,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_blood,
  filename = file.path(opt$outdir, "Correlations_Blood.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.2 Fibroblast correlations --------------------------------------------

p_fibro <- plot_corr(
  corr_fibro, ncols = 5,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_fibro,
  filename = file.path(opt$outdir, "Correlations_Fibro.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.3 GIAB correlations --------------------------------------------------

p_giab <- plot_corr(
  corr_giab, ncols = 2,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_giab,
  filename = file.path(opt$outdir, "Correlations_GIAB.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.4 High-coverage representative samples (Blood3, Fibro4, GIAB2) -------

p_highcov <- plot_corr(
  corr_highcov, ncols = 3,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_highcov,
  filename = file.path(opt$outdir, "Correlations_High_Cov.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.5 ONT vs TWIST extended (0–40x, GIAB1 + GIAB2) ----------------------

p_ont_twist <- ggplot(
  corr_ont_twist,
  aes(x = Coverage2, y = Correlation, color = Sample)
) +
  geom_point(size = 4) +
  geom_line(aes(group = Sample)) +
  scale_x_discrete(
    labels = ifelse(COVERAGES_HIGH == 0, "None", paste0(COVERAGES_HIGH, "x"))
  ) +
  labs(
    title = "Correlation Changes with respect\nto increasing Coverage thresholds in ONT vs. TWIST",
    x     = "Coverage Filter",
    y     = "Pearson Correlation coefficient"
  ) +
  theme_corr()

ggsave(p_ont_twist,
  filename = file.path(opt$outdir, "Correlations_GIAB_High_Cov.png"),
  height = 12, width = 14, dpi = 300
)

cat(sprintf("\n[5/5] Done. 5 figures written to: %s\n", opt$outdir))
