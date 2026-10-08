#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - cross-platform correlation analysis
# =============================================================================
# Pairwise Pearson correlations of CpG methylation between platforms at
# increasing coverage thresholds, overall and per methylation stratum.
# Figure 4; Suppl. Figures 2 and 19; Suppl. Table S7.
#
# Usage:
#   Rscript scripts/R/09_correlation_analysis.R --datadir <matrices>/ --outdir results/figures/
#   (--datadir must contain Blood_, Fibro_ and GIAB_without_EPIC.csv)
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

METHODS_NO_PACBIO <- c("ONT" = "ONT", "WGEC" = "WGEC",
                        "RRBS" = "RRBS", "TWIST" = "TWIST")
METHODS_PACBIO    <- c("ONT" = "ONT", "WGEC" = "WGEC", "RRBS" = "RRBS",
                        "TWIST" = "TWIST", "PacBio" = "PacBio")
METHODS_HIGHCOV   <- c("ONT" = "ONT", "PacBio" = "PacBio", "TWIST" = "TWIST")
HIGHCOV_COMPARISON_COLORS <- c(
  "ONT vs PacBio"   = "dodgerblue2",
  "ONT vs TWIST"    = "green4",
  "PacBio vs TWIST" = "black"
)
COVERAGES         <- c(0, 5, 10, 15)
COVERAGES_HIGH    <- c(0, 5, 10, 15, 20, 25, 30, 35, 40)
COV_MAX_PLOT      <- 15

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

cat("[3/5] Computing ONT vs. PacBio vs. TWIST high-coverage correlations...\n")

corr_ont_pacbio_twist <- computeCorrAcrossCoverages(
  data      = giab,
  samples   = c("GIAB1", "GIAB2"),
  methods   = METHODS_HIGHCOV,
  coverages = COVERAGES_HIGH
)

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

# ---- Figure 4A / Suppl. Figure 2: correlation vs. coverage threshold --------

p_blood <- plot_corr(
  corr_blood, ncols = 5,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_blood,
  filename = file.path(opt$outdir, "Correlations_Blood.png"),
  height = 12, width = 14, dpi = 300
)

p_fibro <- plot_corr(
  corr_fibro, ncols = 5,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_fibro,
  filename = file.path(opt$outdir, "Correlations_Fibro.png"),
  height = 12, width = 14, dpi = 300
)

p_giab <- plot_corr(
  corr_giab, ncols = 2,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_giab,
  filename = file.path(opt$outdir, "Correlations_GIAB.png"),
  height = 12, width = 14, dpi = 300
)

# ---- Representative high-coverage samples (Blood3, Fibro4, GIAB2) -----------

p_highcov <- plot_corr(
  corr_highcov, ncols = 3,
  title = "Correlation Changes with respect\nto increasing Coverage thresholds"
)
ggsave(p_highcov,
  filename = file.path(opt$outdir, "Correlations_High_Cov.png"),
  height = 12, width = 14, dpi = 300
)

# ---- Figure 4B: ONT, PacBio and TWIST up to 40x (GIAB) ----------------------

p_ont_twist <- ggplot(
  corr_ont_pacbio_twist,
  aes(x = Coverage2, y = Correlation, color = Comparison)
) +
  geom_point(size = 4) +
  geom_line(aes(group = Comparison)) +
  scale_color_manual(values = HIGHCOV_COMPARISON_COLORS) +
  facet_wrap(~sample, ncol = 2) +
  scale_x_discrete(
    labels = ifelse(COVERAGES_HIGH == 0, "None", paste0(COVERAGES_HIGH, "x"))
  ) +
  labs(
    title = "Correlation Changes with respect\nto increasing Coverage thresholds in High Coverage samples",
    x     = "Coverage Filter",
    y     = "Pearson Correlation coefficient"
  ) +
  theme_corr()

ggsave(p_ont_twist,
  filename = file.path(opt$outdir, "Correlations_GIAB_High_Cov.png"),
  height = 12, width = 14, dpi = 300
)

# ---- Suppl. Figure 19, Suppl. Table S7: correlation per methylation stratum ----

cat("[5/6] Computing correlations per methylation stratum...\n")

STRATUM_LOWER  <- 0.2
STRATUM_UPPER  <- 0.8
STRATUM_LEVELS <- c("All", "Low", "Intermediate", "High")
STRATUM_LABELS <- c(
  All          = "All CpGs",
  Low          = sprintf("Low (mean beta < %.1f)", STRATUM_LOWER),
  Intermediate = sprintf("Intermediate (%.1f-%.1f)", STRATUM_LOWER, STRATUM_UPPER),
  High         = sprintf("High (mean beta > %.1f)", STRATUM_UPPER)
)

# Strata are assigned from the mean beta of the two compared platforms, so that
# no single platform's measurement error determines stratum membership. The
# "All" stratum reproduces the correlations of Figure 4.
computeStratifiedCorr <- function(data, samples, methods, coverages,
                                  lower = STRATUM_LOWER, upper = STRATUM_UPPER) {
  stopifnot(is.data.table(data))
  out <- list()
  for (smp in samples) {
    for (cov in coverages) {
      if (cov == 0) {
        filtered <- data
      } else {
        cov_cols <- paste0(methods, "_cov_", smp)
        cov_cols <- cov_cols[cov_cols %in% colnames(data)]
        filtered <- extractCovDf(data, threshold = cov, cov_cols = cov_cols)
      }
      for (pair in combn(names(methods), 2, simplify = FALSE)) {
        c1 <- paste0(methods[pair[1]], "_", smp)
        c2 <- paste0(methods[pair[2]], "_", smp)
        if (!all(c(c1, c2) %in% colnames(filtered))) next
        x  <- filtered[[c1]]
        y  <- filtered[[c2]]
        ok <- !is.na(x) & !is.na(y)
        x  <- x[ok]; y <- y[ok]
        if (length(x) == 0) next
        scale_f <- if (max(c(x, y)) > 1.5) 100 else 1
        m <- (x + y) / (2 * scale_f)
        stratum <- ifelse(m < lower, "Low", ifelse(m > upper, "High", "Intermediate"))
        for (s in STRATUM_LEVELS) {
          idx <- if (s == "All") rep(TRUE, length(x)) else stratum == s
          n   <- sum(idx)
          out[[length(out) + 1]] <- data.frame(
            Sample     = smp,
            Coverage   = cov,
            Comparison = paste(pair[1], "vs", pair[2]),
            Stratum    = s,
            n_CpGs     = n,
            frac_CpGs  = n / length(x),
            Pearson    = if (n >= 3) cor(x[idx], y[idx], method = "pearson")  else NA_real_,
            Spearman   = if (n >= 3) cor(x[idx], y[idx], method = "spearman") else NA_real_,
            MAD        = if (n >= 1) mean(abs(x[idx] - y[idx])) / scale_f    else NA_real_
          )
        }
      }
    }
  }
  res <- do.call(rbind, out)
  res$Coverage2 <- factor(ifelse(res$Coverage == 0, "None", paste0(res$Coverage, "x")),
                          levels = ifelse(coverages == 0, "None", paste0(coverages, "x")))
  res$Stratum   <- factor(res$Stratum, levels = STRATUM_LEVELS)
  res
}

strat_blood <- computeStratifiedCorr(blood, paste0("Blood", 1:5),     METHODS_NO_PACBIO, COVERAGES)
strat_fibro <- computeStratifiedCorr(fibro, paste0("Fibro", 1:5),     METHODS_NO_PACBIO, COVERAGES)
strat_giab  <- computeStratifiedCorr(giab,  c("GIAB1", "GIAB2"),      METHODS_PACBIO,    COVERAGES)
strat_blood$Tissue <- "Blood"; strat_fibro$Tissue <- "Fibroblast"; strat_giab$Tissue <- "GIAB"
strat_all <- rbind(strat_blood, strat_fibro, strat_giab)

fwrite(strat_all, file.path(opt$outdir, "Correlation_by_methylation_stratum.tsv"), sep = "\t")

# ---- Suppl. Table S7 (10x) --------------------------------------------------
strat_10x <- as.data.table(strat_all)[Coverage == 10]
table_s7 <- strat_10x[, .(
  n_samples        = .N,
  median_n_CpGs    = as.numeric(median(n_CpGs)),
  mean_frac_pct    = round(100 * mean(frac_CpGs), 1),
  Pearson_mean     = round(mean(Pearson,  na.rm = TRUE), 3),
  Pearson_min      = round(min(Pearson,   na.rm = TRUE), 3),
  Pearson_max      = round(max(Pearson,   na.rm = TRUE), 3),
  Spearman_mean    = round(mean(Spearman, na.rm = TRUE), 3),
  MAD_mean         = round(mean(MAD,      na.rm = TRUE), 3)
), by = .(Tissue, Comparison, Stratum)][order(Tissue, Stratum, Comparison)]
fwrite(table_s7, file.path(opt$outdir, "Table_S7_intermediate_correlation.tsv"), sep = "\t")

# Consistency check: "All" stratum vs. Figure 4 correlations.
chk <- merge(
  as.data.table(strat_blood)[Stratum == "All", .(Sample, Coverage, Comparison, Pearson)],
  as.data.table(corr_blood)[, .(Sample = sub("^_", "", Sample), Coverage, Comparison, Correlation)],
  by = c("Sample", "Coverage", "Comparison")
)
if (nrow(chk) > 0) {
  cat(sprintf("  Check vs. Figure 4 (Blood): max |r_All - r_Fig4| = %.2e (n = %d)\n",
              max(abs(chk$Pearson - chk$Correlation), na.rm = TRUE), nrow(chk)))
}

cat("  Intermediate stratum at 10x (mean Pearson / Spearman / MAD per tissue):\n")
print(table_s7[Stratum == "Intermediate",
               .(Pearson = round(mean(Pearson_mean), 3),
                 Pearson_range = sprintf("%.3f-%.3f", min(Pearson_min), max(Pearson_max)),
                 Spearman = round(mean(Spearman_mean), 3),
                 MAD = round(mean(MAD_mean), 3),
                 frac_pct = round(mean(mean_frac_pct), 1)),
               by = Tissue])

# ---- Suppl. Figure 19 (single panels) ---------------------------------------
cat("[6/6] Plotting stratified correlations...\n")

save_panel <- function(p, name, width, height) {
  ggsave(file.path(opt$outdir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
  ggsave(file.path(opt$outdir, paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf)
}

plot_strat <- function(df) {
  df <- as.data.table(df)[, .(Pearson = mean(Pearson, na.rm = TRUE),
                              lo = min(Pearson, na.rm = TRUE),
                              hi = max(Pearson, na.rm = TRUE)),
                          by = .(Coverage2, Comparison, Stratum)]
  ggplot(df, aes(x = Coverage2, y = Pearson, color = Comparison, group = Comparison)) +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = 0.4, linewidth = 1.2) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = C25) +
    facet_wrap(~Stratum, nrow = 1, scales = "free_y", labeller = as_labeller(STRATUM_LABELS)) +
    labs(x = "Coverage Filter", y = "Pearson correlation coefficient") +
    guides(color = guide_legend(nrow = 2)) +
    theme_corr() +
    theme(strip.text = element_text(size = 16))
}

save_panel(plot_strat(strat_blood), "SupplFig19A_Blood", width = 16, height = 7)
save_panel(plot_strat(strat_fibro), "SupplFig19B_Fibro", width = 16, height = 7)
save_panel(plot_strat(strat_giab),  "SupplFig19C_GIAB",  width = 16, height = 7)

cat(sprintf("\n[6/6] Done. Figures and tables written to: %s\n", opt$outdir))
