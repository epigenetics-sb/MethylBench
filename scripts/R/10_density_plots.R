#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Methylation Density Plots
# =============================================================================
# Description:
#   Generates ridge/density plots of CpG methylation distributions across
#   methods, samples and coverage thresholds. Covers:
#     - Per-sampleset full density (Blood, Fibro, GIAB) – 0 and 10x filter
#     - High-coverage representative samples (Blood3, Fibro4, GIAB2) –
#       0/10/20/30/40x filter, split into low (0–30%) and high (70–100%)
#       methylation windows
#     - GIAB2 high-coverage ONT/TWIST/PacBio only (0–40x)
#     - With-EPIC density at 10x (Blood3, Fibro4, GIAB2)
#
# Input:
#   --datadir    Directory containing merged methylation matrices
#                (Blood_without_EPIC.csv, Fibro_without_EPIC.csv,
#                 GIAB_without_EPIC.csv, and optionally ALL.csv for EPIC)
#   --outdir     Output directory for figures
#   --epic_path  Path to ALL.csv / EPIC-merged matrix (optional)
#
# Usage:
#   Rscript 04_density_plots.R \
#     --datadir  data/matrices/ \
#     --outdir   results/figures/ \
#     --epic_path data/matrices/ALL.csv
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(stringr)
  library(forcats)
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
  ),
  make_option("--epic_path",
    type    = "character",
    default = NULL,
    help    = "Path to ALL.csv (EPIC-merged matrix) for with-EPIC density plot [optional]",
    metavar = "FILE"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$datadir))   stop("ERROR: --datadir is required")
if (!dir.exists(opt$datadir)) stop(paste("Directory not found:", opt$datadir))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

METHODS_NO_PACBIO <- c("ONT" = "ONT", "WGEC" = "WGBS",
                        "RRBS" = "RRBS", "TWIST" = "TWIST")
METHODS_PACBIO    <- c("ONT" = "ONT", "WGEC" = "WGBS", "RRBS" = "RRBS",
                        "TWIST" = "TWIST", "PacBio" = "PacBio")
METHODS_HIGH_COV  <- c("ONT" = "ONT", "TWIST" = "TWIST", "PacBio" = "PacBio")

METHOD_ORDER      <- c("ONT", "PacBio", "RRBS", "TWIST", "WGEC")

COV_COLORS <- c(
  "None" = "black",
  "10x"  = "red",
  "20x"  = "blue",
  "30x"  = "darkgreen",
  "40x"  = "magenta"
)

col.vec <- get_colors()

# Shared theme
theme_density <- function(base_size = 22) {
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = base_size),
    axis.text    = element_text(size = base_size),
    axis.text.x  = element_text(size = base_size),
    axis.title.x = element_text(size = base_size, vjust = -0.5),
    text         = element_text(size = base_size),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank()
  )
}

cat("[1/4] Loading merged methylation matrices...\n")

blood <- fread(file.path(opt$datadir, "Blood_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(file.path(opt$datadir, "Fibro_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")
giab  <- fread(file.path(opt$datadir, "GIAB_without_EPIC.csv"),
               header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  Blood: %d CpGs | Fibro: %d CpGs | GIAB: %d CpGs\n",
  nrow(blood), nrow(fibro), nrow(giab)))

prep_long <- function(dt, method_order = METHOD_ORDER) {
  setDT(dt)
  dt[, SampleGroup := interaction(Sample, Method, sep = "_")]
  dt[, Coverage    := as.factor(Coverage)]
  dt[, Coverage2   := ifelse(Coverage == "0", "None", paste0(Coverage, "x"))]
  dt[, Coverage2   := factor(Coverage2, levels = unique(Coverage2))]
  dt[, Method      := factor(Method, levels = method_order)]
  return(dt)
}

cat("[2/4] Building long-format data...\n")

# ---- 5.1 Blood – full density (0 and 10x) -----------------------------------

tp_blood <- buildMethLong(
  data      = blood,
  samples   = paste0("Blood", 1:5),
  methods   = METHODS_NO_PACBIO,
  coverages = c(0, 10)
)
tp_blood <- prep_long(tp_blood)
tp_blood[, SampleGroup := factor(SampleGroup, levels = c(
  paste0("Blood", 1:5, "_WGEC"),  paste0("Blood", 1:5, "_TWIST"),
  paste0("Blood", 1:5, "_RRBS"),  paste0("Blood", 1:5, "_ONT")
))]

p_blood <- ggplot(tp_blood, aes(x = Methylation, y = SampleGroup, fill = Method)) +
  geom_density_ridges(
    aes(group = interaction(SampleGroup, Coverage2), color = Coverage2),
    scale = 2, alpha = 0.3, from = 0, to = 1
  ) +
  scale_fill_manual(values  = col.vec) +
  scale_color_manual(values = c("None" = "black", "10x" = "red")) +
  scale_y_discrete(labels = rep(paste0("Blood", 1:5), 4)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x     = "Methylation", y = NULL,
    title = "Methylation Density Across Blood Samples\nby Coverage and Method"
  ) +
  guides(color = guide_legend(title = "Coverage Filter")) +
  theme_density(base_size = 25)

ggsave(p_blood,
  filename = file.path(opt$outdir, "Density_Blood.png"),
  height = 12, width = 15, dpi = 300
)

# ---- 5.2 Fibroblast – full density (0 and 10x) ------------------------------

tp_fibro <- buildMethLong(
  data      = fibro,
  samples   = paste0("Fibro", 1:5),
  methods   = METHODS_NO_PACBIO,
  coverages = c(0, 10)
)
tp_fibro <- prep_long(tp_fibro)
tp_fibro[, SampleGroup := factor(SampleGroup, levels = c(
  paste0("Fibro", 1:5, "_WGEC"),  paste0("Fibro", 1:5, "_TWIST"),
  paste0("Fibro", 1:5, "_RRBS"),  paste0("Fibro", 1:5, "_ONT")
))]

p_fibro <- ggplot(tp_fibro, aes(x = Methylation, y = SampleGroup, fill = Method)) +
  geom_density_ridges(
    aes(group = interaction(SampleGroup, Coverage2), color = Coverage2),
    scale = 2, alpha = 0.4, from = 0, to = 1
  ) +
  scale_fill_manual(values  = col.vec) +
  scale_color_manual(values = c("None" = "black", "10x" = "red")) +
  scale_y_discrete(labels = rep(paste0("Fibro", 1:5), 4)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x     = "Methylation", y = NULL,
    title = "Methylation Density Across Fibroblast Samples\nby Coverage and Method"
  ) +
  guides(color = guide_legend(title = "Coverage Filter")) +
  theme_density(base_size = 25)

ggsave(p_fibro,
  filename = file.path(opt$outdir, "Density_Fibro.png"),
  height = 12, width = 15, dpi = 300
)

# ---- 5.3 GIAB – full density (0 and 10x, with PacBio) ----------------------

tp_giab <- buildMethLong(
  data      = giab,
  samples   = c("GIAB1", "GIAB2"),
  methods   = METHODS_PACBIO,
  coverages = c(0, 10)
)
tp_giab <- prep_long(tp_giab)
tp_giab[, SampleGroup := factor(SampleGroup, levels = c(
  "GIAB1_WGEC", "GIAB2_WGEC",
  "GIAB1_TWIST","GIAB2_TWIST",
  "GIAB1_RRBS", "GIAB2_RRBS",
  "GIAB1_PacBio","GIAB2_PacBio",
  "GIAB1_ONT",  "GIAB2_ONT"
))]

p_giab <- ggplot(tp_giab, aes(x = Methylation, y = SampleGroup, fill = Method)) +
  geom_density_ridges(
    aes(group = interaction(SampleGroup, Coverage2), color = Coverage2),
    scale = 2, alpha = 0.4, from = 0, to = 1
  ) +
  scale_fill_manual(values  = col.vec) +
  scale_color_manual(values = c("None" = "black", "10x" = "red")) +
  scale_y_discrete(labels = rep(c("GIAB1", "GIAB2"), 5)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x     = "Methylation", y = NULL,
    title = "Methylation Density Across GIAB Samples\nby Coverage and Method"
  ) +
  guides(color = guide_legend(title = "Coverage Filter")) +
  theme_density(base_size = 25)

ggsave(p_giab,
  filename = file.path(opt$outdir, "Density_GIAB.png"),
  height = 12, width = 15, dpi = 300
)

# ---- 5.4 High-coverage representative samples: Blood3, Fibro4, GIAB2 --------

tp_merged <- rbind(
  buildMethLong(blood, "Blood3", METHODS_NO_PACBIO, c(0, 10, 20, 30, 40)),
  buildMethLong(fibro, "Fibro4", METHODS_NO_PACBIO, c(0, 10, 20, 30, 40)),
  buildMethLong(giab,  "GIAB2",  METHODS_PACBIO,    c(0, 10, 20, 30, 40))
)
tp_merged <- prep_long(tp_merged)
tp_merged[, SampleFacet := fcase(
  grepl("GIAB2",  SampleGroup), "GIAB",
  grepl("Blood3", SampleGroup), "Blood",
  grepl("Fibro4", SampleGroup), "Fibro"
)]
tp_merged[, SampleFacet := factor(SampleFacet, levels = c("GIAB", "Blood", "Fibro"))]

plot_density_window <- function(data, xlim, filename) {
  p <- ggplot(data, aes(x = Methylation, y = fct_rev(SampleGroup), fill = Method)) +
    geom_density_ridges(
      aes(group = interaction(SampleGroup, Coverage2), color = Coverage2),
      linewidth = 1.15, scale = 0.95, alpha = 0.3,
      from = xlim[1], to = xlim[2]
    ) +
    scale_fill_manual(values  = col.vec) +
    scale_color_manual(values = COV_COLORS) +
    coord_cartesian(xlim = xlim) +
    facet_wrap(
      ~SampleFacet, strip.position = "right", scales = "free_y",
      ncol = 1,
      labeller = as_labeller(c("GIAB" = "GIAB", "Blood" = "Blood", "Fibro" = "Fibro"))
    ) +
    labs(
      x     = "Methylation", y = "Density",
      title = "Methylation Density Across Samples\nby Coverage and Method"
    ) +
    guides(color = guide_legend(title = "Coverage Filter")) +
    theme_density()

  ggsave(p, filename = file.path(opt$outdir, filename),
    height = 12, width = 14, dpi = 300)
}

plot_density_window(tp_merged, c(0,   0.3), "Density_High_Cov_0-30.png")
plot_density_window(tp_merged, c(0.7, 1.0), "Density_High_Cov_70-100.png")

# ---- 5.5 GIAB2 high-coverage ONT/TWIST/PacBio only (0–40x) -----------------

tp_giab2 <- buildMethLong(
  data      = giab,
  samples   = c("GIAB1", "GIAB2"),
  methods   = METHODS_HIGH_COV,
  coverages = c(0, 10, 20, 30, 40)
)
tp_giab2 <- prep_long(tp_giab2)
tp_giab2[, Group := factor(SampleGroup, levels = c(
  "GIAB1_TWIST", "GIAB1_PacBio", "GIAB1_ONT",
  "GIAB2_TWIST", "GIAB2_PacBio", "GIAB2_ONT"
))]

plot_giab2_window <- function(xlim, filename) {
  p <- ggplot(
    tp_giab2[Sample == "GIAB2"],
    aes(x = Methylation, y = Group, fill = Method)
  ) +
    geom_density_ridges(
      aes(group = interaction(Group, Coverage2), color = Coverage2),
      linewidth = 1.15, scale = 0.95, alpha = 0.3,
      from = xlim[1], to = xlim[2]
    ) +
    scale_fill_manual(values  = col.vec) +
    scale_color_manual(values = COV_COLORS) +
    scale_y_discrete(labels = rep("GIAB2", 3)) +
    coord_cartesian(xlim = xlim) +
    labs(
      x     = "Methylation", y = "Density",
      title = "Methylation Density Across GIAB2\nby Coverage and Method"
    ) +
    guides(color = guide_legend(title = "Coverage Filter")) +
    theme_density()

  ggsave(p, filename = file.path(opt$outdir, filename),
    height = 12, width = 14, dpi = 300)
}

plot_giab2_window(c(0,   0.3), "Density_GIAB2_0-30.png")
plot_giab2_window(c(0.7, 1.0), "Density_GIAB2_70-100.png")

# ---- 5.6 With-EPIC density at 10x (Blood3, Fibro4, GIAB2) ------------------

if (!is.null(opt$epic_path) && file.exists(opt$epic_path)) {

  cat("[3/4] Building with-EPIC density plot...\n")

  all_mat <- fread(opt$epic_path, header = TRUE, sep = ",", na.strings = "NA")

  # Select only the columns needed: EPIC betas + ONT/WGBS/TWIST/RRBS/PacBio
  # for the three representative samples, then apply joint 10x filter
  keep_meth <- c(
    "Sample3_blood_EPIC", "Sample4_FBK_EPIC", "NA24385_HG002",
    "ONT_Blood3",  "WGBS_Blood3",  "TWIST_Blood3",  "RRBS_Blood3",
    "ONT_Fibro4",  "WGBS_Fibro4",  "TWIST_Fibro4",  "RRBS_Fibro4",
    "ONT_GIAB2",   "WGBS_GIAB2",   "TWIST_GIAB2",   "RRBS_GIAB2",
    "PacBio_GIAB2"
  )
  keep_cov <- c(
    "ONT_cov_Blood3",  "WGBS_cov_Blood3",  "TWIST_cov_Blood3",  "RRBS_cov_Blood3",
    "ONT_cov_Fibro4",  "WGBS_cov_Fibro4",  "TWIST_cov_Fibro4",  "RRBS_cov_Fibro4",
    "ONT_cov_GIAB2",   "WGBS_cov_GIAB2",   "TWIST_cov_GIAB2",   "RRBS_cov_GIAB2",
    "PacBio_cov_GIAB2"
  )

  present_meth <- intersect(keep_meth, colnames(all_mat))
  present_cov  <- intersect(keep_cov,  colnames(all_mat))

  meth_10x <- extractCovDf(all_mat, threshold = 10, cov_cols = present_cov)
  meth_10x <- meth_10x[, ..present_meth]

  # Rename EPIC columns to consistent format
  setnames(meth_10x,
    old = intersect(c("Sample3_blood_EPIC","Sample4_FBK_EPIC","NA24385_HG002"),
                    colnames(meth_10x)),
    new = intersect(c("EPIC_Blood3", "EPIC_Fibro4", "EPIC_GIAB2"),
                    c("EPIC_Blood3", "EPIC_Fibro4", "EPIC_GIAB2")
                    [c("Sample3_blood_EPIC","Sample4_FBK_EPIC","NA24385_HG002")
                      %in% colnames(meth_10x)])
  )

  tp_epic <- melt(meth_10x, measure.vars = colnames(meth_10x),
                  variable.name = "variable", value.name = "value")
  tp_epic <- tp_epic[!is.na(value)]

  # Split variable into Method and Sample
  tp_epic[, c("Method","Sample") := tstrsplit(variable, "_", fixed = TRUE, keep = 1:2)]
  tp_epic[Method == "WGBS", Method := "WGEC"]

  tp_epic[, Group := interaction(Sample, Method, sep = "_")]
  tp_epic[, Group := factor(Group, levels = c(
    "GIAB2_WGEC",  "GIAB2_TWIST", "GIAB2_RRBS",
    "GIAB2_PacBio","GIAB2_ONT",   "GIAB2_EPIC",
    "Fibro4_WGEC", "Fibro4_TWIST","Fibro4_RRBS",
    "Fibro4_ONT",  "Fibro4_EPIC",
    "Blood3_WGEC", "Blood3_TWIST","Blood3_RRBS",
    "Blood3_ONT",  "Blood3_EPIC"
  ))]

  tp_epic[, V2     := factor(Sample, levels = c("GIAB2","Fibro4","Blood3"))]
  tp_epic[, Coverage := factor("10")]

  p_epic <- ggplot(tp_epic, aes(x = value, y = Group, fill = Method)) +
    geom_density_ridges(
      aes(group = interaction(variable, Coverage), color = Coverage),
      linewidth = 1.1, scale = 1.1, alpha = 0.8,
      from = 0, to = 1, panel_scaling = FALSE
    ) +
    scale_fill_manual(values  = col.vec) +
    scale_color_manual(values = c("10" = "black")) +
    scale_y_discrete(labels = c(
      rep("GIAB2", 6), rep("Fibro4", 5), rep("Blood3", 5)
    )) +
    facet_wrap(
      ~V2, scales = "free_y", strip.position = "right", ncol = 1,
      labeller = as_labeller(c(
        "GIAB2" = "GIAB", "Fibro4" = "Fibro", "Blood3" = "Blood"
      ))
    ) +
    labs(
      x     = "Methylation", y = "Density",
      title = "Methylation Density Across Samples\nby Coverage and Method"
    ) +
    guides(color = "none") +
    theme_density()

  ggsave(p_epic,
    filename = file.path(opt$outdir, "Density_EPIC_10x.png"),
    height = 12, width = 14, dpi = 300
  )

} else {
  cat("[3/4] Skipping with-EPIC plot (--epic_path not provided or file not found).\n")
}

cat(sprintf("\n[4/4] Done. Figures written to: %s\n", opt$outdir))
