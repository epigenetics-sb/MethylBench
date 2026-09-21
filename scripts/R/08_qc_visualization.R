#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – QC Visualization
# =============================================================================
# Description:
#   Generates QC figures from cpg_stats and qc_stats summary tables.
#   Produces: Overlapping CpGs, Mean Methylation, Mean Coverage,
#             ONT Genomic Coverage, ONT Bases Sequenced,
#             Mean Coverage per Method, Read Length, Unique Reads
#
# Input:
#   --cpg_stats   cpg_stats_methylbench.tab  (from 01_generate_cpg_stats.R)
#   --qc_stats    qc_stats_methylbench.tab   (manually compiled)
#   --outdir      Output directory for figures
#
# Usage:
#   Rscript 02_qc_visualization.R \
#     --cpg_stats data/stats/cpg_stats_methylbench.tab \
#     --qc_stats  data/stats/qc_stats_methylbench.tab \
#     --outdir    results/figures/
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
  make_option("--cpg_stats",
    type    = "character",
    help    = "Path to cpg_stats_methylbench.tab [required]",
    metavar = "FILE"
  ),
  make_option("--qc_stats",
    type    = "character",
    help    = "Path to qc_stats_methylbench.tab [required]",
    metavar = "FILE"
  ),
  make_option("--outdir",
    type    = "character",
    default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"
  ),
  make_option("--all_path",
    type    = "character",
    default = NULL,
    help    = "Optional: path to the merged CpG-level matrix (e.g. ALL.csv / ALL_without_EPIC.csv) with <Method>_cov_<Sample> columns, for the per-CpG coverage uniformity figure (Reviewer 1, Minor #2). Sample/tissue (Blood/Fibroblast/GIAB) are derived from the <Method>_cov_<Sample> column names themselves, same convention as the rest of the pipeline -- no samplesheet needed. Skipped entirely if not provided.",
    metavar = "FILE"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cpg_stats)) stop("ERROR: --cpg_stats is required")
if (is.null(opt$qc_stats))  stop("ERROR: --qc_stats is required")
if (!file.exists(opt$cpg_stats)) stop(paste("File not found:", opt$cpg_stats))
if (!file.exists(opt$qc_stats))  stop(paste("File not found:", opt$qc_stats))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cpg.stats <- fread(opt$cpg_stats, header = TRUE, sep = "\t")
stats      <- fread(opt$qc_stats,  header = TRUE, sep = "\t",
                    na.strings = c("", "NA"))

col.vec <- get_colors()

stats[, Sampleset_numbers := ifelse(
  Sampleset == "GIAB",
  paste0(Sampleset, "\n(n=2)"),
  paste0(Sampleset, "\n(n=5)")
)]

SAMPLESET_ORDER <- c("Blood\n(n=5)", "Fibroblast\n(n=5)", "GIAB\n(n=2)")
stats[, Sampleset_numbers := factor(Sampleset_numbers, levels = SAMPLESET_ORDER)]

# ---- 3.1 Overlapping CpGs per coverage filter -------------------------------
cat("[1/10] Plotting Overlapping CpGs...\n")

to.plot <- melt(
  cpg.stats[, c(1, 3, 4, 5, 6, 7, 8)],
  id.vars = "Sample"
)
to.plot[, variable := as.character(variable)]
to.plot[, Coverage := substr(variable, 5, nchar(variable))]

x.labels <- c("None", "10x", "15x", "20x", "30x", "40x")

p_cpg <- ggplot(to.plot, aes(x = Coverage, y = value, color = Sample, fill = Sample)) +
  geom_point(size = 4) +
  geom_line(aes(group = Sample)) +
  scale_y_log10(
    labels = c("0", "1000", "10000", "0.1M", "1M", "10M"),
    breaks = c(0, 1e3, 1e4, 1e5, 1e6, 1e7)
  ) +
  scale_color_manual(values = C12) +
  scale_x_discrete(labels = x.labels) +
  labs(
    title = "Number of overlapping CpGs per Coverage filter",
    x     = "Coverage Filter",
    y     = "No. of CpGs"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 22),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22)
  )

ggsave(p_cpg,
  filename = file.path(opt$outdir, "Overlapping_CpGs.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.2 Mean Methylation ------------
cat("[2/10] Plotting Mean Methylation...\n")

p_meth_focused <- ggplot(
  stats,
  aes(x = Method, y = Mean_meth_overlapped, color = Method, fill = Method)
) +
  geom_point(shape = 21, size = 3.5) +
  geom_line(aes(group = Sample), color = "black") +
  facet_wrap(~Sampleset_numbers) +
  scale_fill_manual(values  = col.vec) +
  scale_color_manual(values = col.vec) +
  labs(
    title = "Mean Methylation",
    x     = "",
    y     = "Methylation[%]"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_meth_focused,
  filename = file.path(opt$outdir, "Mean_Methylation_Focused.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.3 Mean Methylation with coverage filter comparison -------------------
cat("[3/10] Plotting Mean Methylation (coverage comparison)...\n")

uno <- cbind(stats[, .(Sample, Method, Mean_meth = Mean_meth_overlapped,
                        Sampleset_numbers)],
             Coverage = "No Filter")
dos <- cbind(stats[, .(Sample, Method, Mean_meth = Mean_meth_10x,
                        Sampleset_numbers)],
             Coverage = "10x Filter")
to.plot <- rbind(uno, dos)
to.plot[, Coverage := factor(Coverage, levels = c("No Filter", "10x Filter"))]
to.plot[, Line     := paste0(Sample, Coverage)]

p_meth_cov <- ggplot(
  to.plot,
  aes(x = Method, y = Mean_meth, color = Method, fill = Method, shape = Coverage)
) +
  geom_point(size = 4) +
  geom_line(aes(group = Line, color = "black")) +
  facet_grid(Coverage ~ Sampleset_numbers) +
  scale_fill_manual(values  = col.vec) +
  scale_color_manual(values = col.vec) +
  ylim(c(0, 1)) +
  labs(
    title = "Mean Methylation",
    x     = "",
    y     = "Methylation Value"
  ) +
  guides(shape = guide_legend(title = "Coverage\nFilter")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_meth_cov,
  filename = file.path(opt$outdir, "Mean_Methylation_Coverage.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.4 Mean CpG Coverage per sample and method ----------------------------
cat("[4/10] Plotting Mean CpG Coverage...\n")

p_cov <- ggplot(stats, aes(x = Sample, y = Mean_Cov, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Mean_Cov, 0), y = 20), size = 6) +
  geom_hline(
    yintercept = c(10, 20, 30, 40, 50, 60),
    linetype   = "dashed",
    alpha      = 0.3
  ) +
  facet_grid(
    rows   = vars(Method),
    cols   = vars(Sampleset_numbers),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = col.vec) +
  labs(
    title = "Mean CpG Coverage",
    x     = "Sample",
    y     = "Mean CpG Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_cov,
  filename = file.path(opt$outdir, "Mean_Coverage.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.5 ONT Mean Genomic Coverage ------------------------------------------
cat("[5/10] Plotting ONT Mean Genomic Coverage...\n")

p_ont_genomic <- ggplot(
  stats[Method == "ONT"],
  aes(x = Sample, y = Mean_Genome_Cov, fill = Method)
) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Mean_Genome_Cov, 0), y = 20), size = 6) +
  geom_hline(
    yintercept = c(10, 20, 30, 40, 50, 60),
    linetype   = "dashed",
    alpha      = 0.3
  ) +
  scale_fill_manual(values = col.vec) +
  labs(
    title = "Mean Genomic Coverage for ONT",
    x     = "Sample",
    y     = "Mean Genomic Coverage"
  ) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_ont_genomic,
  filename = file.path(opt$outdir, "Mean_Genomic_Coverage.png"),
  height = 10, width = 15, dpi = 300
)

# ---- 3.6 ONT Total Bases Sequenced ------------------------------------------
cat("[6/10] Plotting ONT Total Bases Sequenced...\n")

ont_bases <- stats[
  Method == "ONT" & Sampleset %in% c("Blood", "Fibroblast"),
  .(Sample, Method, qc_bases_sequenced = Mean_Genome_Cov)
]

p_ont_bases <- ggplot(
  ont_bases,
  aes(x = Sample, y = qc_bases_sequenced, fill = Method)
) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(qc_bases_sequenced, 0), y = 50), size = 6) +
  geom_hline(
    yintercept = seq(15, 120, by = 15),
    linetype   = "dashed",
    alpha      = 0.3
  ) +
  scale_fill_manual(values = col.vec) +
  labs(
    title = "Total Bases Sequenced (Gigabases) for ONT",
    x     = "Sample",
    y     = "Total Bases Sequenced (Gigabases)"
  ) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_ont_bases,
  filename = file.path(opt$outdir, "QC_Bases_Sequenced.png"),
  height = 10, width = 15, dpi = 300
)

# ---- 3.7 Mean CpG Coverage per method (averaged over all samples) -----------
cat("[7/10] Plotting Mean CpG Coverage per Method...\n")

df_means <- stats[, .(Mean_Cov_overall = mean(Mean_Cov, na.rm = TRUE)), by = Method]

p_mean_mean <- ggplot(df_means, aes(x = Method, y = Mean_Cov_overall, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Mean_Cov_overall, 0), y = 20), size = 6) +
  geom_hline(
    yintercept = c(10, 20, 30, 40, 50, 60),
    linetype   = "dashed",
    alpha      = 0.3
  ) +
  scale_fill_manual(values = col.vec) +
  labs(
    title = "Mean CpG Coverage over all samples",
    x     = "Method",
    y     = "Mean CpG Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    plot.title  = element_text(hjust = 0.5)
  )

ggsave(p_mean_mean,
  filename = file.path(opt$outdir, "Mean_Mean_Coverage.png"),
  height = 10, width = 15, dpi = 300
)

# ---- 3.8 Mean Read Length ---------------------------------------------------
cat("[8/10] Plotting Mean Read Length...\n")

p_readlen <- ggplot(
  stats,
  aes(x = Method, y = Mean_readlength_passed, fill = Sample)
) +
  geom_point(size = 4, shape = 21) +
  facet_wrap(~Sampleset) +
  scale_fill_manual(values = C12) +
  scale_y_log10(
    labels = c("0", "100", "1000", "10000", "20000"),
    breaks = c(0, 100, 1000, 10000, 20000)
  ) +
  labs(
    title = "Mean read length",
    x     = "",
    y     = "Mean read length of unique alignments"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 22),
    axis.text   = element_text(size = 22),
    axis.title  = element_text(size = 22),
    text        = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(p_readlen,
  filename = file.path(opt$outdir, "Mean_Readlength.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.9 Unique Reads (excluding PacBio) ------------------------------------
cat("[9/10] Plotting Unique Reads...\n")

p_unique <- ggplot(
  stats[Method != "PacBio"],
  aes(x = Method, y = Unique_alignments, fill = Sample)
) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = C12) +
  scale_y_log10(
    labels = c("0", "10M", "50M", "100M", "500M"),
    breaks = c(0, 1e7, 5e7, 1e8, 5e8)
  ) +
  labs(
    title = "Unique Reads",
    x     = "",
    y     = "No. of unique reads"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text  = element_text(size = 22),
    axis.title = element_text(size = 22),
    text       = element_text(size = 22)
  )

ggsave(p_unique,
  filename = file.path(opt$outdir, "Unique_Reads.png"),
  height = 10, width = 10, dpi = 300
)

# ---- 3.10 Coverage uniformity across methods (Reviewer 1, Minor #2) --------
# Addresses: "Short-read methods generally achieved a more uniform CpG
# coverage" (Section 2.1) was not directly supported by Figure 3, which
# only shows mean coverage per sample/method, not the per-CpG distribution.
# Optional -- only runs if --all_path was provided.

#' Produce per-CpG coverage density + ECDF plots (normalized, colored by
#' method, faceted by tissue) plus a per-method/tissue coefficient-of-
#' variation (CV) summary table, as a quantitative uniformity metric.
#'
#' Tissue (Blood/Fibroblast/GIAB) and the sample list are derived directly
#' from the <Method>_cov_<Sample> column names in all_dt -- the same
#' Blood1..5/Fibro1..5/GIAB1..2 naming convention used throughout the rest
#' of the pipeline -- rather than requiring a separate samplesheet file.
#'
#' @param all_dt merged CpG-level data.table with <Method>_cov_<Sample> cols.
#' @param methods_vec named character vector of method prefixes, as used
#'   elsewhere in the pipeline, e.g.
#'   c(ONT="ONT", RRBS="RRBS", WGEC="WGEC", TWIST="TWIST", PacBio="PacBio")
#' @param out_dir output directory for the figures/summary table.
plot_coverage_uniformity <- function(all_dt,
                                      methods_vec = c(ONT = "ONT", RRBS = "RRBS",
                                                       WGEC = "WGEC", TWIST = "TWIST",
                                                       PacBio = "PacBio"),
                                      out_dir = "results/figures/") {

  # ---- Derive samples + tissue from <Method>_cov_<Sample> column names -----
  cov_cols <- grep("_cov_", colnames(all_dt), value = TRUE)
  samples  <- unique(sub(".*_cov_", "", cov_cols))
  if (length(samples) == 0) {
    warning("plot_coverage_uniformity: no <Method>_cov_<Sample> columns found in all_dt -- skipping.")
    return(invisible(NULL))
  }

  tissue_of <- function(s) {
    if (grepl("^Blood", s)) return("Blood")
    if (grepl("^Fibro", s)) return("Fibroblast")
    if (grepl("^GIAB",  s)) return("GIAB")
    NA_character_
  }

  # ---- Long-format per-CpG coverage, via helpers.R's buildCovLong() --------
  cov_long <- buildCovLong(all_dt, samples = samples, methods = methods_vec)

  cov_long[, Tissue := vapply(Sample, tissue_of, character(1))]
  tissue_levels <- intersect(c("Blood", "Fibroblast", "GIAB"), unique(cov_long$Tissue))
  cov_long[, Tissue := factor(Tissue, levels = tissue_levels)]

  # ---- Normalize per Sample x Method: methods differ hugely in mean depth
  #      (e.g. WGEC 28-30M CpGs vs. RRBS/TWIST 1-4M), so per-CpG coverage
  #      must be compared relative to each sample/method's own mean to
  #      assess UNIFORMITY rather than raw depth -----------------------------
  cov_long[, mean_cov_sample := mean(Coverage), by = .(Sample, Method)]
  cov_long[, NormCoverage    := Coverage / mean_cov_sample]

  # ---- Quantitative uniformity metric: CV per Sample x Method, averaged
  #      per Tissue x Method (lower CV = more uniform per-CpG coverage) -----
  cv_summary <- cov_long[
    , .(cv = sd(Coverage) / mean(Coverage)),
    by = .(Sample, Tissue, Method)
  ][
    , .(mean_CV = mean(cv), sd_CV = sd(cv), n_samples = .N),
    by = .(Tissue, Method)
  ][order(Tissue, Method)]

  fwrite(cv_summary, file.path(out_dir, "Coverage_Uniformity_CV_summary.tab"), sep = "\t")

  method_colors <- get_colors()
  breaks_x <- c(0.1, 0.5, 1, 2, 8)
  labels_x <- c("0.1x", "0.5x", "1x", "2x", "8x")

  p_density <- ggplot(cov_long, aes(x = NormCoverage, color = Method, fill = Method)) +
    geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2) +
    scale_x_log10(breaks = breaks_x, labels = labels_x) +
    coord_cartesian(xlim = c(0.03, 10)) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    facet_wrap(~ Tissue, nrow = 1) +
    labs(
      x = "Per-CpG coverage (normalized to sample mean, log10 scale)",
      y = "Density",
      color = "Method", fill = "Method",
      title = "Per-CpG coverage uniformity across methylation profiling methods"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom"
    )

  p_ecdf <- ggplot(cov_long, aes(x = NormCoverage, color = Method)) +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    scale_x_log10(breaks = breaks_x, labels = labels_x) +
    coord_cartesian(xlim = c(0.03, 10)) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    facet_wrap(~ Tissue, nrow = 1) +
    labs(
      x = "Per-CpG coverage (normalized to sample mean, log10 scale)",
      y = "Cumulative fraction of CpGs",
      color = "Method",
      title = "ECDF of per-CpG coverage across methylation profiling methods"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom"
    )

  ggsave(file.path(out_dir, "Coverage_Uniformity_Density.png"),
         p_density, width = 13, height = 4.5, units = "in", dpi = 300)
  ggsave(file.path(out_dir, "Coverage_Uniformity_ECDF.png"),
         p_ecdf, width = 13, height = 4.5, units = "in", dpi = 300)

  invisible(list(density = p_density, ecdf = p_ecdf,
                  cv_summary = cv_summary, long = cov_long))
}

if (!is.null(opt$all_path)) {
  if (!file.exists(opt$all_path)) stop(paste("File not found (--all_path):", opt$all_path))
  cat("[10/10] Plotting Coverage Uniformity...\n")
  all_dt <- fread(opt$all_path, header = TRUE, sep = ",", na.strings = "NA")
  plot_coverage_uniformity(all_dt, out_dir = opt$outdir)
} else {
  cat("[10/10] Coverage Uniformity skipped (--all_path not provided)\n")
}
