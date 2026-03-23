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

to.plot <- melt(
  cpg.stats[, c(1, 3, 4, 5, 6, 7)],
  id.vars = "Sample"
)
to.plot[, variable := as.character(variable)]
to.plot[, Coverage := substr(variable, 5, nchar(variable))]

x.labels <- c("None", "10x", "15x", "20x", "40x")

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

ont_bases <- stats[
  Method == "ONT" & Sampleset %in% c("Blood", "Fibroblast"),
  .(Sample, Method, qc_bases_sequenced = Mean_Cov)
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
