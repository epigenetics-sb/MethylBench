#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – DMC Annotation (annotatr)
# =============================================================================
# Description:
#   Annotates significant differentially methylated CpGs (DMCs) from both
#   limma and Wilcoxon analyses using the annotatr package against the hg38
#   genome. Generates Figure S14:
#     - Left panel:  raw annotation distribution (with redundancy across
#                    overlapping annotation categories)
#     - Right panel: collapsed to unique annotations per CpG (one annotation
#                    per CpG, priority-based deduplication)
#
#   Additionally exports per-method BED files of significant DMCs for use
#   as annotation input, and saves the full annotated GRanges as RDS.
#
# Annotation categories (hg38):
#   Genes    : 1to5kb, 3UTRs, 5UTRs, CDS, ExonIntronBoundaries, Exons,
#               FirstExon, Intergenic, IntronExonBoundaries, Introns, Promoters
#   lncRNA   : Gencode
#   Enhancers: Fantom5
#   CpG      : Inter, Islands, Shelves, Shores
#
# Input:
#   --limma_combined   Path to all_methods_limma_combined.csv
#                      (from 06a_limma_diff_meth.R)
#   --wilcoxon         Path to Wilcoxon_results.csv
#                      (from 06_differential_methylation.R)
#   --outdir           Output directory for figures
#   --datadir          Output directory for intermediate annotation files
#   --genome           Genome build [default: hg38]
#   --fdr_cutoff       FDR threshold for significance [default: 0.05]
#   --delta_cutoff     Absolute delta-beta threshold [default: 0.1]
#
# Usage:
#   Rscript 07_annotation.R \
#     --limma_combined  data/diff_meth/all_methods_limma_combined.csv \
#     --wilcoxon        data/diff_meth/Wilcoxon_results.csv \
#     --outdir          results/figures/ \
#     --datadir         data/diff_meth/annotation/
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(annotatr)
  library(GenomicRanges)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--limma_combined",
    type    = "character",
    help    = "Path to all_methods_limma_combined.csv [required]",
    metavar = "FILE"
  ),
  make_option("--wilcoxon",
    type    = "character",
    help    = "Path to Wilcoxon_results.csv [required]",
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
    default = "data/diff_meth/annotation/",
    help    = "Output directory for annotation data files [default: data/diff_meth/annotation/]",
    metavar = "DIR"
  ),
  make_option("--genome",
    type    = "character",
    default = "hg38",
    help    = "Genome build for annotatr [default: hg38]",
    metavar = "STRING"
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
    help    = "Absolute delta-beta threshold [default: 0.1]",
    metavar = "FLOAT"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$limma_combined)) stop("ERROR: --limma_combined is required")
if (is.null(opt$wilcoxon))       stop("ERROR: --wilcoxon is required")
if (!file.exists(opt$limma_combined)) stop(paste("File not found:", opt$limma_combined))
if (!file.exists(opt$wilcoxon))       stop(paste("File not found:", opt$wilcoxon))

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

GENOME       <- opt$genome
FDR_CUTOFF   <- opt$fdr_cutoff
DELTA_CUTOFF <- opt$delta_cutoff

VALID_CHRS <- paste0("chr", c(1:22, "X", "Y"))

ANNOT_PRIORITY <- c(
  "hg38_genes_promoters",
  "hg38_genes_firstexons",
  "hg38_genes_5UTRs",
  "hg38_genes_cds",
  "hg38_genes_exons",
  "hg38_genes_exonintronboundaries",
  "hg38_genes_intronexonboundaries",
  "hg38_genes_introns",
  "hg38_genes_3UTRs",
  "hg38_genes_1to5kb",
  "hg38_genes_intergenic",
  "hg38_lncrna_gencode",
  "hg38_enhancers_fantom",
  "hg38_cpg_islands",
  "hg38_cpg_shores",
  "hg38_cpg_shelves",
  "hg38_cpg_inter"
)

col.vec <- get_colors()

cat("[1/5] Loading DMC results...\n")

limma_res <- fread(opt$limma_combined, header = TRUE, sep = ",",
                   na.strings = "NA")
wilcox_res <- fread(opt$wilcoxon, header = TRUE, sep = ",",
                    na.strings = "NA")

# Ensure WGEC label consistency
limma_res[Method  == "WGBS", Method := "WGEC"]
wilcox_res[Method == "WGBS", Method := "WGEC"]

cat(sprintf("  Limma   : %d rows, %d methods\n",
  nrow(limma_res), length(unique(limma_res$Method))))
cat(sprintf("  Wilcoxon: %d rows, %d methods\n",
  nrow(wilcox_res), length(unique(wilcox_res$Method))))

cat("[2/5] Extracting significant DMCs and writing BED files...\n")

parse_cpg_coords <- function(cpg_ids) {
  # Handles both "chr1:12345" and "CpG_12345" formats
  coords <- str_match(cpg_ids, "^(chr[^:]+):([0-9]+)$")
  valid  <- !is.na(coords[, 1])
  dt <- data.table(
    chr   = ifelse(valid, coords[, 2], NA_character_),
    pos   = ifelse(valid, as.integer(coords[, 3]), NA_integer_),
    valid = valid
  )
  return(dt)
}

write_bed <- function(df, path) {
  coords <- parse_cpg_coords(df$CpG)
  bed <- data.table(
    chr        = coords$chr,
    chromStart = coords$pos - 1L,    # 0-based BED
    chromEnd   = coords$pos,
    delta_beta = df$delta_beta
  )
  bed <- bed[coords$valid & chr %in% VALID_CHRS]
  fwrite(bed, path, sep = "\t", col.names = FALSE, quote = FALSE)
  return(invisible(bed))
}

# Limma significant DMCs
limma_sig <- limma_res[
  !is.na(adj.P.Val) &
  adj.P.Val < FDR_CUTOFF &
  abs(delta_beta) > DELTA_CUTOFF
]

cat(sprintf("  Limma significant DMCs: %d (across all methods)\n", nrow(limma_sig)))

for (m in unique(limma_sig$Method)) {
  df_m   <- limma_sig[Method == m]
  prefix <- if (m == "WGEC") "WGBS" else m
  out    <- file.path(opt$datadir, paste0(prefix, "_sig_DMCs.bed"))
  write_bed(df_m, out)
  cat(sprintf("    %s: %d DMCs → %s\n", m, nrow(df_m), basename(out)))
}

# Wilcoxon significant DMCs
wilcox_sig <- wilcox_res[
  !is.na(FDR) &
  FDR < FDR_CUTOFF &
  abs(delta_beta) > DELTA_CUTOFF
]

cat(sprintf("  Wilcoxon significant DMCs: %d\n", nrow(wilcox_sig)))

for (m in unique(wilcox_sig$Method)) {
  df_m   <- wilcox_sig[Method == m]
  prefix <- if (m == "WGEC") "WGBS" else m
  out    <- file.path(opt$datadir, paste0(prefix, "_sig_DMCs_wilcoxon.bed"))
  write_bed(df_m, out)
  cat(sprintf("    %s: %d DMCs → %s\n", m, nrow(df_m), basename(out)))
}

cat("\n[3/5] Building annotatr annotations...\n")

all_annots  <- builtin_annotations()
hg38_annots <- all_annots[grep(GENOME, all_annots)]
cat(sprintf("  Available hg38 annotations: %d\n", length(hg38_annots)))

annot_db <- build_annotations(genome = GENOME, annotations = hg38_annots)

annotateAnnotatr <- function(path, genome = "hg38", annot_db = NULL) {

  df <- fread(path, sep = "\t", header = FALSE)
  df <- df[df$V1 %in% VALID_CHRS]

  colnames(df) <- c("chr", "chromStart", "chromEnd", "Diff")

  dfGR <- GenomicRanges::makeGRangesFromDataFrame(
    as.data.frame(df),
    keep.extra.columns = TRUE,
    seqnames.field     = "chr",
    start.field        = "chromStart",
    end.field          = "chromEnd"
  )

  if (is.null(annot_db)) {
    all_a  <- builtin_annotations()
    hg38_a <- all_a[grep(genome, all_a)]
    annot_db <- build_annotations(genome = genome, annotations = hg38_a)
  }

  result <- annotate_regions(
    regions     = dfGR,
    annotations = annot_db,
    ignore.strand = TRUE,
    quiet         = FALSE
  )

  return(result)
}

cat("\n[4/5] Annotating DMCs...\n")

methods_to_annotate <- unique(limma_sig$Method)
annot_results       <- list()

for (m in methods_to_annotate) {
  prefix <- if (m == "WGEC") "WGBS" else m
  bed_path <- file.path(opt$datadir, paste0(prefix, "_sig_DMCs.bed"))

  if (!file.exists(bed_path)) {
    warning(sprintf("BED file not found for %s – skipping", m))
    next
  }

  cat(sprintf("  Annotating %s...\n", m))

  gr_annot <- annotateAnnotatr(bed_path, genome = GENOME, annot_db = annot_db)

  rds_path <- file.path(opt$datadir, paste0(prefix, "_annotated.rds"))
  saveRDS(gr_annot, rds_path)

  df_annot        <- as.data.frame(gr_annot)
  df_annot$Method <- m
  annot_results[[m]] <- df_annot
}

if (length(annot_results) == 0) {
  stop("No annotation results produced – check BED files and genome build.")
}

annot_combined <- bind_rows(annot_results)
annot_combined$Method[annot_combined$Method == "WGBS"] <- "WGEC"

fwrite(as.data.table(annot_combined),
  file.path(opt$datadir, "all_methods_annotated.csv"),
  sep = ",", quote = FALSE
)

cat("\n[5/5] Generating annotation figures...\n")

# Clean up annotation type labels for plotting
clean_annot_label <- function(x) {
  x <- sub(paste0(GENOME, "_"), "", x)
  x <- sub("genes_", "", x)
  x <- sub("cpg_", "CpG ", x)
  x <- sub("lncrna_gencode", "lncRNA (Gencode)", x)
  x <- sub("enhancers_fantom", "Enhancers (Fantom)", x)
  x <- tools::toTitleCase(gsub("_", " ", x))
  return(x)
}

annot_plot <- annot_combined %>%
  filter(!is.na(annot.type)) %>%
  mutate(annot_label = clean_annot_label(annot.type))

# Compute proportions per method
annot_freq <- annot_plot %>%
  group_by(Method, annot_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Define consistent category order and color palette
annot_order <- unique(annot_freq$annot_label[
  order(match(
    sub(paste0(GENOME, "_"), "", unique(annot_plot$annot.type)),
    sub(paste0(GENOME, "_"), "", ANNOT_PRIORITY)
  ))
])

annot_freq$annot_label <- factor(annot_freq$annot_label, levels = rev(annot_order))

# Color palette for annotation categories (grouped by type)
n_annot     <- length(unique(annot_freq$annot_label))
annot_colors <- setNames(
  colorRampPalette(c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#8CD17D","#86BCB6","#499894","#F1CE63","#D4A6C8",
    "#FFBE7D","#A0CBE8"
  ))(n_annot),
  levels(annot_freq$annot_label)
)

p_raw <- ggplot(annot_freq, aes(x = Method, y = pct, fill = annot_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = annot_colors, name = "Annotation") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Annotation of significant DMCs",
    x     = NULL,
    y     = "Proportion (%)"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 18),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16),
    text         = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

ggsave(p_raw,
  filename = file.path(opt$outdir, "Annotation_DMCs_raw.png"),
  height = 12, width = 14, dpi = 300
)
cat("  Saved: Annotation_DMCs_raw.png\n")

annot_collapsed <- annot_plot %>%
  mutate(
    priority = match(annot.type,
                     sub(paste0("^"), paste0(GENOME, "_"), ANNOT_PRIORITY))
  ) %>%
  group_by(Method, seqnames, start, end) %>%
  slice_min(priority, n = 1, with_ties = FALSE) %>%
  ungroup()

annot_freq_collapsed <- annot_collapsed %>%
  group_by(Method, annot_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Method) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

annot_freq_collapsed$annot_label <- factor(
  annot_freq_collapsed$annot_label,
  levels = levels(annot_freq$annot_label)
)

p_collapsed <- ggplot(
  annot_freq_collapsed,
  aes(x = Method, y = pct, fill = annot_label)
) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = annot_colors, name = "Annotation") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Annotation (unique per CpG)",
    x     = NULL,
    y     = "Proportion (%)"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 18),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16),
    text         = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

ggsave(p_collapsed,
  filename = file.path(opt$outdir, "Annotation_DMCs_collapsed.png"),
  height = 12, width = 14, dpi = 300
)
cat("  Saved: Annotation_DMCs_collapsed.png\n")

annot_freq$Panel         <- "Raw (with redundancy)"
annot_freq_collapsed$Panel <- "Collapsed (unique per CpG)"

annot_combined_plot <- bind_rows(annot_freq, annot_freq_collapsed)
annot_combined_plot$Panel <- factor(
  annot_combined_plot$Panel,
  levels = c("Raw (with redundancy)", "Collapsed (unique per CpG)")
)
annot_combined_plot$annot_label <- factor(
  annot_combined_plot$annot_label,
  levels = levels(annot_freq$annot_label)
)

p_s14 <- ggplot(
  annot_combined_plot,
  aes(x = Method, y = pct, fill = annot_label)
) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = annot_colors, name = "Annotation") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  facet_wrap(~Panel, ncol = 2) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "Proportion (%)"
  ) +
  theme_bw() +
  theme(
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16),
    text         = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

ggsave(p_s14,
  filename = file.path(opt$outdir, "Annotation_DMCs_S14.png"),
  height = 12, width = 20, dpi = 300
)
