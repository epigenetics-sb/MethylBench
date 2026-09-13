#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Background-Corrected Annotation Enrichment (Gene-centric vs.
# CpG-structural, assessed separately)
# =============================================================================
# Description:
#   Addresses Reviewer 1, Major Comment:
#   "The annotation analysis does not fully support the enrichment claims.
#    Figures 9C and Supplementary Figure 9 report annotation proportions
#    among called DMRs without comparison to the corresponding tested CpG/
#    region background. In addition, in Supplementary Figure 9, the apparent
#    disappearance of CpG-island-related categories after collapsing
#    annotations to unique features may depend on the annotation priority
#    hierarchy. It may therefore be more informative to assess CpG structural
#    and gene-centric annotations separately."
#
#   This script replaces the single, mixed priority hierarchy used previously
#   (13_annotation.R: ANNOT_PRIORITY, which ranks gene features above CpG-
#   structural features) with TWO independent, internally mutually-exclusive
#   hierarchies:
#
#     Gene-centric  (each CpG assigned to exactly one class):
#       Promoter > 5'UTR > CDS/Exon > Intron > 3'UTR > 1-5kb upstream >
#       Intergenic
#
#     CpG-structural (each CpG assigned to exactly one class):
#       CpG Island > CpG Shore > CpG Shelf > Open Sea
#
#   Because these two hierarchies do not compete with one another, a CpG can
#   (and does) carry one gene-centric AND one CpG-structural label
#   simultaneously; no CGI signal is discarded by gene-centric prioritization
#   and vice versa.
#
#   Critically, both DMC/DMR annotation profiles and the annotation profile
#   of the FULL TESTED BACKGROUND (the Tier-1/Tier-2 consensus CpG set that
#   entered DSS/limma testing, i.e. the set at risk of being called, not the
#   full genome) are computed. Enrichment is then reported as:
#     - the DMC/DMR proportion in each category,
#     - the background proportion in each category,
#     - log2 fold enrichment = log2(DMC%/background%),
#     - a one-sided Fisher's exact test / odds ratio with 95% CI per
#       category (DMC/DMR-in-category vs. not, among tested vs. not tested).
#
# Input:
#   --sig_bed_dir     Directory containing per-method significant DMC BED
#                      files (as produced by 13_annotation.R / 14 or 15).
#   --background_tsv  TSV with the tested CpG background, one row per CpG,
#                      columns: Chr, Pos (e.g. Tier1_consensus_CpGs.tsv /
#                      Tier2_consensus_CpGs.tsv from 14_DMR_DSS_analysis.R,
#                      or the DMR-testing background from 15_DMR_DMRcate).
#   --outdir          Output directory for figures
#   --datadir         Output directory for tables
#   --genome          Genome build [default: hg38]
#   --level           "CpG" or "DMR" (only affects labeling) [default: CpG]
#
# Output:
#   - annotation_enrichment_genecentric.tsv
#   - annotation_enrichment_cpgstructural.tsv
#   - Fig_annotation_enrichment_genecentric.png
#   - Fig_annotation_enrichment_cpgstructural.png
#
# Usage:
#   Rscript scripts/R/17_annotation_enrichment_background.R \
#     --sig_bed_dir     data/diff_meth/annotation/ \
#     --background_tsv  results/dmr_dss/Tier1_consensus_CpGs.tsv \
#     --outdir          results/figures/ \
#     --datadir         data/diff_meth/annotation_enrichment/
#
# Author: MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(annotatr)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--sig_bed_dir",    type = "character", metavar = "DIR"),
  make_option("--background_tsv", type = "character", metavar = "FILE"),
  make_option("--outdir",  type = "character", default = "results/figures/"),
  make_option("--datadir", type = "character", default = "data/diff_meth/annotation_enrichment/"),
  make_option("--genome",  type = "character", default = "hg38"),
  make_option("--level",   type = "character", default = "CpG")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$sig_bed_dir))    stop("ERROR: --sig_bed_dir is required")
if (is.null(opt$background_tsv)) stop("ERROR: --background_tsv is required")

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

GENOME <- opt$genome
VALID_CHRS <- paste0("chr", c(1:22, "X", "Y"))

# -------------------------------------------------------------------------
# Two independent, mutually-exclusive annotation hierarchies
# -------------------------------------------------------------------------
GENE_PRIORITY <- c(
  paste0(GENOME, "_genes_promoters"),
  paste0(GENOME, "_genes_5UTRs"),
  paste0(GENOME, "_genes_cds"),
  paste0(GENOME, "_genes_firstexons"),
  paste0(GENOME, "_genes_exons"),
  paste0(GENOME, "_genes_exonintronboundaries"),
  paste0(GENOME, "_genes_intronexonboundaries"),
  paste0(GENOME, "_genes_introns"),
  paste0(GENOME, "_genes_3UTRs"),
  paste0(GENOME, "_genes_1to5kb"),
  paste0(GENOME, "_genes_intergenic")
)

CPG_PRIORITY <- c(
  paste0(GENOME, "_cpg_islands"),
  paste0(GENOME, "_cpg_shores"),
  paste0(GENOME, "_cpg_shelves"),
  paste0(GENOME, "_cpg_inter")
)

clean_label <- function(x, genome) {
  x <- sub(paste0(genome, "_"), "", x)
  x <- sub("genes_", "", x)
  x <- sub("cpg_", "CpG ", x)
  x <- tools::toTitleCase(gsub("_", " ", x))
  x
}

cat("[1/6] Building annotatr databases (gene-centric & CpG-structural)...\n")
gene_db <- build_annotations(genome = GENOME, annotations = GENE_PRIORITY)
cpg_db  <- build_annotations(genome = GENOME, annotations = CPG_PRIORITY)

annotate_exclusive <- function(gr, annot_db, priority) {
  res <- annotate_regions(regions = gr, annotations = annot_db,
                           ignore.strand = TRUE, quiet = TRUE)
  df  <- as.data.frame(res)
  if (nrow(df) == 0) return(df)
  df$priority <- match(df$annot.type, priority)
  df <- df %>%
    group_by(seqnames, start, end) %>%
    slice_min(priority, n = 1, with_ties = FALSE) %>%
    ungroup()
  df
}

read_bed_as_gr <- function(path) {
  df <- fread(path, sep = "\t", header = FALSE)
  colnames(df) <- c("chr", "chromStart", "chromEnd", "Diff")[seq_len(ncol(df))]
  df <- df[df$chr %in% VALID_CHRS]
  GenomicRanges::makeGRangesFromDataFrame(
    as.data.frame(df), keep.extra.columns = TRUE,
    seqnames.field = "chr", start.field = "chromStart", end.field = "chromEnd"
  )
}

# -------------------------------------------------------------------------
# 2. Annotate the tested BACKGROUND once per hierarchy
# -------------------------------------------------------------------------
cat("[2/6] Annotating tested background CpG set...\n")

bg <- fread(opt$background_tsv)
if (all(c("Chr", "Pos") %in% colnames(bg))) {
  bg_gr <- GRanges(seqnames = bg$Chr, ranges = IRanges(bg$Pos, bg$Pos))
} else if ("cpg_id" %in% colnames(bg)) {
  parts <- str_match(bg$cpg_id, "^(chr[^:]+):([0-9]+)$")
  bg_gr <- GRanges(seqnames = parts[, 2], ranges = IRanges(as.integer(parts[, 3]),
                                                            as.integer(parts[, 3])))
} else {
  stop("--background_tsv must contain either Chr/Pos or cpg_id columns")
}
bg_gr <- bg_gr[as.character(seqnames(bg_gr)) %in% VALID_CHRS]

bg_gene <- annotate_exclusive(bg_gr, gene_db, GENE_PRIORITY)
bg_cpg  <- annotate_exclusive(bg_gr, cpg_db,  CPG_PRIORITY)

bg_gene_freq <- bg_gene %>% count(annot.type, name = "n") %>%
  mutate(pct = n / sum(n) * 100, Set = "Background (tested)")
bg_cpg_freq  <- bg_cpg  %>% count(annot.type, name = "n") %>%
  mutate(pct = n / sum(n) * 100, Set = "Background (tested)")

cat(sprintf("  Background CpGs annotated: %d (gene-centric), %d (CpG-structural)\n",
            nrow(bg_gene), nrow(bg_cpg)))

# -------------------------------------------------------------------------
# 3. Annotate each method's significant DMCs/DMRs, both hierarchies
# -------------------------------------------------------------------------
cat("[3/6] Annotating per-method significant sites...\n")

bed_files <- list.files(opt$sig_bed_dir, pattern = "_sig_DMCs\\.bed$|_sig_DMRs\\.bed$",
                         full.names = TRUE)
if (length(bed_files) == 0) stop("No significant-site BED files found in --sig_bed_dir")

gene_list <- list()
cpg_list  <- list()

for (f in bed_files) {
  method <- gsub("_sig_DMCs\\.bed$|_sig_DMRs\\.bed$", "", basename(f))
  method <- ifelse(method == "WGBS", "WGEC", method)
  gr <- read_bed_as_gr(f)
  if (length(gr) == 0) next

  gene_ann <- annotate_exclusive(gr, gene_db, GENE_PRIORITY)
  cpg_ann  <- annotate_exclusive(gr, cpg_db,  CPG_PRIORITY)

  gene_list[[method]] <- gene_ann %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = method)
  cpg_list[[method]]  <- cpg_ann %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = method)
}

gene_freq <- bind_rows(gene_list, .id = NULL)
cpg_freq  <- bind_rows(cpg_list,  .id = NULL)

# -------------------------------------------------------------------------
# 4. Enrichment statistics: log2FE + Fisher's exact test vs. background
# -------------------------------------------------------------------------
cat("[4/6] Computing enrichment vs. background (log2FE + Fisher's exact test)...\n")

compute_enrichment <- function(sig_freq, bg_freq, sig_totals, bg_total) {
  categories <- union(sig_freq$annot.type, bg_freq$annot.type)
  out <- list()
  for (m in unique(sig_freq$Set)) {
    n_sig_total <- sig_totals[[m]]
    for (cat_i in categories) {
      n_sig_cat <- sig_freq$n[sig_freq$Set == m & sig_freq$annot.type == cat_i]
      n_sig_cat <- if (length(n_sig_cat) == 0) 0 else n_sig_cat
      n_bg_cat  <- bg_freq$n[bg_freq$annot.type == cat_i]
      n_bg_cat  <- if (length(n_bg_cat) == 0) 0 else n_bg_cat

      tab <- matrix(c(
        n_sig_cat, n_sig_total - n_sig_cat,
        n_bg_cat,  bg_total - n_bg_cat
      ), nrow = 2)

      ft <- tryCatch(fisher.test(tab), error = function(e) NULL)

      pct_sig <- 100 * n_sig_cat / n_sig_total
      pct_bg  <- 100 * n_bg_cat  / bg_total
      log2fe  <- log2((pct_sig + 1e-6) / (pct_bg + 1e-6))

      out[[length(out) + 1]] <- data.table(
        Method = m, Category = clean_label(cat_i, GENOME),
        pct_sig = pct_sig, pct_background = pct_bg, log2FE = log2fe,
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        ci_low  = if (!is.null(ft)) ft$conf.int[1] else NA_real_,
        ci_high = if (!is.null(ft)) ft$conf.int[2] else NA_real_,
        p_value = if (!is.null(ft)) ft$p.value else NA_real_
      )
    }
  }
  dt <- rbindlist(out)
  dt[, FDR := p.adjust(p_value, method = "BH")]
  dt
}

gene_totals <- gene_freq %>% group_by(Set) %>% summarise(n = sum(n)) %>% deframe()
cpg_totals  <- cpg_freq  %>% group_by(Set) %>% summarise(n = sum(n)) %>% deframe()

gene_enrich <- compute_enrichment(gene_freq, bg_gene_freq, gene_totals, sum(bg_gene_freq$n))
cpg_enrich  <- compute_enrichment(cpg_freq,  bg_cpg_freq,  cpg_totals,  sum(bg_cpg_freq$n))

fwrite(gene_enrich, file.path(opt$datadir, "annotation_enrichment_genecentric.tsv"), sep = "\t")
fwrite(cpg_enrich,  file.path(opt$datadir, "annotation_enrichment_cpgstructural.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# 5. Figures: proportion (sig vs. background) + log2FE, per hierarchy
# -------------------------------------------------------------------------
cat("[5/6] Generating figures...\n")

make_prop_plot <- function(sig_freq, bg_freq, title) {
  sig_freq2 <- sig_freq %>% mutate(annot_label = clean_label(annot.type, GENOME))
  bg_freq2  <- bg_freq  %>% mutate(annot_label = clean_label(annot.type, GENOME))
  combined <- bind_rows(sig_freq2, bg_freq2)
  ggplot(combined, aes(x = Set, y = pct, fill = annot_label)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(title = title, x = NULL, y = "Proportion (%)", fill = "Annotation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_gene <- make_prop_plot(gene_freq, bg_gene_freq,
  sprintf("Gene-centric annotation of significant %ss vs. tested background", opt$level))
p_cpg  <- make_prop_plot(cpg_freq, bg_cpg_freq,
  sprintf("CpG-structural annotation of significant %ss vs. tested background", opt$level))

ggsave(file.path(opt$outdir, "Fig_annotation_proportion_genecentric.png"),
       p_gene, width = 10, height = 7, dpi = 300)
ggsave(file.path(opt$outdir, "Fig_annotation_proportion_cpgstructural.png"),
       p_cpg, width = 10, height = 7, dpi = 300)

p_gene_fe <- ggplot(gene_enrich, aes(x = Category, y = log2FE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Gene-centric enrichment vs. tested background",
       x = NULL, y = expression(log[2]~"fold enrichment (sig. / background)")) +
  theme_bw()

p_cpg_fe <- ggplot(cpg_enrich, aes(x = Category, y = log2FE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "CpG-structural enrichment vs. tested background",
       x = NULL, y = expression(log[2]~"fold enrichment (sig. / background)")) +
  theme_bw()

ggsave(file.path(opt$outdir, "Fig_annotation_enrichment_genecentric.png"),
       p_gene_fe, width = 9, height = 6, dpi = 300)
ggsave(file.path(opt$outdir, "Fig_annotation_enrichment_cpgstructural.png"),
       p_cpg_fe, width = 8, height = 5, dpi = 300)

cat("[6/6] Done. Key outputs:\n")
cat("  - annotation_enrichment_genecentric.tsv / .png\n")
cat("  - annotation_enrichment_cpgstructural.tsv / .png\n")
cat("  - Fig_annotation_proportion_genecentric.png / _cpgstructural.png\n")
