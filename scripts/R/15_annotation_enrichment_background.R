#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Background-Corrected Annotation Enrichment (Gene-centric vs.
# CpG-structural, assessed separately), Tier 1 and Tier 2
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
#   Replaces the single, mixed priority hierarchy used by the now-retired
#   13_annotation.R (ANNOT_PRIORITY, which ranked gene features above CpG-
#   structural features) with TWO independent, internally mutually-exclusive
#   hierarchies:
#
#     Gene-centric  (each CpG/DMR assigned to exactly one class):
#       Promoter > 5'UTR > CDS/Exon > Intron > 3'UTR > 1-5kb upstream >
#       Intergenic
#
#     CpG-structural (each CpG/DMR assigned to exactly one class):
#       CpG Island > CpG Shore > CpG Shelf > Open Sea
#
#   Because these two hierarchies do not compete with one another, a
#   CpG/DMR can (and does) carry one gene-centric AND one CpG-structural
#   label simultaneously; no CGI signal is discarded by gene-centric
#   prioritization and vice versa.
#
#   Both the significant-site annotation profile and the annotation profile
#   of the FULL TESTED BACKGROUND (the Tier 1/Tier 2 consensus CpG set that
#   entered testing, i.e. the set at risk of being called, not the full
#   genome) are computed, for both tiers. Enrichment is reported as: the
#   significant-set proportion in each category, the background proportion,
#   log2 fold enrichment, and a Fisher's exact test / odds ratio with 95% CI
#   per category.
#
# CHANGELOG (post-review, renumbered from 17_annotation_enrichment_
# background.R):
#   FIXED a scope mismatch that was not part of the reviewer's original
#   comment but was found while addressing it: this script previously took
#   its significant sites from --sig_bed_dir, i.e. the BED files written by
#   13_annotation.R from the EXPLORATORY-stage limma/Wilcoxon results
#   (Section 2.5.1, all 146,704 overlapping CpGs, all five platforms), while
#   its background came from --background_tsv, i.e. 14_DMR_DSS_analysis.R's
#   Tier 1/Tier 2 CONSENSUS set (Section 2.5.2/2.5.3, coverage-matched,
#   sequencing-only or +EPIC). Comparing an exploratory-stage significant
#   set against a primary-stage background mixes two different CpG
#   universes and statistical models. Since Reviewer 1's comment explicitly
#   ties this to Figure 9C (now Tier 1/Tier 2-based, in
#   13_DMR_DSS_analysis.R), this script now reads BOTH the significant
#   sites (Tier{1,2}_DML_significant_<method>.tsv for --level CpG,
#   Tier{1,2}_DMR_<method>.tsv for --level DMR) and the background
#   (Tier{1,2}_consensus_CpGs.tsv) from the SAME source
#   (13_DMR_DSS_analysis.R's --datadir), consistently scoped to Tier 1/2,
#   and runs both tiers automatically. 13_annotation.R and its BED-file
#   output are no longer used anywhere and have been removed; RENUMBERED
#   17 -> 15 accordingly.
#
# Input:
#   --dss_dir   Path to 13_DMR_DSS_analysis.R's --datadir (contains
#               Tier{1,2}_consensus_CpGs.tsv, Tier{1,2}_DML_significant_
#               <method>.tsv, Tier{1,2}_DMR_<method>.tsv) [required]
#   --outdir    Output directory for figures
#   --datadir   Output directory for tables
#   --genome    Genome build [default: hg38]
#   --level     "CpG" (per-CpG significant DMCs) or "DMR" (called DMRs)
#               [default: CpG]
#
# Output (per tier, Tier1 and Tier2):
#   - annotation_enrichment_genecentric_<Tier>.tsv
#   - annotation_enrichment_cpgstructural_<Tier>.tsv
#   - Fig_annotation_proportion_genecentric_<Tier>.png
#   - Fig_annotation_proportion_cpgstructural_<Tier>.png
#   - Fig_annotation_enrichment_genecentric_<Tier>.png
#   - Fig_annotation_enrichment_cpgstructural_<Tier>.png
#
# Usage:
#   Rscript scripts/R/15_annotation_enrichment_background.R \
#     --dss_dir  results/dmr_dss/ \
#     --outdir   results/figures/ \
#     --datadir  results/annotation_enrichment/ \
#     --level    CpG
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
  make_option("--dss_dir", type = "character", metavar = "DIR",
              help = "Path to 13_DMR_DSS_analysis.R's --datadir [required]"),
  make_option("--outdir",  type = "character", default = "results/figures/"),
  make_option("--datadir", type = "character", default = "results/annotation_enrichment/"),
  make_option("--genome",  type = "character", default = "hg38"),
  make_option("--level",   type = "character", default = "CpG",
              help = "'CpG' (significant DMCs) or 'DMR' (called DMRs) [default: CpG]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$dss_dir)) stop("ERROR: --dss_dir is required")
if (!dir.exists(opt$dss_dir)) stop(paste("Directory not found:", opt$dss_dir))
if (!opt$level %in% c("CpG", "DMR")) stop("--level must be 'CpG' or 'DMR'")

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

GENOME <- opt$genome
VALID_CHRS <- paste0("chr", c(1:22, "X", "Y"))
TIERS <- c("Tier1", "Tier2")

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

cat("[1/5] Building annotatr databases (gene-centric & CpG-structural)...\n")
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

# --- Readers: consensus background (CpG-level) and per-method significant
#     sites, at either CpG (DML) or region (DMR) granularity -----------------

read_consensus_as_gr <- function(tier) {
  path <- file.path(opt$dss_dir, paste0(tier, "_consensus_CpGs.tsv"))
  if (!file.exists(path)) stop("Background file not found: ", path)
  bg <- fread(path)
  parts <- str_match(bg$cpg_id, "^(chr[^:]+):([0-9]+)$")
  gr <- GRanges(seqnames = parts[, 2], ranges = IRanges(as.integer(parts[, 3]), as.integer(parts[, 3])))
  gr[as.character(seqnames(gr)) %in% VALID_CHRS]
}

discover_methods <- function(tier, level) {
  pattern <- if (level == "CpG") {
    paste0("^", tier, "_DML_significant_(.+)\\.tsv$")
  } else {
    paste0("^", tier, "_DMR_(.+)\\.tsv$")
  }
  files <- list.files(opt$dss_dir, pattern = pattern, full.names = TRUE)
  methods <- sub(pattern, "\\1", basename(files))
  setNames(files, methods)
}

read_sig_as_gr <- function(path, level) {
  df <- fread(path)
  if (nrow(df) == 0) return(GRanges())
  if (level == "CpG") {
    # Tier{1,2}_DML_significant_<method>.tsv: chr, pos, diff, fdr
    gr <- GRanges(seqnames = df$chr, ranges = IRanges(df$pos, df$pos))
  } else {
    # Tier{1,2}_DMR_<method>.tsv: seqnames, start, end, ... (as.data.frame(GRanges))
    seq_col <- intersect(c("seqnames", "chr"), colnames(df))[1]
    gr <- GRanges(seqnames = df[[seq_col]], ranges = IRanges(df$start, df$end))
  }
  gr[as.character(seqnames(gr)) %in% VALID_CHRS]
}

# -------------------------------------------------------------------------
# Enrichment computation (shared across tiers)
# -------------------------------------------------------------------------

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

make_fe_plot <- function(enrich_df, title) {
  ggplot(enrich_df, aes(x = Category, y = log2FE, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(log[2]~"fold enrichment (sig. / background)")) +
    theme_bw()
}

# -------------------------------------------------------------------------
# Run per tier
# -------------------------------------------------------------------------

for (tier in TIERS) {

  cat(sprintf("[2/5] [%s] Annotating tested background CpG set...\n", tier))
  bg_gr <- read_consensus_as_gr(tier)

  bg_gene <- annotate_exclusive(bg_gr, gene_db, GENE_PRIORITY)
  bg_cpg  <- annotate_exclusive(bg_gr, cpg_db,  CPG_PRIORITY)

  bg_gene_freq <- bg_gene %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = "Background (tested)")
  bg_cpg_freq  <- bg_cpg  %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = "Background (tested)")

  cat(sprintf("  [%s] Background CpGs annotated: %d (gene-centric), %d (CpG-structural)\n",
              tier, nrow(bg_gene), nrow(bg_cpg)))

  cat(sprintf("[3/5] [%s] Annotating per-method significant %ss...\n", tier, opt$level))

  method_files <- discover_methods(tier, opt$level)
  if (length(method_files) == 0) {
    cat(sprintf("  [%s] No %s files found for this tier -- skipping.\n", tier, opt$level))
    next
  }

  gene_list <- list()
  cpg_list  <- list()

  for (method in names(method_files)) {
    gr <- read_sig_as_gr(method_files[[method]], opt$level)
    if (length(gr) == 0) next

    gene_ann <- annotate_exclusive(gr, gene_db, GENE_PRIORITY)
    cpg_ann  <- annotate_exclusive(gr, cpg_db,  CPG_PRIORITY)

    gene_list[[method]] <- gene_ann %>% count(annot.type, name = "n") %>%
      mutate(pct = n / sum(n) * 100, Set = method)
    cpg_list[[method]]  <- cpg_ann %>% count(annot.type, name = "n") %>%
      mutate(pct = n / sum(n) * 100, Set = method)
  }

  gene_freq <- bind_rows(gene_list)
  cpg_freq  <- bind_rows(cpg_list)

  if (nrow(gene_freq) == 0) {
    cat(sprintf("  [%s] No significant sites survived annotation -- skipping enrichment.\n", tier))
    next
  }

  cat(sprintf("[4/5] [%s] Computing enrichment vs. background (log2FE + Fisher's exact test)...\n", tier))

  gene_totals <- gene_freq %>% group_by(Set) %>% summarise(n = sum(n)) %>% deframe()
  cpg_totals  <- cpg_freq  %>% group_by(Set) %>% summarise(n = sum(n)) %>% deframe()

  gene_enrich <- compute_enrichment(gene_freq, bg_gene_freq, gene_totals, sum(bg_gene_freq$n))
  cpg_enrich  <- compute_enrichment(cpg_freq,  bg_cpg_freq,  cpg_totals,  sum(bg_cpg_freq$n))

  fwrite(gene_enrich, file.path(opt$datadir, paste0("annotation_enrichment_genecentric_", tier, ".tsv")), sep = "\t")
  fwrite(cpg_enrich,  file.path(opt$datadir, paste0("annotation_enrichment_cpgstructural_", tier, ".tsv")), sep = "\t")

  cat(sprintf("[5/5] [%s] Generating figures...\n", tier))

  p_gene <- make_prop_plot(gene_freq, bg_gene_freq,
    sprintf("Gene-centric annotation of significant %ss vs. tested background – %s", opt$level, tier))
  p_cpg  <- make_prop_plot(cpg_freq, bg_cpg_freq,
    sprintf("CpG-structural annotation of significant %ss vs. tested background – %s", opt$level, tier))

  ggsave(file.path(opt$outdir, paste0("Fig_annotation_proportion_genecentric_", tier, ".png")),
         p_gene, width = 10, height = 7, dpi = 300)
  ggsave(file.path(opt$outdir, paste0("Fig_annotation_proportion_cpgstructural_", tier, ".png")),
         p_cpg, width = 10, height = 7, dpi = 300)

  ggsave(file.path(opt$outdir, paste0("Fig_annotation_enrichment_genecentric_", tier, ".png")),
         make_fe_plot(gene_enrich, paste("Gene-centric enrichment vs. tested background –", tier)),
         width = 9, height = 6, dpi = 300)
  ggsave(file.path(opt$outdir, paste0("Fig_annotation_enrichment_cpgstructural_", tier, ".png")),
         make_fe_plot(cpg_enrich, paste("CpG-structural enrichment vs. tested background –", tier)),
         width = 8, height = 5, dpi = 300)
}

cat("\nDone. Key outputs (per tier):\n")
cat("  - annotation_enrichment_genecentric_<Tier>.tsv / cpgstructural_<Tier>.tsv\n")
cat("  - Fig_annotation_proportion_genecentric_<Tier>.png / _cpgstructural_<Tier>.png\n")
cat("  - Fig_annotation_enrichment_genecentric_<Tier>.png / _cpgstructural_<Tier>.png\n")
