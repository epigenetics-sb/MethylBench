#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - annotation enrichment relative to the tested background
# =============================================================================
# Annotates significant sites (--level CpG: DMCs; --level DMR: tested CpGs
# within DMRs) and the tested Tier 1/Tier 2 consensus background with two
# independent annotation hierarchies (gene-centric, CpG-structural) and tests
# enrichment with Fisher's exact test.
# Figure 9C, Suppl. Figures 9 and 18, Suppl. Tables S3/S4.
#
# Usage:
#   Rscript scripts/R/15_annotation_enrichment_background.R \
#     --dss_dir results/dmr_dss/ --outdir results/figures/ \
#     --datadir results/annotation_enrichment/ --level CpG|DMR
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(annotatr)
  library(GenomicRanges)
})

source("scripts/R/utils/helpers.R")

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
col.vec      <- get_colors()
METHOD_ORDER <- c("ONT", "TWIST", "WGEC", "RRBS", "EPIC")

# ---- Annotation hierarchies (independent, each mutually exclusive) ----------
GENE_PRIORITY <- c(
  paste0(GENOME, "_genes_promoters"),
  paste0(GENOME, "_genes_5UTRs"),
  paste0(GENOME, "_genes_cds"),
  paste0(GENOME, "_genes_3UTRs"),
  paste0(GENOME, "_genes_firstexons"),
  paste0(GENOME, "_genes_exons"),
  paste0(GENOME, "_genes_exonintronboundaries"),
  paste0(GENOME, "_genes_intronexonboundaries"),
  paste0(GENOME, "_genes_introns"),
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

ANNOT_LABELS <- c(
  promoters            = "Promoters",
  `5utrs`              = "5' UTRs",
  firstexons           = "First exons",
  cds                  = "CDS",
  exons                = "Exons",
  exonintronboundaries = "Exon-intron boundaries",
  intronexonboundaries = "Intron-exon boundaries",
  introns              = "Introns",
  `3utrs`              = "3' UTRs",
  `1to5kb`             = "1-5 kb upstream",
  intergenic           = "Intergenic",
  islands              = "CpG islands",
  shores               = "CpG shores",
  shelves              = "CpG shelves",
  inter                = "Open sea"
)

ANNOT_COLORS <- c(
  promoters            = "#B2182B",
  `5utrs`              = "#F28E2B",
  firstexons           = "#FFC685",
  cds                  = "#1F4E89",
  exons                = "#4E9BD6",
  exonintronboundaries = "#BFDDF2",
  intronexonboundaries = "#8C6BB1",
  introns              = "#D9C8EA",
  `3utrs`              = "#4C9A5B",
  `1to5kb`             = "#CFC8BB",
  intergenic           = "#6B6B6B",
  islands              = "#8C1C13",
  shores               = "#E36414",
  shelves              = "#F2C14E",
  inter                = "#9BB1C4"
)

annot_key <- function(x) tolower(sub("^.*_(genes|cpg)_", "", x))

cat("[1/5] Building annotatr databases (gene-centric & CpG-structural)...\n")
gene_db <- build_annotations(genome = GENOME, annotations = GENE_PRIORITY)
cpg_db  <- build_annotations(genome = GENOME, annotations = CPG_PRIORITY)

# One label per CpG and hierarchy: the highest-priority overlapping category.
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
    gr <- GRanges(seqnames = df$chr, ranges = IRanges(df$pos, df$pos))
  } else {
    seq_col <- intersect(c("seqnames", "chr"), colnames(df))[1]
    gr <- GRanges(seqnames = df[[seq_col]], ranges = IRanges(df$start, df$end))
  }
  gr[as.character(seqnames(gr)) %in% VALID_CHRS]
}

# ---- Enrichment computation (shared across tiers) ---------------------------

# Fisher's exact test per category: significant vs. non-significant tested CpGs
# (significant sites are a subset of the background). BH correction across all
# categories and platforms of one hierarchy and tier.
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

      n_nonsig_cat   <- n_bg_cat - n_sig_cat
      n_nonsig_total <- bg_total - n_sig_total
      stopifnot(n_nonsig_cat >= 0, n_nonsig_total >= 0)
      tab <- matrix(c(
        n_sig_cat,    n_sig_total    - n_sig_cat,
        n_nonsig_cat, n_nonsig_total - n_nonsig_cat
      ), nrow = 2)

      ft <- tryCatch(fisher.test(tab), error = function(e) NULL)

      pct_sig <- 100 * n_sig_cat / n_sig_total
      pct_bg  <- 100 * n_bg_cat  / bg_total
      log2fe  <- log2((pct_sig + 1e-6) / (pct_bg + 1e-6))

      out[[length(out) + 1]] <- data.table(
        Method = m, Key = annot_key(cat_i), Category = clean_label(cat_i, GENOME),
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
  combined <- bind_rows(sig_freq, bg_freq) %>% mutate(key = annot_key(annot.type))

  unknown <- setdiff(unique(combined$key), names(ANNOT_LABELS))
  if (length(unknown) > 0)
    stop("No label/colour defined for annotation(s): ", paste(unknown, collapse = ", "))

  keys <- names(ANNOT_LABELS)[names(ANNOT_LABELS) %in% combined$key]
  combined$annot_label <- factor(unname(ANNOT_LABELS[combined$key]),
                                 levels = unname(ANNOT_LABELS[keys]))
  fill_values <- setNames(unname(ANNOT_COLORS[keys]), unname(ANNOT_LABELS[keys]))

  bg_name <- "Background"
  present <- intersect(METHOD_ORDER, unique(combined$Set))
  combined$Set <- factor(combined$Set, levels = c(bg_name, present))

  ggplot(combined, aes(x = Set, y = pct, fill = annot_label)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.85,
             colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = fill_values) +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(title = title, x = NULL, y = "Proportion (%)", fill = "Annotation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 26, colour = c("black", unname(col.vec[present]))),
          plot.title = element_text(hjust = 0.5, size = 26),
          plot.title.position = "plot",
          axis.text = element_text(size = 26),
          axis.title = element_text(size = 26),
          text = element_text(size = 26),
          panel.grid.major.x = element_blank())
}

make_fe_plot <- function(enrich_df, title) {
  enrich_df <- as.data.frame(enrich_df)
  enrich_df$Method <- factor(enrich_df$Method,
                             levels = intersect(METHOD_ORDER, unique(enrich_df$Method)))

  unknown <- setdiff(unique(enrich_df$Key), names(ANNOT_LABELS))
  if (length(unknown) > 0)
    stop("No label defined for annotation(s): ", paste(unknown, collapse = ", "))

  keys <- names(ANNOT_LABELS)[names(ANNOT_LABELS) %in% enrich_df$Key]
  enrich_df$Category <- factor(unname(ANNOT_LABELS[enrich_df$Key]),
                               levels = rev(unname(ANNOT_LABELS[keys])))

  ggplot(enrich_df, aes(x = Category, y = log2FE, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = col.vec) +
    guides(fill = guide_legend(reverse = TRUE)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(log[2]~"fold enrichment (sig. / background)")) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 26),
          plot.title = element_text(hjust = 0.5, size = 26),
          plot.title.position = "plot",
          axis.text = element_text(size = 26),
          axis.title = element_text(size = 26),
          text = element_text(size = 26),
          panel.grid.major.y = element_blank())
}

# ---- Run per tier -----------------------------------------------------------

for (tier in TIERS) {

  cat(sprintf("[2/5] [%s] Annotating tested background CpG set...\n", tier))
  bg_gr <- read_consensus_as_gr(tier)

  bg_gene <- annotate_exclusive(bg_gr, gene_db, GENE_PRIORITY)
  bg_cpg  <- annotate_exclusive(bg_gr, cpg_db,  CPG_PRIORITY)

  bg_gene_freq <- bg_gene %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = "Background")
  bg_cpg_freq  <- bg_cpg  %>% count(annot.type, name = "n") %>%
    mutate(pct = n / sum(n) * 100, Set = "Background")

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
    # Region level: each DMR is represented by the tested CpGs it contains, so
    # that significant sites and background share the same unit.
    if (opt$level == "DMR") gr <- subsetByOverlaps(bg_gr, gr, ignore.strand = TRUE)
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

  lv <- opt$level
  fwrite(gene_enrich, file.path(opt$datadir, sprintf("annotation_enrichment_genecentric_%s_%s.tsv",   lv, tier)), sep = "\t")
  fwrite(cpg_enrich,  file.path(opt$datadir, sprintf("annotation_enrichment_cpgstructural_%s_%s.tsv", lv, tier)), sep = "\t")

  cat(sprintf("[5/5] [%s] Generating figures...\n", tier))

  sig_label <- if (lv == "CpG") "sig. CpGs" else "CpGs in DMRs"

  p_gene <- make_prop_plot(gene_freq, bg_gene_freq,
                           sprintf("Gene-centric annotation of %s vs. tested background: %s", sig_label, tier))
  p_cpg  <- make_prop_plot(cpg_freq, bg_cpg_freq,
                           sprintf("CpG-structural annotation of %s vs. tested background: %s", sig_label, tier))

  ggsave(file.path(opt$outdir, sprintf("Annotation_proportion_genecentric_%s_%s.png",   lv, tier)),
         p_gene, width = 14, height = 12, dpi = 300)
  ggsave(file.path(opt$outdir, sprintf("Annotation_proportion_cpgstructural_%s_%s.png", lv, tier)),
         p_cpg, width = 14, height = 12, dpi = 300)

  ggsave(file.path(opt$outdir, sprintf("Annotation_enrichment_genecentric_%s_%s.png", lv, tier)),
         make_fe_plot(gene_enrich, sprintf("Gene-centric enrichment (%s) vs. tested background: %s", sig_label, tier)),
         width = 14, height = 12, dpi = 300)
  ggsave(file.path(opt$outdir, sprintf("Annotation_enrichment_cpgstructural_%s_%s.png", lv, tier)),
         make_fe_plot(cpg_enrich, sprintf("CpG-structural enrichment (%s) vs. tested background: %s", sig_label, tier)),
         width = 14, height = 12, dpi = 300)
}

cat("\nDone.\n")
cat("  Enrichment tables : ", normalizePath(opt$datadir), "\n", sep = "")
cat("  Figures             : ", normalizePath(opt$outdir), "\n", sep = "")
