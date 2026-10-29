#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - coverage by genomic context and its relation to DMC/DMR
# sensitivity (Reviewer 2, MC 7)
# =============================================================================
# Part A (genome-wide, all samples): per-CpG coverage of every platform and
#   sample, stratified by CpG-structural category (CpG island > shore > shelf >
#   open sea), gene-centric category (same hierarchy as script 15) and local CpG
#   density (CpGs within +/- window bp). Reported per category: coverage relative
#   to the mean coverage of that sample and platform over all reference CpGs,
#   median coverage, and the fraction of CpGs with >= min_cov reads. Reference
#   set: all CpGs of ALL_without_EPIC.csv, i.e. the CpGs measured by the
#   sequencing-based platforms (with the merged matrices used here: by all
#   platforms in all samples), so that contexts are compared on identical CpGs;
#   CpGs without a call on a platform count as coverage 0.
# Part B (Tier 1 consensus set, identical CpGs on all four sequencing-based
#   platforms): per category and platform, median per-CpG coverage, fraction of
#   tested CpGs called as DMCs and located in DMRs, and the recall of a
#   leave-one-out reference: CpGs called as DMCs (or located in DMRs) by all
#   three other platforms. The leave-one-out reference removes the biological
#   differences between categories (e.g. few tissue differences in CpG
#   islands), so that differences in recall between categories reflect
#   platform-specific sensitivity.
# Suppl. Figures 1E,F and 12E; Suppl. Tables S4B and S11F.
#
# Usage:
#   Rscript scripts/R/18_coverage_by_context.R \
#     --all_path ALL_without_EPIC.csv --dss_dir results/dmr_dss/ \
#     --outdir results/figures/ --datadir results/coverage_context/
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(annotatr)
  library(bsseq)
  library(ggplot2)
  library(patchwork)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path", type = "character", metavar = "FILE",
              help = "ALL_without_EPIC.csv (<Method>_cov_<Sample> columns; its CpGs form the reference set, see header) [required]"),
  make_option("--dss_dir",  type = "character", default = NULL, metavar = "DIR",
              help = "--datadir of the paired run of 13_DMR_DSS_analysis.R; enables part B"),
  make_option("--outdir",   type = "character", default = "results/figures/", metavar = "DIR"),
  make_option("--datadir",  type = "character", default = "results/coverage_context/", metavar = "DIR"),
  make_option("--genome",   type = "character", default = "hg38"),
  make_option("--min_cov",  type = "integer",   default = 10L),
  make_option("--window",   type = "integer",   default = 500L,
              help = "Half-width (bp) of the window for local CpG density [default: 500]"),
  make_option("--fdr_cutoff",   type = "double", default = 0.05),
  make_option("--delta_cutoff", type = "double", default = 0.1),
  make_option("--plot_only", action = "store_true", default = FALSE,
              help = "Only redraw the figures from the tables in --datadir")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$all_path) && !opt$plot_only) stop("ERROR: --all_path is required")
dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

GENOME     <- opt$genome
VALID_CHRS <- paste0("chr", c(1:22, "X", "Y"))
PLATFORMS  <- c("ONT", "PacBio", "WGEC", "TWIST", "RRBS")
SEQ_T1     <- c("ONT", "TWIST", "WGEC", "RRBS")
col.vec    <- get_colors()

# ---- Annotation: one label per CpG and hierarchy (as in script 15) ----------
GENE_PRIORITY <- paste0(GENOME, "_genes_", c("promoters", "5UTRs", "cds", "3UTRs", "firstexons",
                                             "exons", "exonintronboundaries", "intronexonboundaries",
                                             "introns", "1to5kb", "intergenic"))
CPG_PRIORITY  <- paste0(GENOME, "_cpg_", c("islands", "shores", "shelves", "inter"))
LABELS <- c(promoters = "Promoter", `5utrs` = "5' UTR", cds = "CDS", `3utrs` = "3' UTR",
            firstexons = "First exon", exons = "Exon", exonintronboundaries = "Exon-intron boundary",
            intronexonboundaries = "Intron-exon boundary", introns = "Intron",
            `1to5kb` = "1-5 kb upstream", intergenic = "Intergenic",
            islands = "CpG island", shores = "CpG shore", shelves = "CpG shelf", inter = "Open sea")
annot_key <- function(x) tolower(sub("^.*_(genes|cpg)_", "", x))
DENS_BREAKS <- c(0, 5, 10, 20, 40, 80, Inf)
DENS_LABELS <- c("1-5", "6-10", "11-20", "21-40", "41-80", ">80")

# Highest-priority overlapping category; equivalent to the slice_min() step of
# script 15, but vectorized for genome-wide CpG sets.
label_by_priority <- function(gr, db, priority, fallback) {
  lab <- rep(NA_character_, length(gr))
  for (p in priority) {
    todo <- which(is.na(lab))
    if (length(todo) == 0) break
    hit <- overlapsAny(gr[todo], db[db$type == p], ignore.strand = TRUE)
    lab[todo[hit]] <- annot_key(p)
  }
  lab[is.na(lab)] <- fallback
  unname(LABELS[lab])
}

# Number of reference CpGs within +/- window bp (including the CpG itself).
local_density <- function(chr, pos, ref_chr, ref_pos, w) {
  out <- integer(length(pos))
  for (ch in unique(chr)) {
    i <- which(chr == ch)
    r <- sort(ref_pos[ref_chr == ch])
    out[i] <- findInterval(pos[i] + w, r) - findInterval(pos[i] - w - 1, r)
  }
  out
}

STRATA <- c(cpg_structural = "CpG-structural", gene_centric = "Gene-centric",
            cpg_density = paste0("CpG density (CpGs per ", 2 * opt$window / 1000, " kb)"))
cat_levels <- c(LABELS[annot_key(CPG_PRIORITY)], LABELS[annot_key(GENE_PRIORITY)], DENS_LABELS)

if (!opt$plot_only) {
cat("[1/5] Building annotatr databases...\n")
gene_db <- build_annotations(genome = GENOME, annotations = GENE_PRIORITY)
cpg_db  <- build_annotations(genome = GENOME, annotations = CPG_PRIORITY)

annotate_cpgs <- function(dt, ref) {
  gr <- GRanges(dt$chr, IRanges(dt$start, width = 1))
  dt[, cpg_structural := label_by_priority(gr, cpg_db,  CPG_PRIORITY,  "inter")]
  dt[, gene_centric   := label_by_priority(gr, gene_db, GENE_PRIORITY, "intergenic")]
  dens <- local_density(dt$chr, dt$start, ref$chr, ref$start, opt$window)
  dt[, cpg_density := as.character(cut(dens, DENS_BREAKS, DENS_LABELS))]
  invisible(dt)
}

# ---- 2. Part A: coverage by genomic context, genome-wide ----------------------
cat("[2/5] Reading coverage columns of", opt$all_path, "...\n")
hdr      <- names(fread(opt$all_path, nrows = 0))
cov_cols <- grep(paste0("^(", paste(PLATFORMS, collapse = "|"), ")_cov_"), hdr, value = TRUE)
cov_dt   <- fread(opt$all_path, select = c("chr", "start", cov_cols))
cov_dt   <- cov_dt[chr %in% VALID_CHRS]
setkey(cov_dt, chr, start)
cat(sprintf("  Reference set: %s CpGs (rows of --all_path)\n",
            format(nrow(cov_dt), big.mark = ",")))
ref_pos <- cov_dt[, .(chr, start)]

cat("[3/5] Annotating the reference CpGs...\n")
annotate_cpgs(cov_dt, ref_pos)

tissue_of <- function(s) fifelse(grepl("^Blood", s), "Blood",
                         fifelse(grepl("^Fibro", s), "Fibroblast",
                         fifelse(grepl("^GIAB", s), "GIAB", NA_character_)))

per_sample <- rbindlist(lapply(cov_cols, function(cc) {
  v <- cov_dt[[cc]]; v[is.na(v)] <- 0
  if (sum(v > 0) == 0) return(NULL)
  m_all <- mean(v)
  rbindlist(lapply(names(STRATA), function(s) {
    data.table(category = cov_dt[[s]], v = v)[, .(
      n_CpGs      = .N,
      rel_cov     = mean(v) / m_all,
      median_cov  = as.numeric(median(v)),
      pct_ge_min  = 100 * mean(v >= opt$min_cov),
      pct_covered = 100 * mean(v > 0)), by = category][, Stratification := STRATA[[s]]]
  }))[, `:=`(Platform = sub("_cov_.*", "", cc), Sample = sub(".*_cov_", "", cc),
             genome_mean_cov = m_all)]
}))
per_sample[, Tissue := tissue_of(Sample)]
fwrite(per_sample, file.path(opt$datadir, "coverage_by_context_per_sample.tsv"), sep = "\t")

ref_share <- rbindlist(lapply(names(STRATA), function(s)
  data.table(category = cov_dt[[s]])[, .(pct_of_reference = 100 * .N / nrow(cov_dt)), by = category][
    , Stratification := STRATA[[s]]]))
tabA <- per_sample[, .(n_samples = .N,
                       rel_cov = mean(rel_cov), rel_cov_min = min(rel_cov), rel_cov_max = max(rel_cov),
                       median_cov = mean(median_cov), pct_ge_min = mean(pct_ge_min),
                       pct_covered = mean(pct_covered), genome_mean_cov = mean(genome_mean_cov)),
                   by = .(Stratification, category, Tissue, Platform)]
tabA <- merge(tabA, ref_share, by = c("Stratification", "category"))
tabA[, category := factor(category, levels = unname(cat_levels))]
tabA[, Platform := factor(Platform, levels = PLATFORMS)]
setorder(tabA, Stratification, category, Tissue, Platform)
fwrite(tabA, file.path(opt$datadir, "TableS6B_coverage_by_context.tsv"), sep = "\t")

# ---- 4. Part B: sensitivity by genomic context on Tier 1 ---------------------
if (!is.null(opt$dss_dir)) {
  cat("[4/5] Part B: DMC/DMR sensitivity by genomic context on Tier 1...\n")
  dml   <- readRDS(file.path(opt$dss_dir, "Tier1_DML.rds"))[SEQ_T1]
  dmr   <- readRDS(file.path(opt$dss_dir, "Tier1_DMR.rds"))[SEQ_T1]
  bs_t1 <- readRDS(file.path(opt$dss_dir, "BSseq_Tier1.rds"))

  key <- function(chr, pos) paste0(chr, ":", pos)
  t1 <- as.data.table(dml[[1]])[, .(chr = as.character(chr), start = as.integer(pos))]
  t1 <- unique(t1)[chr %in% VALID_CHRS]
  annotate_cpgs(t1, ref_pos)
  t1[, id := key(chr, start)]
  gr_t1 <- GRanges(t1$chr, IRanges(t1$start, width = 1))

  for (m in SEQ_T1) {
    d   <- as.data.table(dml[[m]])
    sig <- !is.na(d$fdrs) & d$fdrs < opt$fdr_cutoff & abs(d$delta_beta) >= opt$delta_cutoff
    t1[, (paste0("dmc_", m)) := id %in% key(d$chr[sig], d$pos[sig])]
    t1[, (paste0("dmr_", m)) := overlapsAny(gr_t1, dmr[[m]], ignore.strand = TRUE)]
    # Median coverage per CpG over the blood and fibroblast samples.
    cv <- do.call(cbind, lapply(c("Blood", "Fibro"), function(tis) {
      b <- bs_t1[[paste0(m, "_", tis)]]
      x <- as.matrix(getCoverage(b, type = "Cov"))
      x[match(t1$id, key(as.character(seqnames(b)), start(b))), , drop = FALSE]
    }))
    t1[, (paste0("cov_", m)) := matrixStats::rowMedians(cv, na.rm = TRUE)]
  }

  tabB <- rbindlist(lapply(names(STRATA), function(s) rbindlist(lapply(SEQ_T1, function(m) {
    others <- setdiff(SEQ_T1, m)
    ref_dmc <- Reduce(`&`, t1[, paste0("dmc_", others), with = FALSE])
    ref_dmr <- Reduce(`&`, t1[, paste0("dmr_", others), with = FALSE])
    data.table(category = t1[[s]], dmc = t1[[paste0("dmc_", m)]], dmr = t1[[paste0("dmr_", m)]],
               cov = t1[[paste0("cov_", m)]], ref_dmc = ref_dmc, ref_dmr = ref_dmr)[, .(
      n_tested         = .N,
      median_cov       = as.numeric(median(cov, na.rm = TRUE)),
      pct_DMC          = 100 * mean(dmc),
      pct_in_DMR       = 100 * mean(dmr),
      n_ref_DMC        = sum(ref_dmc),
      recall_ref_DMC   = 100 * sum(dmc & ref_dmc) / sum(ref_dmc),
      n_ref_DMR_CpGs   = sum(ref_dmr),
      recall_ref_DMR   = 100 * sum(dmr & ref_dmr) / sum(ref_dmr)), by = category][
        , `:=`(Platform = m, Stratification = STRATA[[s]])]
  }))))
  tabB[, category := factor(category, levels = unname(cat_levels))]
  tabB[, Platform := factor(Platform, levels = SEQ_T1)]
  setorder(tabB, Stratification, category, Platform)
  fwrite(tabB, file.path(opt$datadir, "TableS11F_sensitivity_by_context.tsv"), sep = "\t")
  fwrite(t1, file.path(opt$datadir, "Tier1_CpGs_context_calls.tsv.gz"), sep = "\t")
  print(tabB[Stratification != STRATA[["gene_centric"]]])

}
} # end of computation

# ---- Figures (also with --plot_only) -----------------------------------------
fA <- file.path(opt$datadir, "TableS6B_coverage_by_context.tsv")
fB <- file.path(opt$datadir, "TableS11F_sensitivity_by_context.tsv")
fac <- function(d, plats) d[, `:=`(category = factor(category, levels = unname(cat_levels)),
                                   Platform = factor(Platform, levels = plats))]
tabA <- fac(fread(fA), PLATFORMS)
# Figure: relative coverage and fraction >= min_cov per category.
plot_A <- function(stratum, metric, ylab, log_y = FALSE) {
  d <- tabA[Stratification == stratum]
  p <- ggplot(d, aes(x = category, y = .data[[metric]], colour = Platform, shape = Tissue)) +
    geom_point(position = position_dodge(width = 0.7), size = 2.4) +
    scale_colour_manual(values = col.vec) +
    labs(x = NULL, y = ylab, title = stratum) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  if (log_y) p <- p + scale_y_log10() + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40")
  p
}
pE <- wrap_plots(lapply(unname(STRATA), plot_A, metric = "rel_cov",
                        ylab = "Coverage relative to mean over all CpGs", log_y = TRUE), nrow = 1)
pF <- wrap_plots(lapply(unname(STRATA), plot_A, metric = "pct_ge_min",
                        ylab = sprintf("CpGs with >= %dx (%%)", opt$min_cov)), nrow = 1)
figA <- (pE / pF) + plot_layout(guides = "collect") + plot_annotation(tag_levels = list(c("E", "", "", "F", "", "")))
ggsave(file.path(opt$outdir, "SuppFig1EF_coverage_by_context.png"), figA,
       width = 20, height = 11, dpi = 400)

if (file.exists(fB)) {
  tabB <- fac(fread(fB), SEQ_T1)
  plot_B <- function(stratum, metric, ylab) {
    ggplot(tabB[Stratification == stratum & n_ref_DMC >= 100],
           aes(x = category, y = .data[[metric]], colour = Platform, group = Platform)) +
      geom_point(size = 2.4) +
      (if (grepl("density", stratum)) geom_line(alpha = 0.5, show.legend = FALSE) else NULL) +
      scale_colour_manual(values = col.vec) +
      coord_cartesian(ylim = c(50, 100)) +
      labs(x = NULL, y = ylab, title = stratum) +
      theme_bw(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  }
  strata_B <- unname(STRATA[c("cpg_structural", "cpg_density")])
  figB <- wrap_plots(c(lapply(strata_B, plot_B, metric = "recall_ref_DMC",
                              ylab = "Recall of reference DMCs (%)"),
                       lapply(strata_B, plot_B, metric = "recall_ref_DMR",
                              ylab = "Recall of reference DMR CpGs (%)")), nrow = 2) +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = list(c("E", "", "", "")))
  ggsave(file.path(opt$outdir, "SuppFig12E_sensitivity_by_context.png"), figB,
         width = 14, height = 10, dpi = 400)
}

cat("[5/5] Done. Tables in", normalizePath(opt$datadir), "\n")
