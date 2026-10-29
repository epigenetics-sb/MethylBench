#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - detection of imprinted germline DMRs (Reviewer 2, MC 9)
# =============================================================================
# Uses known imprinted germline DMRs (gDMRs) as an internal reference for
# allele-specific methylation: in somatic tissue, one parental allele is
# methylated and the other is not, so the expected bulk methylation level is
# ~0.5 and should not differ between blood and fibroblasts.
# For every platform and sample, the per-platform methylation calls (not the
# merged matrices, which contain only CpGs shared by all platforms) are
# intersected with the gDMRs to quantify
#   (1) detectability: CpGs/probes per gDMR, fraction of gDMRs with >= 3 CpGs
#       at >= min_cov (EPIC: >= 3 probes; also >= 1 probe);
#   (2) accuracy/precision: mean beta per gDMR (expected ~0.5), spread of
#       per-CpG beta values within gDMRs;
#   (3) cross-platform concordance of gDMR methylation within each sample;
#   (4) optionally, DMCs located in gDMRs in the primary analysis (negative
#       control: no blood-fibroblast difference expected).
# Suppl. Figure 4 (formerly 21), Suppl. Table S12.
#
# Usage:
#   Rscript scripts/R/17_imprinted_regions.R \
#     --samplesheet samplesheet.csv --calls_dir /path/to/calls \
#     --epic_path EPIC_blood_fibro.tsv,EPIC_GIAB.tsv \
#     --gdmr_bed imprinted_gDMRs_hg38.bed \
#     --dss_dir results/dmr_dss/ --outdir results/figures/ --datadir results/imprinting/
#
# --calls_dir layout: <calls_dir>/<METHOD>/<Sample><suffix>; the suffixes are
#   options (defaults in brackets). Each raw format is read with its own reader:
#   ONT     modkit pileup bedMethyl, --cpg --combine-strands, 0-based [_modkit_pileup.bed]
#   WGEC, TWIST, RRBS (GIAB)
#           Bismark coverage [.bismark.cov.gz]; per-strand (.cov, end == start)
#           or coverage2cytosine --merge_CpG output is detected automatically
#   RRBS (blood/fibroblast)
#           lab-specific space-separated table: chr pos coverage n_unmeth fraction,
#           1-based, per strand [.methylation.tab]; used whenever present,
#           otherwise the Bismark file of the sample is read
#   PacBio  pb-CpG-tools combined BED, 0-based [.combined.bed]
# All positions are converted to 1-based CpG positions (C of the + strand).
# --epic_path: EPIC beta matrices (tab or comma) with a 'coord' column (chr:pos)
#   or 'chr' and 'start' columns (1-based C position),
#   and one column per sample named EPIC_<Sample>, as joined in buildMergedMatrix().
# --gdmr_bed: BED (hg38, 0-based) with columns chr, start, end, name.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--samplesheet", type = "character", help = "Sample sheet (columns Sample, Sampleset) [required]"),
  make_option("--calls_dir",   type = "character", help = "Directory with per-platform methylation calls [required]"),
  make_option("--epic_path",   type = "character", default = NULL, help = "Comma-separated EPIC beta matrices (optional)"),
  make_option("--gdmr_bed",    type = "character", help = "BED file of imprinted gDMRs, hg38 [required]"),
  make_option("--dss_dir",     type = "character", default = NULL, help = "--datadir of 13_DMR_DSS_analysis.R (optional)"),
  make_option("--outdir",      type = "character", default = "results/figures/"),
  make_option("--datadir",     type = "character", default = "results/imprinting/"),
  make_option("--ont_suffix",     type = "character", default = "_modkit_pileup.bed"),
  make_option("--bismark_suffix", type = "character", default = ".bismark.cov.gz"),
  make_option("--rrbs_suffix",    type = "character", default = ".methylation.tab"),
  make_option("--pacbio_suffix",  type = "character", default = ".combined.bed"),
  make_option("--min_cov",     type = "double",  default = 10),
  make_option("--min_cpgs",    type = "integer", default = 3),
  make_option("--fdr_cutoff",  type = "double",  default = 0.05),
  make_option("--delta_cutoff",type = "double",  default = 0.1),
  make_option("--plot_only",   action = "store_true", default = FALSE,
              help = "Only redraw the figure from the tables of a previous run in --datadir")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$plot_only) for (a in c("samplesheet", "calls_dir", "gdmr_bed")) if (is.null(opt[[a]])) stop("--", a, " is required")

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

PLATFORMS <- c("ONT", "PacBio", "RRBS", "WGEC", "TWIST", "EPIC")
MCOL <- METHOD_COLORS[intersect(PLATFORMS, names(METHOD_COLORS))]

# ---- Supplementary Figure 4 (formerly 21) ---------------------------------------
make_figure <- function(summary_tab, per_gdmr, cpg_level) {
  th <- theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5))
  tissue_lab <- c(Blood = "Blood", Fibro = "Fibroblast", GIAB = "GIAB")

  pA <- summary_tab |>
    mutate(Platform = factor(Platform, levels = PLATFORMS)) |>
    ggplot(aes(Platform, pct_detectable, fill = Platform)) +
    geom_col(width = 0.7) + scale_fill_manual(values = MCOL, guide = "none") +
    facet_grid(~ Tissue, scales = "free_x", space = "free_x", labeller = as_labeller(tissue_lab)) +
    labs(title = sprintf("gDMRs with >= %d CpGs at >= %gx (EPIC: probes)", opt$min_cpgs, opt$min_cov),
         x = NULL, y = "% of gDMRs") + th +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pB <- per_gdmr[detectable == TRUE] |>
    mutate(Platform = factor(Platform, levels = PLATFORMS)) |>
    ggplot(aes(Platform, beta_mean, fill = Platform)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_violin(scale = "width", alpha = 0.7) + geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white") +
    scale_fill_manual(values = MCOL, guide = "none") + coord_cartesian(ylim = c(0, 1)) +
    facet_grid(~ Tissue, scales = "free_x", space = "free_x", labeller = as_labeller(tissue_lab)) +
    labs(title = "Mean methylation per gDMR", x = NULL, y = expression("Mean " * beta)) + th +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pC <- cpg_level |>
    mutate(Platform = factor(Platform, levels = PLATFORMS)) |>
    ggplot(aes(beta, colour = Platform)) +
    geom_density(adjust = 1.5, linewidth = 0.9) + geom_vline(xintercept = 0.5, linetype = "dashed") +
    scale_colour_manual(values = MCOL) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    facet_wrap(~ Tissue, labeller = as_labeller(tissue_lab)) +
    labs(title = "Per-CpG methylation within gDMRs", x = expression(beta), y = "Density", colour = NULL) +
    th + theme(legend.position = "bottom")

  heat <- per_gdmr[, .(beta = mean(beta_mean, na.rm = TRUE)), by = .(Tissue, Platform, gdmr)]
  heat[is.nan(beta), beta := NA_real_]
  pD <- heat |>
    mutate(Platform = factor(Platform, levels = PLATFORMS)) |>
    ggplot(aes(Platform, gdmr, fill = beta)) + geom_tile() +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0.5,
                         limits = c(0, 1), na.value = "grey85", name = expression(beta)) +
    facet_grid(~ Tissue, scales = "free_x", space = "free_x", labeller = as_labeller(tissue_lab)) +
    labs(title = "Mean methylation per gDMR (grey: not detectable)", x = NULL, y = NULL) + th +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8))

  p_all <- ((pA / pB / pC) | pD) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 28, face = "bold"))
  ggsave(file.path(opt$outdir, "SuppFig4_imprinted_gDMRs.png"),
         p_all, width = 22, height = 18, dpi = 300, bg = "white")
  invisible(NULL)
}

# --plot_only: redraw the figure from the tables of a previous run in --datadir
if (isTRUE(opt$plot_only)) {
  calls <- fread(file.path(opt$datadir, "gDMR_calls_long.tsv.gz"))
  calls[, passes := fifelse(Platform == "EPIC", !is.na(beta), !is.na(cov) & cov >= opt$min_cov)]
  make_figure(summary_tab = fread(file.path(opt$datadir, "TableS12A_gDMR_detection.tsv")),
              per_gdmr    = fread(file.path(opt$datadir, "gDMR_per_sample_platform.tsv")),
              cpg_level   = calls[passes == TRUE])
  cat("Figure redrawn from ", normalizePath(opt$datadir), "\n", sep = "")
  quit(save = "no", status = 0)
}

# ---- 1. gDMRs ------------------------------------------------------------------
gdmr <- fread(opt$gdmr_bed, header = FALSE, select = 1:4,
              col.names = c("chr", "start", "end", "gdmr"))
gdmr[, start := start + 1L]                    # BED 0-based -> 1-based
if (anyDuplicated(gdmr$gdmr)) gdmr[, gdmr := make.unique(gdmr)]
setkey(gdmr, chr, start, end)
cat(sprintf("[1/5] %d imprinted gDMRs (median width %d bp)\n", nrow(gdmr), as.integer(median(gdmr$end - gdmr$start + 1))))

in_gdmr <- function(dt) {
  dt <- dt[chr %in% gdmr$chr]
  dt[, `:=`(start = pos, end = pos)]
  ov <- foverlaps(dt, gdmr, by.x = c("chr", "start", "end"), nomatch = NULL)
  ov[, .(gdmr, chr, pos, cov, beta)]
}

# ---- 2. Per-platform calls restricted to gDMRs ---------------------------------
# One reader per raw format; output: chr, pos (1-based CpG position), cov, beta.

call_path <- function(platform, smp, suffix) file.path(opt$calls_dir, platform, paste0(smp, suffix))

# Bismark coverage files: per-strand .cov (end == start, 1-based) are merged to
# CpGs as in buildMergedMatrix(); coverage2cytosine --merge_CpG output
# (end - start == 1: 1-based; == 2: 0-based) is already merged.
read_bismark <- function(p) {
  d <- fread(p, header = FALSE, sep = "\t")
  if (ncol(d) != 6L) stop("Bismark file with ", ncol(d), " columns (expected 6): ", p)
  setnames(d, c("chr", "start", "end", "percentage", "count_meth", "count_unmeth"))
  w <- unique(d$end - d$start)
  if (identical(w, 0L) || all(w == 0)) {
    d <- mergeCpGStrands(d[, .(chr, start, count_meth, count_unmeth)])
    d[, pos := start]
  } else if (all(w == 1)) {
    d[, pos := start]
  } else if (all(w == 2)) {
    d[, pos := start + 1L]
  } else stop("Unexpected interval widths in Bismark file: ", p)
  d[, .(chr, pos, cov = count_meth + count_unmeth,
        beta = fifelse(count_meth + count_unmeth > 0, count_meth / (count_meth + count_unmeth), NA_real_))]
}

read_calls <- function(platform, smp) {
  d <- NULL; src <- NA_character_
  if (platform == "ONT") {
    p <- call_path("ONT", smp, opt$ont_suffix)
    if (file.exists(p)) {
      d <- readBedMethyl(p)[modbase == "m"]
      d <- d[, .(chr, pos = as.integer(start) + 1L, cov = Nvalid_cov, beta = fraction_mod / 100)]
      src <- p
    }
  } else if (platform == "RRBS" && file.exists(call_path("RRBS", smp, opt$rrbs_suffix))) {
    # Lab-specific RRBS table (as in 07_generate_cpg_stats.R): 1-based, per strand.
    p <- call_path("RRBS", smp, opt$rrbs_suffix)
    d <- mergeCpGStrands(readRRBSTab(p))
    d <- d[, .(chr, pos = as.integer(start), cov = count_meth + count_unmeth,
               beta = fifelse(count_meth + count_unmeth > 0, count_meth / (count_meth + count_unmeth), NA_real_))]
    src <- p
  } else if (platform %in% c("RRBS", "WGEC", "TWIST")) {
    p <- call_path(platform, smp, opt$bismark_suffix)
    if (file.exists(p)) { d <- read_bismark(p); src <- p }
  } else if (platform == "PacBio") {
    p <- call_path("PacBio", smp, opt$pacbio_suffix)
    if (file.exists(p)) {
      d <- readPacBio(p)[, .(chr, pos = as.integer(start) + 1L, cov = coverage, beta = percentage / 100)]
      src <- p
    }
  }
  if (is.null(d)) return(NULL)
  # Harmonize chromosome names to UCSC style (e.g. Ensembl "1"/"MT" -> "chr1"/"chrM").
  d[, chr := as.character(chr)]
  d[!startsWith(chr, "chr"), chr := paste0("chr", sub("^MT$", "M", chr))]
  d <- d[cov > 0][!duplicated(paste(chr, pos))]
  attr_src <- src
  out <- in_gdmr(d)
  attr(out, "src") <- attr_src
  out
}

ss <- as.data.table(readSampleSheet(opt$samplesheet))
ss[, Tissue := fifelse(Sampleset == "Fibroblast", "Fibro", Sampleset)]

cat("[2/5] Reading per-platform calls (sequencing)...\n")
calls <- rbindlist(lapply(seq_len(nrow(ss)), function(i) {
  smp <- ss$Sample[i]
  rbindlist(lapply(setdiff(PLATFORMS, "EPIC"), function(m) {
    d <- read_calls(m, smp)
    if (is.null(d)) {
      if (m != "PacBio" || ss$Sampleset[i] == "GIAB") cat(sprintf("  %-7s %-8s no file found\n", m, smp))
      return(NULL)
    }
    cat(sprintf("  %-7s %-8s %6d CpGs in gDMRs  [%s]\n", m, smp, nrow(d), basename(attr(d, "src"))))
    if (nrow(d) == 0) return(NULL)
    d[, `:=`(Platform = m, Sample = smp, Tissue = ss$Tissue[i])]
  }))
}))

if (!is.null(opt$epic_path)) {
  cat("  Reading EPIC matrices...\n")
  for (f in strsplit(opt$epic_path, ",")[[1]]) {
    e <- fread(f, header = TRUE)                       # tab- or comma-separated
    if (!"coord" %in% names(e)) {
      if (all(c("chr", "start") %in% names(e))) e[, coord := paste0(chr, ":", start)]
      else stop("EPIC matrix needs a 'coord' column or 'chr' and 'start' columns: ", f)
    }
    e[, `:=`(chr = sub(":.*", "", coord), pos = as.integer(sub(".*:", "", coord)))]
    for (smp in intersect(ss$Sample, sub("^EPIC_", "", grep("^EPIC_", names(e), value = TRUE)))) {
      tissue_smp <- ss$Tissue[ss$Sample == smp]
      d <- e[, .(chr, pos, cov = NA_real_, beta = get(paste0("EPIC_", smp)))][!is.na(beta)]
      d <- in_gdmr(d)
      d[, `:=`(Platform = "EPIC", Sample = smp, Tissue = tissue_smp)]
      calls <- rbind(calls, d)
    }
  }
}
PLAT_ORDER <- intersect(PLATFORMS, unique(calls$Platform))
fwrite(calls, file.path(opt$datadir, "gDMR_calls_long.tsv.gz"), sep = "\t")

# Reference CpGs per gDMR: union of positions observed by the whole-genome
# platforms (WGEC, ONT, PacBio) in any sample.
ref_cpgs <- unique(calls[Platform %in% c("WGEC", "ONT", "PacBio"), .(gdmr, chr, pos)])
n_ref <- ref_cpgs[, .(n_ref_cpgs = .N), by = gdmr]

# ---- 3. Per sample x platform x gDMR summaries ------------------------------------
cat("[3/5] Summarizing per gDMR...\n")
calls[, passes := fifelse(Platform == "EPIC", !is.na(beta), !is.na(cov) & cov >= opt$min_cov)]
per_gdmr <- calls[, .(n_cpgs_any = .N, n_cpgs_pass = sum(passes),
                      beta_mean = if (sum(passes) > 0) mean(beta[passes]) else NA_real_),
                  by = .(Platform, Tissue, Sample, gdmr)]
# gDMRs without any call on a platform/sample are added with 0 CpGs.
grid <- CJ(gdmr = gdmr$gdmr, Sample = unique(calls$Sample), Platform = PLAT_ORDER)
grid <- merge(grid, unique(calls[, .(Sample, Platform, Tissue)]), by = c("Sample", "Platform"))
per_gdmr <- merge(grid, per_gdmr, by = c("gdmr", "Sample", "Platform", "Tissue"), all.x = TRUE)
per_gdmr[is.na(n_cpgs_any), `:=`(n_cpgs_any = 0L, n_cpgs_pass = 0L)]
per_gdmr <- merge(per_gdmr, n_ref, by = "gdmr", all.x = TRUE)
per_gdmr[, detectable := n_cpgs_pass >= opt$min_cpgs]
per_gdmr[detectable == FALSE, beta_mean := NA_real_]
fwrite(per_gdmr, file.path(opt$datadir, "gDMR_per_sample_platform.tsv"), sep = "\t")

# Per-CpG precision at the observed coverage: binomial SE at beta = 0.5.
cpg_level <- calls[passes == TRUE]
summary_tab <- per_gdmr[, .(
  gDMRs = uniqueN(gdmr),
  pct_detectable        = 100 * mean(detectable),
  pct_with_any_call     = 100 * mean(n_cpgs_pass >= 1),
  median_cpgs_per_gDMR  = as.numeric(median(n_cpgs_pass)),
  median_pct_ref_cpgs   = median(100 * pmin(n_cpgs_pass / n_ref_cpgs, 1), na.rm = TRUE),
  median_gDMR_beta      = median(beta_mean, na.rm = TRUE),
  IQR_gDMR_beta         = IQR(beta_mean, na.rm = TRUE),
  median_abs_dev_0.5    = median(abs(beta_mean - 0.5), na.rm = TRUE)
), by = .(Platform, Tissue)]
cpg_tab <- cpg_level[, .(
  median_cov_in_gDMRs   = as.numeric(median(cov)),
  median_SE_at_0.5      = median(0.5 / sqrt(cov)),
  pct_cpg_beta_0.3_0.7  = 100 * mean(beta >= 0.3 & beta <= 0.7),
  SD_cpg_beta           = sd(beta)
), by = .(Platform, Tissue)]
summary_tab <- merge(summary_tab, cpg_tab, by = c("Platform", "Tissue"), all.x = TRUE)
setorder(summary_tab, Tissue, Platform)
fwrite(summary_tab, file.path(opt$datadir, "TableS12A_gDMR_detection.tsv"), sep = "\t")
print(summary_tab)

# ---- 4. Cross-platform concordance of gDMR methylation -----------------------------
cat("[4/5] Cross-platform concordance...\n")
wide <- dcast(per_gdmr[detectable == TRUE], Tissue + Sample + gdmr ~ Platform, value.var = "beta_mean")
plats <- intersect(PLAT_ORDER, names(wide))
conc <- rbindlist(lapply(combn(plats, 2, simplify = FALSE), function(p) {
  w <- wide[!is.na(get(p[1])) & !is.na(get(p[2]))]
  if (nrow(w) == 0) return(NULL)
  w[, .(n_gDMRs = .N,
        pearson_r = if (.N >= 3) cor(get(p[1]), get(p[2])) else NA_real_,
        MAD = mean(abs(get(p[1]) - get(p[2])))),
    by = .(Tissue, Sample)][, .(Pair = paste(p[1], "vs.", p[2]),
                                n_gDMRs = round(mean(n_gDMRs)),
                                pearson_r = mean(pearson_r, na.rm = TRUE),
                                MAD = mean(MAD)), by = Tissue]
}))
fwrite(conc, file.path(opt$datadir, "TableS12B_gDMR_concordance.tsv"), sep = "\t")
print(conc)

# Blood vs. fibroblast difference per gDMR (same individual, same platform).
bf <- per_gdmr[Tissue %in% c("Blood", "Fibro") & detectable == TRUE]
bf[, subject := sub("^(Blood|Fibro)", "", Sample)]
tis <- dcast(bf, Platform + gdmr + subject ~ Tissue, value.var = "beta_mean")
if (!all(c("Blood", "Fibro") %in% names(tis))) tis[, setdiff(c("Blood", "Fibro"), names(tis)) := NA_real_]
tis_tab <- tis[!is.na(Blood) & !is.na(Fibro),
               .(n_pairs = .N, mean_delta = mean(Blood - Fibro),
                 pct_abs_delta_ge_0.1 = 100 * mean(abs(Blood - Fibro) >= 0.1)), by = Platform]
fwrite(tis_tab, file.path(opt$datadir, "TableS12C_gDMR_blood_vs_fibro.tsv"), sep = "\t")
print(tis_tab)

# ---- 5. DMCs in gDMRs (primary analysis) -------------------------------------------
if (!is.null(opt$dss_dir)) {
  cat("[5/5] DMCs located in gDMRs (primary analysis)...\n")
  gdmr_gr <- GRanges(gdmr$chr, IRanges(gdmr$start, gdmr$end))
  load_dml <- function(tier) {
    f <- file.path(opt$dss_dir, paste0(tier, "_DML", if (tier == "Tier2") "_withEPIC" else "", ".rds"))
    if (!file.exists(f)) f <- file.path(opt$dss_dir, paste0(tier, "_DML.rds"))
    readRDS(f)
  }
  dmc_tab <- rbindlist(lapply(c("Tier1", "Tier2"), function(tier) {
    dl <- load_dml(tier)
    rbindlist(lapply(names(dl), function(m) {
      t <- dl[[m]]
      g <- GRanges(t$chr, IRanges(as.integer(t$pos), width = 1))
      inside <- overlapsAny(g, gdmr_gr)
      sig <- !is.na(t$fdrs) & t$fdrs < opt$fdr_cutoff & abs(t$delta_beta) >= opt$delta_cutoff
      data.table(Tier = tier, Platform = m, tested_in_gDMRs = sum(inside),
                 gDMRs_tested = sum(overlapsAny(gdmr_gr, g)),
                 DMCs_in_gDMRs = sum(inside & sig),
                 pct_DMC_in_gDMRs = 100 * mean(sig[inside]),
                 pct_DMC_genomewide = 100 * mean(sig))
    }))
  }))
  fwrite(dmc_tab, file.path(opt$datadir, "TableS12D_DMCs_in_gDMRs.tsv"), sep = "\t")
  print(dmc_tab)
}

make_figure(summary_tab, per_gdmr, cpg_level)

writeLines(capture.output(sessionInfo()), file.path(opt$datadir, "sessionInfo.txt"))
cat("\nDone. Tables: ", normalizePath(opt$datadir), "\n", sep = "")
