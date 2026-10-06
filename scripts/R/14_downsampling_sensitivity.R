#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Depth-Matched Downsampling Sensitivity Analysis (Tier 1)
# Supplementary Figure 16, Supplementary Table S5
# =============================================================================
# Description:
#   Addresses Reviewer 1, Major Comment 3: "Given that sequencing depth remains
#   substantially different across platforms within Tier 1 (Figure 8), please
#   clarify whether a controlled downsampling or depth-matched sensitivity
#   analysis was considered."
#
#   Starting from the Tier 1 BSseq objects and the paired results of
#   13_DMR_DSS_analysis.R (same --datadir), this script:
#
#     1. Determines a common target coverage T: the lowest of the
#        platform-wise median per-CpG coverages on Tier 1 (or --target_cov).
#     2. Downsamples every sample of every platform: at each CpG with
#        coverage > T, exactly T reads are drawn without replacement, and the
#        methylated count is drawn from a hypergeometric distribution
#        (M ~ Hypergeom(M_orig, Cov_orig - M_orig, T)). CpGs with coverage
#        <= T are left unchanged. Per-CpG coverage is thus capped at T on all
#        platforms, and the methylation proportion is preserved in
#        expectation.
#     3. Re-runs the IDENTICAL primary workflow of 13_DMR_DSS_analysis.R on
#        the downsampled data: paired Beta-Binomial model
#        (DSS::DMLfit.multiFactor, ~ subject + group, smoothing span 500 bp),
#        Wald test of the group coefficient, DMC calling at FDR < 0.05 and
#        |delta beta| >= 0.1, and DMRcate on the paired per-CpG statistics.
#     4. Repeats steps 2-3 for --n_reps independent replicates (fixed seeds).
#     5. Compares the downsampled results with the non-downsampled ("before")
#        results, which are read directly from 13's output (Tier1_DML.rds,
#        Tier1_DMR.rds), so that "before" is exactly Figures 8/9:
#          - number of DMCs and DMRs per platform,
#          - pairwise Pearson r of CpG-level delta beta (all Tier 1 CpGs),
#          - pairwise DMC Jaccard index,
#          - pairwise base-pair-weighted DMR Jaccard index.
#
#   CHANGELOG (post-review): replaces the earlier version of this script,
#   which (a) used the UNPAIRED DSS::DMLtest() + callDMR() instead of the
#   paired DMLfit.multiFactor() + DMRcate workflow of the primary analysis,
#   (b) redrew M binomially from the original methylation proportion instead
#   of subsampling the observed reads, and (c) did not report DMC/DMR counts.
#
# Input:
#   --dss_dir       --datadir of 13_DMR_DSS_analysis.R (paired run); must
#                   contain BSseq_Tier1.rds, Tier1_DML.rds, Tier1_DMR.rds
#   --outdir        Output directory for figures
#   --datadir       Output directory for tables / RDS
#   --target_cov    Fixed target coverage [default: auto, see step 1]
#   --n_reps        Number of downsampling replicates [default: 3]
#   --seed          RNG seed [default: 42]
#   --cores         Platforms processed in parallel (fork) [default: 1]
#   --plot_only     Only redraw the figure panels from the tables in --datadir
#                   (no downsampling/re-analysis; needs a previous full run)
#
# Output (--datadir):
#   - downsampling_coverage.tsv        median/mean coverage per platform, before/after, target T
#   - downsampling_counts.tsv          DMCs/DMRs per platform, before vs. after (mean, SD, per replicate)
#   - downsampling_pairwise.tsv        Pearson r (delta beta), DMC Jaccard, DMR Jaccard, before vs. after
#   - Table_S5_downsampling.tsv        combined summary used for Supplementary Table S5
#   - downsampling_results_rep<i>.rds  DML/DMC/DMR of each replicate
# Output (--outdir):
#   - SupplFig16A_coverage, SupplFig16B_DMC_counts, SupplFig16C_DMR_counts,
#     SupplFig16D_deltabeta_r, SupplFig16E_DMR_jaccard, SupplFig16F_DMC_jaccard
#     -- one panel per file, each as .png (300 dpi) and .pdf (vector, for Inkscape)
#
# Usage:
#   Rscript scripts/R/14_downsampling_sensitivity.R \
#     --dss_dir  results/dmr_dss/ \
#     --outdir   results/figures/ \
#     --datadir  results/dmr_dss/downsampling/ \
#     --n_reps 3 --cores 4
#
# Author: MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(bsseq)
  library(BiocGenerics)
  library(GenomicRanges)
  library(DMRcate)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(parallel)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--dss_dir", type = "character", metavar = "DIR",
              help = "--datadir of 13_DMR_DSS_analysis.R (paired run) [required]"),
  make_option("--outdir",  type = "character", default = "results/figures/", metavar = "DIR"),
  make_option("--datadir", type = "character", default = "results/dmr_dss/downsampling/", metavar = "DIR"),
  make_option("--target_cov", type = "double", default = NA,
              help = "Fixed target coverage; default = lowest platform median on Tier 1"),
  make_option("--n_reps", type = "integer", default = 3, help = "Downsampling replicates [default: 3]"),
  make_option("--seed", type = "integer", default = 42),
  make_option("--cores", type = "integer", default = 1, help = "Platforms in parallel [default: 1]"),
  make_option("--fdr_cutoff", type = "double", default = 0.05),
  make_option("--delta_cutoff", type = "double", default = 0.1),
  make_option("--smoothing_span", type = "integer", default = 500),
  make_option("--lambda", type = "double", default = 1000),
  make_option("--C", type = "double", default = 2),
  make_option("--min_cpgs", type = "integer", default = 3),
  make_option("--genome", type = "character", default = "hg38"),
  make_option("--plot_only", action = "store_true", default = FALSE,
              help = "Skip downsampling/re-analysis; redraw figures from the tables in --datadir")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$dss_dir)) stop("ERROR: --dss_dir is required")
need <- file.path(opt$dss_dir, if (opt$plot_only) "BSseq_Tier1.rds" else c("BSseq_Tier1.rds", "Tier1_DML.rds", "Tier1_DMR.rds"))
if (!all(file.exists(need))) stop("Missing in --dss_dir: ", paste(basename(need[!file.exists(need)]), collapse = ", "))

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

METHODS   <- c("ONT", "TWIST", "WGEC", "RRBS")
TISSUES   <- c("Blood", "Fibro")
MCOL      <- METHOD_COLORS[METHODS]            # central palette from helpers.R
COND_COL  <- c(Before = "#9DB4C0", After = "#E76F51")

# -------------------------------------------------------------------------
# Helpers (identical model and thresholds as 13_DMR_DSS_analysis.R)
# -------------------------------------------------------------------------

get_subject_id <- function(sample_names) {
  m  <- regmatches(sample_names, regexec("(?:Blood|Fibro)([0-9]+)$", sample_names))
  ok <- lengths(m) == 2
  if (!all(ok)) stop("Cannot parse subject ID from: ", paste(sample_names[!ok], collapse = ", "))
  vapply(m, `[[`, character(1), 2)
}

compute_delta_beta <- function(bs, blood_samples, fibro_samples) {
  M <- as.matrix(getCoverage(bs, type = "M")); Cov <- as.matrix(getCoverage(bs, type = "Cov"))
  beta <- M / Cov; beta[Cov == 0] <- NA
  rowMeans(beta[, blood_samples, drop = FALSE], na.rm = TRUE) -
    rowMeans(beta[, fibro_samples, drop = FALSE], na.rm = TRUE)
}

run_dss_paired <- function(bs_blood, bs_fibro) {
  blood <- sampleNames(bs_blood); fibro <- sampleNames(bs_fibro)
  bs <- BiocGenerics::combine(bs_blood, bs_fibro)
  subj <- get_subject_id(sampleNames(bs))
  design <- data.frame(
    subject = factor(subj, levels = sort(unique(subj))),
    group   = factor(c(rep("Blood", length(blood)), rep("Fibro", length(fibro))),
                     levels = c("Fibro", "Blood"))
  )
  fit  <- DMLfit.multiFactor(bs, design = design, formula = ~ subject + group,
                             smoothing = TRUE, smoothing.span = opt$smoothing_span)
  test <- DMLtest.multiFactor(fit, coef = "groupBlood")
  db <- compute_delta_beta(bs, blood, fibro)
  stopifnot(nrow(test) == length(db))
  test$delta_beta <- db
  test
}

is_dmc <- function(test) !is.na(test$fdrs) & test$fdrs < opt$fdr_cutoff &
  abs(test$delta_beta) >= opt$delta_cutoff

run_dmrcate <- function(test, label) {
  pcol <- intersect(c("pvals", "pval"), colnames(test))[1]
  gr <- GRanges(test$chr, IRanges(test$pos, width = 1),
                stat = test$stat, rawpval = test[[pcol]], diff = test$delta_beta,
                ind.fdr = test$fdrs, is.sig = is_dmc(test))
  names(gr) <- paste0(test$chr, ":", test$pos)
  res <- tryCatch({
    dmr <- dmrcate(new("CpGannotated", ranges = gr), lambda = opt$lambda, C = opt$C,
                   min.cpgs = opt$min_cpgs)
    extractRanges(dmr, genome = opt$genome)
  }, error = function(e) { cat(sprintf("    [%s] DMRcate: %s\n", label, conditionMessage(e))); NULL })
  if (is.null(res)) GRanges() else res
}

downsample_bsseq <- function(bs, target, seed) {
  Cov <- as.matrix(getCoverage(bs, type = "Cov")); M <- as.matrix(getCoverage(bs, type = "M"))
  storage.mode(Cov) <- "numeric"; storage.mode(M) <- "numeric"
  Cov[is.na(Cov)] <- 0; M[is.na(M)] <- 0
  M <- pmin(M, Cov)
  above <- which(Cov > target)
  newCov <- Cov; newM <- M
  if (length(above) > 0) {
    set.seed(seed)
    newCov[above] <- target
    # draw exactly `target` reads without replacement from the observed reads
    newM[above] <- rhyper(length(above), m = M[above], n = Cov[above] - M[above], k = target)
  }
  BSseq(chr = as.character(seqnames(bs)), pos = start(bs),
        M = unname(newM), Cov = unname(newCov), sampleNames = sampleNames(bs))
}

jaccard_ids <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_dmr <- function(g1, g2) {
  if (length(g1) == 0 || length(g2) == 0) return(NA_real_)
  g1 <- GenomicRanges::reduce(g1); g2 <- GenomicRanges::reduce(g2)
  u <- GenomicRanges::reduce(c(g1, g2)); if (length(u) == 0) return(NA_real_)
  sum(width(GenomicRanges::intersect(g1, g2))) / sum(width(u))
}

pairwise_metrics <- function(dml, dmr) {
  ids <- lapply(dml, function(t) paste0(t$chr, ":", t$pos)[is_dmc(t)])
  rbindlist(lapply(combn(METHODS, 2, simplify = FALSE), function(p) {
    a <- data.table(id = paste0(dml[[p[1]]]$chr, ":", dml[[p[1]]]$pos), d1 = dml[[p[1]]]$delta_beta)
    b <- data.table(id = paste0(dml[[p[2]]]$chr, ":", dml[[p[2]]]$pos), d2 = dml[[p[2]]]$delta_beta)
    m <- merge(a, b, by = "id")
    data.table(Pair = paste(p[1], "vs.", p[2]),
               Pearson_r   = suppressWarnings(cor(m$d1, m$d2, use = "complete.obs")),
               DMC_Jaccard = jaccard_ids(ids[[p[1]]], ids[[p[2]]]),
               DMR_Jaccard = jaccard_dmr(dmr[[p[1]]], dmr[[p[2]]]))
  }))
}

cov_values <- function(bs) { v <- as.numeric(as.matrix(getCoverage(bs, type = "Cov"))); v[!is.na(v) & v > 0] }

# -------------------------------------------------------------------------
# 1. Load Tier 1 data and "before" results from 13_DMR_DSS_analysis.R
# -------------------------------------------------------------------------
cat("[1/5] Loading Tier 1 BSseq objects and paired results of 13_DMR_DSS_analysis.R...\n")
bsseq_t1   <- readRDS(file.path(opt$dss_dir, "BSseq_Tier1.rds"))
if (!opt$plot_only) {
  dml_before <- readRDS(file.path(opt$dss_dir, "Tier1_DML.rds"))[METHODS]
  dmr_before <- readRDS(file.path(opt$dss_dir, "Tier1_DMR.rds"))[METHODS]
}
cat(sprintf("  Tier 1: %d CpGs\n", nrow(bsseq_t1[[paste0(METHODS[1], "_Blood")]])))

# -------------------------------------------------------------------------
# 2. Target coverage
# -------------------------------------------------------------------------
cat("[2/5] Determining target coverage...\n")
cov_before <- lapply(setNames(METHODS, METHODS), function(m)
  unlist(lapply(TISSUES, function(t) cov_values(bsseq_t1[[paste0(m, "_", t)]]))))
med_before <- sapply(cov_before, median)
print(round(med_before, 1))
TARGET <- if (!is.na(opt$target_cov)) opt$target_cov else floor(min(med_before))
cat(sprintf("  -> Target coverage T = %dx (lowest platform median: %s)\n",
            as.integer(TARGET), names(which.min(med_before))))

# -------------------------------------------------------------------------
# 3. Downsampling replicates: paired DSS + DMRcate per platform
# -------------------------------------------------------------------------
if (!opt$plot_only) {
cat(sprintf("[3/5] Downsampling and re-analysis (%d replicates, %d core(s))...\n", opt$n_reps, opt$cores))

run_platform <- function(m, rep_i) {
  mi <- match(m, METHODS)
  b <- downsample_bsseq(bsseq_t1[[paste0(m, "_Blood")]], TARGET, opt$seed + 1000 * rep_i + 10 * mi)
  f <- downsample_bsseq(bsseq_t1[[paste0(m, "_Fibro")]], TARGET, opt$seed + 1000 * rep_i + 10 * mi + 1)
  cat(sprintf("  [rep %d | %s] DMLfit.multiFactor (paired)...\n", rep_i, m))
  test <- run_dss_paired(b, f)
  dmr  <- run_dmrcate(test, paste0("rep", rep_i, "|", m))
  cat(sprintf("  [rep %d | %s] %d DMCs, %d DMRs\n", rep_i, m, sum(is_dmc(test)), length(dmr)))
  v <- c(cov_values(b), cov_values(f))
  list(test = test, dmr = dmr,
       cov = c(median = median(v), mean = mean(v)),
       cov_sample = if (rep_i == 1) sample(v, min(1e5, length(v))) else NULL)
}

reps <- lapply(seq_len(opt$n_reps), function(i) {
  res <- if (opt$cores > 1) mclapply(METHODS, run_platform, rep_i = i, mc.cores = opt$cores)
         else lapply(METHODS, run_platform, rep_i = i)
  names(res) <- METHODS
  saveRDS(lapply(res, function(x) list(dml = x$test, dmr = x$dmr)),
          file.path(opt$datadir, sprintf("downsampling_results_rep%d.rds", i)))
  res
})

# -------------------------------------------------------------------------
# 4. Tables
# -------------------------------------------------------------------------
cat("[4/5] Writing tables...\n")

cov_tab <- rbindlist(lapply(METHODS, function(m) data.table(
  Platform = m, Target_Cov = TARGET,
  Median_Cov_before = median(cov_before[[m]]), Mean_Cov_before = mean(cov_before[[m]]),
  Median_Cov_after  = mean(sapply(reps, function(r) r[[m]]$cov["median"])),
  Mean_Cov_after    = mean(sapply(reps, function(r) r[[m]]$cov["mean"])))))
fwrite(cov_tab, file.path(opt$datadir, "downsampling_coverage.tsv"), sep = "\t")

count_rep <- rbindlist(lapply(seq_along(reps), function(i) rbindlist(lapply(METHODS, function(m)
  data.table(Replicate = i, Platform = m, DMCs = sum(is_dmc(reps[[i]][[m]]$test)),
             DMRs = length(reps[[i]][[m]]$dmr))))))
count_tab <- merge(
  data.table(Platform = METHODS,
             DMCs_before = sapply(dml_before, function(t) sum(is_dmc(t))),
             DMRs_before = sapply(dmr_before, length)),
  count_rep[, .(DMCs_after = mean(DMCs), DMCs_after_sd = sd(DMCs),
                DMRs_after = mean(DMRs), DMRs_after_sd = sd(DMRs)), by = Platform],
  by = "Platform")
count_tab[, `:=`(DMCs_change_pct = 100 * (DMCs_after - DMCs_before) / DMCs_before,
                 DMRs_change_pct = 100 * (DMRs_after - DMRs_before) / DMRs_before)]
count_tab <- count_tab[match(METHODS, Platform)]
fwrite(count_tab, file.path(opt$datadir, "downsampling_counts.tsv"), sep = "\t")
fwrite(count_rep, file.path(opt$datadir, "downsampling_counts_per_replicate.tsv"), sep = "\t")

pw_before <- pairwise_metrics(dml_before, dmr_before)[, Condition := "Before"]
pw_after_rep <- rbindlist(lapply(seq_along(reps), function(i)
  pairwise_metrics(lapply(reps[[i]], `[[`, "test"), lapply(reps[[i]], `[[`, "dmr"))[, Replicate := i]))
pw_after <- pw_after_rep[, .(Pearson_r = mean(Pearson_r), Pearson_r_sd = sd(Pearson_r),
                             DMC_Jaccard = mean(DMC_Jaccard), DMC_Jaccard_sd = sd(DMC_Jaccard),
                             DMR_Jaccard = mean(DMR_Jaccard, na.rm = TRUE),
                             DMR_Jaccard_sd = sd(DMR_Jaccard, na.rm = TRUE)), by = Pair][, Condition := "After"]
pw_tab <- rbindlist(list(pw_before, pw_after), fill = TRUE)
fwrite(pw_tab, file.path(opt$datadir, "downsampling_pairwise.tsv"), sep = "\t")

# Combined summary for Supplementary Table S5
s5_platform <- merge(cov_tab[, .(Platform, Median_Cov_before, Median_Cov_after)], count_tab, by = "Platform")
s5_pairs <- merge(pw_before[, .(Pair, r_before = Pearson_r, DMC_J_before = DMC_Jaccard, DMR_J_before = DMR_Jaccard)],
                  pw_after[, .(Pair, r_after = Pearson_r, DMC_J_after = DMC_Jaccard, DMR_J_after = DMR_Jaccard)],
                  by = "Pair")
fwrite(s5_platform[match(METHODS, Platform)], file.path(opt$datadir, "Table_S5_downsampling_platforms.tsv"), sep = "\t")
fwrite(s5_pairs, file.path(opt$datadir, "Table_S5_downsampling_pairs.tsv"), sep = "\t")

cat("\n  Target coverage T =", TARGET, "\n")
print(count_tab, digits = 3)
print(s5_pairs, digits = 3)

} else {
  cat("[3-4/5] --plot_only: reading tables from --datadir...\n")
  count_tab <- fread(file.path(opt$datadir, "downsampling_counts.tsv"))
  pw_tab    <- fread(file.path(opt$datadir, "downsampling_pairwise.tsv"))
}

# -------------------------------------------------------------------------
# 5. Supplementary Figure 16 -- one file per panel (PNG + PDF for Inkscape)
# -------------------------------------------------------------------------
cat("[5/5] Generating Supplementary Figure 16 panels...\n")

COND_LEVELS <- c("Before", "After")
COND_LABELS <- c(Before = "Before depth matching",
                 After  = sprintf("After depth matching (n = %d)", opt$n_reps))
PAIR_LEVELS <- vapply(combn(METHODS, 2, simplify = FALSE), paste, character(1), collapse = " vs. ")

panel_theme <- theme_bw() +
  theme(text = element_text(size = 26), axis.text = element_text(size = 26),
        axis.title = element_text(size = 26), plot.title = element_text(hjust = 0.5, size = 26),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 22), panel.grid.major.x = element_blank())

save_panel <- function(p, name, width = 14, height = 12) {
  ggsave(file.path(opt$outdir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(opt$outdir, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  cat("  Saved:", name, "(.png, .pdf)\n")
}

# A: per-CpG coverage before vs. after. "After" uses the replicate-1 seeds;
#    downsampling alone is cheap, so this also works with --plot_only.
set.seed(opt$seed)
cov_after <- lapply(setNames(METHODS, METHODS), function(m) {
  mi <- match(m, METHODS)
  b <- downsample_bsseq(bsseq_t1[[paste0(m, "_Blood")]], TARGET, opt$seed + 1000 + 10 * mi)
  f <- downsample_bsseq(bsseq_t1[[paste0(m, "_Fibro")]], TARGET, opt$seed + 1000 + 10 * mi + 1)
  c(cov_values(b), cov_values(f))
})
sub_n <- function(v) sample(v, min(1e5, length(v)))
cov_df <- rbindlist(lapply(METHODS, function(m) rbind(
  data.table(Platform = m, Condition = "Before", Coverage = sub_n(cov_before[[m]])),
  data.table(Platform = m, Condition = "After",  Coverage = sub_n(cov_after[[m]])))))
cov_df[, `:=`(Platform = factor(Platform, METHODS), Condition = factor(Condition, COND_LEVELS))]

pA <- ggplot(cov_df, aes(Platform, Coverage, fill = Condition)) +
  geom_violin(scale = "width", position = position_dodge(0.85), width = 0.8, colour = "grey25") +
  geom_hline(yintercept = TARGET, linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = 0.5, y = TARGET * 1.25, hjust = 0, colour = "red", size = 7,
           label = sprintf("T = %dx", as.integer(TARGET))) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = COND_COL, labels = COND_LABELS) +
  labs(title = "Per-CpG coverage", x = NULL, y = "Coverage (log10 scale)") +
  panel_theme + theme(axis.text.x = element_text(colour = MCOL[METHODS]))
save_panel(pA, "SupplFig16A_coverage")

# B1/B2: number of DMCs and DMRs (separate files, so that they can be placed freely)
count_long <- function(before, after, after_sd) {
  d <- rbind(count_tab[, .(Platform, Condition = "Before", Value = get(before), SD = NA_real_)],
             count_tab[, .(Platform, Condition = "After",  Value = get(after),  SD = get(after_sd))])
  d[, `:=`(Platform = factor(Platform, METHODS), Condition = factor(Condition, COND_LEVELS))]
  d
}
count_plot <- function(d, ylab, title) {
  ggplot(d, aes(Platform, Value, fill = Condition)) +
    geom_col(position = position_dodge(0.85), width = 0.8, colour = "grey25") +
    geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD), position = position_dodge(0.85),
                  width = 0.25, linewidth = 0.8, na.rm = TRUE) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = COND_COL, labels = COND_LABELS) +
    labs(title = title, x = NULL, y = ylab) +
    panel_theme + theme(axis.text.x = element_text(colour = MCOL[METHODS]))
}
save_panel(count_plot(count_long("DMCs_before", "DMCs_after", "DMCs_after_sd"), "#DMCs", "Number of DMCs"),
           "SupplFig16B_DMC_counts")
save_panel(count_plot(count_long("DMRs_before", "DMRs_after", "DMRs_after_sd"), "#DMRs", "Number of DMRs"),
           "SupplFig16C_DMR_counts")

# Pairwise metrics, before -> after (one fill scale, so a single legend)
pair_plot <- function(metric, sd_col, ylab, title) {
  d <- pw_tab[, .(Pair, Condition, Value = get(metric),
                  SD = if (sd_col %in% names(pw_tab)) get(sd_col) else NA_real_)]
  d[, `:=`(Pair = factor(Pair, PAIR_LEVELS), Condition = factor(Condition, COND_LEVELS))]
  ggplot(d, aes(Pair, Value)) +
    geom_line(aes(group = Pair), colour = "grey50", linewidth = 0.8) +
    geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD), width = 0.15, na.rm = TRUE) +
    geom_point(aes(fill = Condition), shape = 21, size = 6, colour = "grey25") +
    scale_fill_manual(values = COND_COL, labels = COND_LABELS) +
    labs(title = title, x = NULL, y = ylab) +
    panel_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 22))
}
save_panel(pair_plot("Pearson_r", "Pearson_r_sd",
                     expression(paste("Pearson ", italic(r), " (", Delta*beta, ")")),
                     "CpG-level effect-size concordance"), "SupplFig16D_deltabeta_r")
save_panel(pair_plot("DMR_Jaccard", "DMR_Jaccard_sd", "Base-pair-weighted Jaccard index",
                     "DMR-level overlap"), "SupplFig16E_DMR_jaccard")
save_panel(pair_plot("DMC_Jaccard", "DMC_Jaccard_sd", "Jaccard index",
                     "DMC-level overlap"), "SupplFig16F_DMC_jaccard")

if (!opt$plot_only) writeLines(capture.output(sessionInfo()), file.path(opt$datadir, "sessionInfo.txt"))
cat("\nDone.\n")
cat("  Tables  : ", normalizePath(opt$datadir), "\n", sep = "")
cat("  Figures : ", normalizePath(opt$outdir), "\n", sep = "")
