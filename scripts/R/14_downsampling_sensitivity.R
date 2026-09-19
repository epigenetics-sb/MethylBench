#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Depth-Matched Downsampling Sensitivity Analysis (Tier 1)
# =============================================================================
# Description:
#   Addresses Reviewer 1, Major Comment: "Given that sequencing depth remains
#   substantially different across platforms within Tier 1 (Figure 8), please
#   clarify whether a controlled downsampling or depth-matched sensitivity
#   analysis was considered."
#
#   Starting from the existing Tier 1 BSseq objects (output of
#   13_DMR_DSS_analysis.R: BSseq_Tier1.rds), this script:
#
#     1. Determines a common target coverage T_target (default: the median
#        per-CpG coverage of the lowest-depth platform, i.e. ONT, unless
#        overridden via --target_cov).
#     2. For every sample/platform with mean coverage > T_target, performs
#        binomial read thinning of BOTH total coverage (Cov) and methylated
#        counts (M) down to T_target, preserving the per-CpG methylation
#        fraction in expectation. CpGs already at or below T_target are left
#        unchanged (cannot be "upsampled").
#     3. Re-runs the identical DSS::DMLtest() + callDMR() workflow used in
#        13_DMR_DSS_analysis.R on the depth-matched data.
#     4. Recomputes pairwise CpG-level Pearson correlation of delta-beta and
#        pairwise DMR Jaccard indices under depth-matched conditions.
#     5. Reports, for each platform pair, the concordance metric BEFORE vs.
#        AFTER depth matching, so that any residual platform effect that is
#        NOT attributable to coverage differences becomes directly visible.
#
#   This is a binomial-thinning depth-matching approach (as commonly used in
#   ChIP-seq/RNA-seq/bisulfite-seq benchmarking, e.g. via `rbinom`), applied
#   directly to the already-tabulated per-CpG M/Cov counts. It does not
#   require re-alignment of raw reads.
#
# Input:
#   --bsseq_tier1   Path to BSseq_Tier1.rds produced by 13_DMR_DSS_analysis.R
#   --outdir        Output directory for figures
#   --datadir       Output directory for tables / intermediate RDS
#   --target_cov    Optional fixed target coverage (default: NULL -> auto)
#   --n_reps        Number of independent thinning replicates [default: 5]
#   --seed          RNG seed [default: 42]
#
# Output:
#   - downsampling_target_coverage.tsv          chosen target depth + rationale
#   - downsampling_dml_TierD.rds                DMLtest results per replicate
#   - downsampling_dmr_TierD.rds                callDMR results per replicate
#   - downsampling_deltabeta_correlation.tsv    pairwise Pearson r, pre/post
#   - downsampling_dmr_jaccard.tsv              pairwise DMR Jaccard, pre/post
#   - Fig_S_downsampling_concordance.png        before/after summary figure
#
# Usage:
#   Rscript scripts/R/14_downsampling_sensitivity.R \
#     --bsseq_tier1  results/dmr_dss/BSseq_Tier1.rds \
#     --outdir       results/figures/ \
#     --datadir      results/dmr_dss/downsampling/
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--bsseq_tier1", type = "character", metavar = "FILE",
              help = "Path to BSseq_Tier1.rds [required]"),
  make_option("--outdir",  type = "character", default = "results/figures/",
              metavar = "DIR"),
  make_option("--datadir", type = "character", default = "results/dmr_dss/downsampling/",
              metavar = "DIR"),
  make_option("--target_cov", type = "double", default = NA,
              help = "Fixed target coverage; default = auto (min per-platform median)"),
  make_option("--n_reps", type = "integer", default = 5,
              help = "Number of independent thinning replicates [default: 5]"),
  make_option("--seed", type = "integer", default = 42),
  make_option("--smoothing_span", type = "integer", default = 500),
  make_option("--delta_cutoff", type = "double", default = 0.1),
  make_option("--p_threshold", type = "double", default = 1e-5),
  make_option("--minlen", type = "integer", default = 50),
  make_option("--min_cpgs", type = "integer", default = 3),
  make_option("--dis_merge", type = "integer", default = 50),
  make_option("--pct_sig", type = "double", default = 0.5)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$bsseq_tier1)) stop("ERROR: --bsseq_tier1 is required")
if (!file.exists(opt$bsseq_tier1)) stop(paste("File not found:", opt$bsseq_tier1))

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

set.seed(opt$seed)

METHOD_COLORS <- c(ONT = "#E69F00", TWIST = "#009E73", WGEC = "#5654E9", RRBS = "#0072B2")
TISSUES <- c("Blood", "Fibro")
METHODS <- names(METHOD_COLORS)

cat("[1/6] Loading Tier 1 BSseq objects...\n")
bsseq_t1 <- readRDS(opt$bsseq_tier1)

# -------------------------------------------------------------------------
# 1. Determine common target coverage
# -------------------------------------------------------------------------
cat("[2/6] Determining common target coverage...\n")

median_cov_per_method <- sapply(METHODS, function(m) {
  covs <- unlist(lapply(TISSUES, function(t) {
    key <- paste0(m, "_", t)
    as.numeric(getCoverage(bsseq_t1[[key]]))
  }))
  covs <- covs[covs > 0]
  median(covs)
})

cat("  Median per-CpG coverage on Tier 1 CpGs (>0x only):\n")
print(round(median_cov_per_method, 1))

TARGET_COV <- if (!is.na(opt$target_cov)) {
  opt$target_cov
} else {
  floor(min(median_cov_per_method))
}

cat(sprintf("  -> Target coverage for depth matching: %.0fx\n", TARGET_COV))

fwrite(
  data.table(
    Method       = names(median_cov_per_method),
    Median_Cov   = round(median_cov_per_method, 2),
    Target_Cov   = TARGET_COV
  ),
  file.path(opt$datadir, "downsampling_target_coverage.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 2. Binomial thinning of Cov and M to the target coverage
# -------------------------------------------------------------------------
cat("[3/6] Performing binomial depth-matching (thinning)...\n")

thin_bsseq <- function(bs, target_cov, seed_offset = 0) {

  Cov <- as.matrix(getCoverage(bs))
  M   <- as.matrix(getCoverage(bs, type = "M"))

  storage.mode(Cov) <- "numeric"
  storage.mode(M)   <- "numeric"

  new_Cov <- Cov
  new_M   <- M

  above <- which(Cov > target_cov, arr.ind = FALSE)

  if (length(above) > 0) {
    p <- target_cov / Cov[above]
    p[p > 1] <- 1

    set.seed(opt$seed + seed_offset)
    new_Cov[above] <- rbinom(length(above), size = Cov[above], prob = p)

    # thin methylated counts proportionally (hypergeometric approximation
    # via binomial with the original methylation fraction, conditioned on
    # the newly drawn coverage)
    beta_orig <- ifelse(Cov[above] > 0, M[above] / Cov[above], 0)
    new_M[above] <- rbinom(length(above), size = new_Cov[above], prob = beta_orig)
  }

  # never exceed the newly drawn coverage
  new_M <- pmin(new_M, new_Cov)

  BSseq(
    chr         = as.character(seqnames(bs)),
    pos         = start(bs),
    M           = unname(new_M),
    Cov         = unname(new_Cov),
    sampleNames = sampleNames(bs)
  )
}

run_dss_tissue <- function(bs_blood, bs_fibro, smoothing_span) {
  bs_combined <- BiocGenerics::combine(bs_blood, bs_fibro)
  DMLtest(
    bs_combined,
    group1 = sampleNames(bs_blood),
    group2 = sampleNames(bs_fibro),
    smoothing = TRUE,
    smoothing.span = smoothing_span
  )
}

call_dss_dmr <- function(dml_result) {
  callDMR(
    dml_result,
    delta       = opt$delta_cutoff,
    p.threshold = opt$p_threshold,
    minlen      = opt$minlen,
    minCG       = opt$min_cpgs,
    dis.merge   = opt$dis_merge,
    pct.sig     = opt$pct_sig
  )
}

jaccard_dmr <- function(gr1, gr2) {
  if (is.null(gr1) || is.null(gr2) || length(gr1) == 0 || length(gr2) == 0) return(NA_real_)
  gr1 <- GenomicRanges::reduce(gr1); gr2 <- GenomicRanges::reduce(gr2)
  inter <- GenomicRanges::intersect(gr1, gr2)
  union <- GenomicRanges::reduce(c(gr1, gr2))
  if (length(union) == 0) return(NA_real_)
  sum(width(inter)) / sum(width(union))
}

# -------------------------------------------------------------------------
# 3. Baseline (as-is, "before") DML/DMR per method — for reference
# -------------------------------------------------------------------------
cat("[4/6] Running baseline (non-depth-matched) DSS per method...\n")

dml_before <- list()
dmr_before <- list()
for (m in METHODS) {
  dml_before[[m]] <- run_dss_tissue(bsseq_t1[[paste0(m, "_Blood")]],
                                     bsseq_t1[[paste0(m, "_Fibro")]],
                                     opt$smoothing_span)
  dmr_before[[m]] <- call_dss_dmr(dml_before[[m]])
}

# -------------------------------------------------------------------------
# 4. Depth-matched ("after") DML/DMR, replicated n_reps times
# -------------------------------------------------------------------------
cat(sprintf("[5/6] Running depth-matched DSS per method (%d replicates)...\n", opt$n_reps))

dml_after  <- vector("list", opt$n_reps)
dmr_after  <- vector("list", opt$n_reps)

for (rep_i in seq_len(opt$n_reps)) {
  cat(sprintf("  Replicate %d/%d\n", rep_i, opt$n_reps))
  dml_after[[rep_i]] <- list()
  dmr_after[[rep_i]] <- list()
  for (m in METHODS) {
    bs_blood_thin <- thin_bsseq(bsseq_t1[[paste0(m, "_Blood")]], TARGET_COV, seed_offset = rep_i * 100)
    bs_fibro_thin <- thin_bsseq(bsseq_t1[[paste0(m, "_Fibro")]], TARGET_COV, seed_offset = rep_i * 100 + 1)
    dml_after[[rep_i]][[m]] <- run_dss_tissue(bs_blood_thin, bs_fibro_thin, opt$smoothing_span)
    dmr_after[[rep_i]][[m]] <- call_dss_dmr(dml_after[[rep_i]][[m]])
  }
}

saveRDS(dml_before, file.path(opt$datadir, "downsampling_dml_before.rds"))
saveRDS(dmr_before, file.path(opt$datadir, "downsampling_dmr_before.rds"))
saveRDS(dml_after,  file.path(opt$datadir, "downsampling_dml_after.rds"))
saveRDS(dmr_after,  file.path(opt$datadir, "downsampling_dmr_after.rds"))

# -------------------------------------------------------------------------
# 5. Pairwise concordance: delta-beta correlation + DMR Jaccard, pre/post
# -------------------------------------------------------------------------
cat("[6/6] Computing pairwise concordance before vs. after depth matching...\n")

pairs <- combn(METHODS, 2, simplify = FALSE)

get_deltabeta_cor <- function(dml_list, pairs) {
  rbindlist(lapply(pairs, function(pr) {
    d1 <- as.data.table(dml_list[[pr[1]]])
    d2 <- as.data.table(dml_list[[pr[2]]])
    setnames(d1, c("chr", "pos"), c("chr1", "pos1"))
    setnames(d2, c("chr", "pos"), c("chr2", "pos2"))
    d1[, cpg_id := paste0(chr1, ":", pos1)]
    d2[, cpg_id := paste0(chr2, ":", pos2)]
    m <- merge(d1[, .(cpg_id, diff)], d2[, .(cpg_id, diff)], by = "cpg_id",
               suffixes = c("_1", "_2"))
    r <- suppressWarnings(cor(m$diff_1, m$diff_2, use = "complete.obs"))
    data.table(Method1 = pr[1], Method2 = pr[2], n_CpGs = nrow(m), Pearson_r = r)
  }))
}

cor_before <- get_deltabeta_cor(dml_before, pairs)
cor_before[, Condition := "Before depth-matching"]

cor_after_list <- lapply(seq_len(opt$n_reps), function(i) {
  d <- get_deltabeta_cor(dml_after[[i]], pairs)
  d[, Replicate := i]
  d
})
cor_after <- rbindlist(cor_after_list)
cor_after_summary <- cor_after[, .(
  n_CpGs    = round(mean(n_CpGs)),
  Pearson_r = mean(Pearson_r),
  Pearson_r_sd = sd(Pearson_r)
), by = .(Method1, Method2)]
cor_after_summary[, Condition := "After depth-matching (mean of replicates)"]

cor_combined <- rbindlist(list(
  cor_before[, .(Method1, Method2, n_CpGs, Pearson_r, Condition)],
  cor_after_summary[, .(Method1, Method2, n_CpGs, Pearson_r, Condition)]
))

fwrite(cor_combined, file.path(opt$datadir, "downsampling_deltabeta_correlation.tsv"), sep = "\t")

jaccard_before <- rbindlist(lapply(pairs, function(pr) {
  data.table(Method1 = pr[1], Method2 = pr[2],
             Jaccard = jaccard_dmr(dmr_before[[pr[1]]], dmr_before[[pr[2]]]),
             Condition = "Before depth-matching")
}))

jaccard_after_list <- lapply(seq_len(opt$n_reps), function(i) {
  rbindlist(lapply(pairs, function(pr) {
    data.table(Method1 = pr[1], Method2 = pr[2],
               Jaccard = jaccard_dmr(dmr_after[[i]][[pr[1]]], dmr_after[[i]][[pr[2]]]),
               Replicate = i)
  }))
})
jaccard_after <- rbindlist(jaccard_after_list)
jaccard_after_summary <- jaccard_after[, .(
  Jaccard = mean(Jaccard, na.rm = TRUE)
), by = .(Method1, Method2)]
jaccard_after_summary[, Condition := "After depth-matching (mean of replicates)"]

jaccard_combined <- rbindlist(list(jaccard_before, jaccard_after_summary), fill = TRUE)
fwrite(jaccard_combined, file.path(opt$datadir, "downsampling_dmr_jaccard.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# 6. Summary figure
# -------------------------------------------------------------------------
cor_combined[, Pair := paste(Method1, "vs.", Method2)]
jaccard_combined[, Pair := paste(Method1, "vs.", Method2)]

p_cor <- ggplot(cor_combined, aes(x = Pair, y = Pearson_r, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(title = "CpG-level delta-beta concordance: before vs. after depth-matching",
       x = NULL, y = expression(paste("Pearson ", italic(r), " (", Delta*beta, ")"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

p_jac <- ggplot(jaccard_combined, aes(x = Pair, y = Jaccard, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(title = "DMR-level Jaccard concordance: before vs. after depth-matching",
       x = NULL, y = "Base-pair-weighted Jaccard index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(file.path(opt$outdir, "Fig_S_downsampling_deltabeta_concordance.png"),
       p_cor, width = 9, height = 6, dpi = 300)
ggsave(file.path(opt$outdir, "Fig_S_downsampling_dmr_jaccard.png"),
       p_jac, width = 9, height = 6, dpi = 300)

cat("\nDone. Key outputs:\n")
cat("  - downsampling_target_coverage.tsv\n")
cat("  - downsampling_deltabeta_correlation.tsv\n")
cat("  - downsampling_dmr_jaccard.tsv\n")
cat("  - Fig_S_downsampling_deltabeta_concordance.png\n")
cat("  - Fig_S_downsampling_dmr_jaccard.png\n")
