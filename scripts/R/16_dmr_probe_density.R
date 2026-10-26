#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - probe density and DMRcate region calling (Reviewer 2, MC 5)
# =============================================================================
# Quantifies how the spacing of the tested CpGs interacts with the DMRcate
# parameters (lambda, C, min.cpgs) used in 13_DMR_DSS_analysis.R:
#   (1) inter-CpG spacing and "DMR-eligible" CpG clusters (consecutive tested
#       CpGs separated by gaps < lambda, with >= min.cpgs members) on the Tier 1
#       set, the Tier 2 (EPIC-restricted) set and, optionally, all QC-passing
#       EPIC probes;
#   (2) per platform and tier: fraction of DMCs that lie in a DMR-eligible
#       cluster (upper bound for DMR membership) and fraction actually contained
#       in a DMR;
#   (3) density-only loss: Tier 1 DMRs of each sequencing platform that remain
#       detectable when the identical data are restricted to EPIC positions
#       (Tier 2 DMRs of the same platform);
#   (4) DMRcate parameter sensitivity on Tier 2 (all five platforms, identical
#       CpG positions): DMR counts, CpGs per DMR, direction and pairwise
#       base-pair- and CpG-weighted Jaccard indices.
# Suppl. Figure 20, Suppl. Table S11.
#
# Usage:
#   Rscript scripts/R/16_dmr_probe_density.R --dss_dir results/dmr_dss/ \
#     --outdir results/figures/ --datadir results/dmr_density/ \
#     [--all_path ALL.csv]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(DMRcate)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--dss_dir", type = "character", metavar = "DIR",
              help = "--datadir of the paired run of 13_DMR_DSS_analysis.R [required]"),
  make_option("--all_path", type = "character", default = NULL, metavar = "FILE",
              help = "Optional ALL.csv; adds the spacing of all QC-passing EPIC probes (blood/fibroblast cohort)"),
  make_option("--outdir",  type = "character", default = "results/figures/", metavar = "DIR"),
  make_option("--datadir", type = "character", default = "results/dmr_density/", metavar = "DIR"),
  make_option("--fdr_cutoff",   type = "double",  default = 0.05),
  make_option("--delta_cutoff", type = "double",  default = 0.1),
  make_option("--lambda",       type = "double",  default = 1000),
  make_option("--C",            type = "double",  default = 2),
  make_option("--min_cpgs",     type = "integer", default = 3),
  make_option("--min_samples",  type = "integer", default = 4),
  make_option("--genome",       type = "character", default = "hg38")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$dss_dir)) stop("ERROR: --dss_dir is required")
need <- file.path(opt$dss_dir, c("Tier1_consensus_CpGs.tsv", "Tier2_consensus_CpGs.tsv",
                                 "Tier1_DML.rds", "Tier1_DMR.rds", "Tier2_DML.rds", "Tier2_DMR.rds",
                                 "Tier2_DML_EPIC.tsv"))
if (!all(file.exists(need))) stop("Missing in --dss_dir: ", paste(basename(need[!file.exists(need)]), collapse = ", "))

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

SEQ_METHODS <- c("ONT", "TWIST", "WGEC", "RRBS")
ALL_METHODS <- c(SEQ_METHODS, "EPIC")
MCOL        <- METHOD_COLORS[ALL_METHODS]
LAMBDA      <- opt$lambda
MIN_CPGS    <- opt$min_cpgs
SIGMA       <- opt$lambda / opt$C

# DMRcate settings evaluated on Tier 2. "default" is the setting of the primary
# analysis; lambda = 2000 with C = 4 keeps the kernel (sigma = 500 bp) but
# allows larger gaps; lambda = 2000 with C = 2 doubles the kernel width;
# "consec" treats CpGs as equally spaced (bandwidth in CpGs, density-independent).
PARAM_GRID <- data.frame(
  setting      = c("default", "min.cpgs = 2", "lambda 2000, C 4", "lambda 2000, C 2", "consecutive (10 CpGs)"),
  lambda       = c(LAMBDA, LAMBDA, 2000, 2000, LAMBDA),
  C            = c(opt$C,  opt$C,  4,    2,    opt$C),
  min_cpgs     = c(MIN_CPGS, 2L, MIN_CPGS, MIN_CPGS, MIN_CPGS),
  consec       = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  conseclambda = c(10, 10, 10, 10, 10),
  stringsAsFactors = FALSE
)

# ---- Helpers ----------------------------------------------------------------

ids_to_gr <- function(ids) {
  GRanges(sub(":.*", "", ids), IRanges(as.integer(sub(".*:", "", ids)), width = 1))
}

is_dmc <- function(t) !is.na(t$fdrs) & t$fdrs < opt$fdr_cutoff & abs(t$delta_beta) >= opt$delta_cutoff

# Identical CpGannotated construction as in 13_DMR_DSS_analysis.R.
to_cpgannotated <- function(t) {
  pcol <- intersect(c("pvals", "pval", "P.Value"), colnames(t))[1]
  gr <- GRanges(t$chr, IRanges(as.integer(t$pos), width = 1),
                stat = t$stat, rawpval = t[[pcol]], diff = t$delta_beta,
                ind.fdr = t$fdrs, is.sig = is_dmc(t))
  names(gr) <- paste0(t$chr, ":", t$pos)
  new("CpGannotated", ranges = sort(gr))
}

run_dmrcate_setting <- function(t, p, label) {
  res <- tryCatch({
    d <- dmrcate(to_cpgannotated(t), lambda = p$lambda, C = p$C, min.cpgs = p$min_cpgs,
                 consec = p$consec, conseclambda = p$conseclambda)
    extractRanges(d, genome = opt$genome)
  }, error = function(e) { cat(sprintf("    [%s] DMRcate: %s\n", label, conditionMessage(e))); NULL })
  if (is.null(res)) GRanges() else res
}

# Spacing diagnostics of a CpG set. A CpG is "DMR-eligible" if it belongs to a
# run of consecutive tested CpGs whose gaps are all < lambda and which contains
# >= min.cpgs CpGs; only such CpGs can be members of a DMRcate DMR.
spacing_stats <- function(gr, set_label, lambda = LAMBDA, min_cpgs = MIN_CPGS) {
  gr <- sort(gr)
  dt <- data.table(chr = as.character(seqnames(gr)), pos = start(gr))
  dt[, `:=`(gap_prev = c(NA_integer_, diff(pos)), gap_next = c(diff(pos), NA_integer_)), by = chr]
  dt[, nn_dist := pmin(gap_prev, gap_next, na.rm = TRUE)]
  dt[, cluster := cumsum(is.na(gap_prev) | gap_prev >= lambda), by = chr]
  dt[, cluster_size := .N, by = .(chr, cluster)]
  # Tested CpGs within +-2 sigma (the effective kernel support) of each CpG.
  dt[, n_within_2sigma := {
    lo <- findInterval(pos - 2 * SIGMA, pos, left.open = TRUE)
    hi <- findInterval(pos + 2 * SIGMA, pos)
    hi - lo - 1L
  }, by = chr]
  dt[, eligible := cluster_size >= min_cpgs]
  list(
    per_cpg = dt,
    summary = data.table(
      Set = set_label, n_CpGs = nrow(dt),
      median_nn_dist_bp       = median(dt$nn_dist, na.rm = TRUE),
      pct_nn_le_sigma         = 100 * mean(dt$nn_dist <= SIGMA, na.rm = TRUE),
      pct_nn_ge_lambda        = 100 * mean(dt$nn_dist >= lambda, na.rm = TRUE),
      median_n_within_2sigma  = median(dt$n_within_2sigma),
      pct_isolated_2sigma     = 100 * mean(dt$n_within_2sigma == 0),
      pct_DMR_eligible        = 100 * mean(dt$eligible),
      n_eligible_clusters     = dt[eligible == TRUE, uniqueN(paste(chr, cluster))]
    )
  )
}

# Base-pair- and CpG-weighted Jaccard index of two DMR sets.
jaccard_bp <- function(g1, g2) {
  if (length(g1) == 0 || length(g2) == 0) return(NA_real_)
  g1 <- reduce(g1); g2 <- reduce(g2); u <- reduce(c(g1, g2))
  sum(width(GenomicRanges::intersect(g1, g2))) / sum(width(u))
}
jaccard_cpg <- function(g1, g2, cpgs) {
  if (length(g1) == 0 || length(g2) == 0) return(NA_real_)
  a <- overlapsAny(cpgs, g1); b <- overlapsAny(cpgs, g2)
  if (!any(a | b)) return(NA_real_)
  sum(a & b) / sum(a | b)
}

load_dmr_tsv <- function(path) {
  if (!file.exists(path)) return(GRanges())
  makeGRangesFromDataFrame(fread(path), keep.extra.columns = TRUE)
}

# ---- 1. Load results of 13_DMR_DSS_analysis.R --------------------------------
cat("[1/5] Loading consensus sets and paired DML/DMR results...\n")
t1_ids <- fread(file.path(opt$dss_dir, "Tier1_consensus_CpGs.tsv"))$cpg_id
t2_ids <- fread(file.path(opt$dss_dir, "Tier2_consensus_CpGs.tsv"))$cpg_id
t1_gr  <- sort(ids_to_gr(t1_ids)); t2_gr <- sort(ids_to_gr(t2_ids))

pick_rds <- function(stem) {
  f <- file.path(opt$dss_dir, paste0(stem, "_withEPIC.rds"))
  if (file.exists(f)) f else file.path(opt$dss_dir, paste0(stem, ".rds"))
}
dml <- list(Tier1 = readRDS(file.path(opt$dss_dir, "Tier1_DML.rds"))[SEQ_METHODS],
            Tier2 = readRDS(pick_rds("Tier2_DML")))
dmr <- list(Tier1 = readRDS(file.path(opt$dss_dir, "Tier1_DMR.rds"))[SEQ_METHODS],
            Tier2 = readRDS(pick_rds("Tier2_DMR")))

# Tier2_DML.rds / Tier2_DMR.rds contain EPIC only if written after step 5 of
# 13_DMR_DSS_analysis.R; otherwise EPIC is read from the exported tables.
if (!"EPIC" %in% names(dml$Tier2)) {
  e <- as.data.frame(fread(file.path(opt$dss_dir, "Tier2_DML_EPIC.tsv")))
  dml$Tier2[["EPIC"]] <- e
}
if (!"EPIC" %in% names(dmr$Tier2)) {
  dmr$Tier2[["EPIC"]] <- load_dmr_tsv(file.path(opt$dss_dir, "Tier2_DMR_EPIC.tsv"))
}
dml$Tier2 <- dml$Tier2[ALL_METHODS]; dmr$Tier2 <- dmr$Tier2[ALL_METHODS]

# ---- 2. Spacing of the tested CpG sets ---------------------------------------
cat("[2/5] Spacing diagnostics...\n")
sp <- list(Tier1 = spacing_stats(t1_gr, "Tier 1 (sequencing consensus)"),
           Tier2 = spacing_stats(t2_gr, "Tier 2 (EPIC-restricted consensus)"))

if (!is.null(opt$all_path)) {
  # ALL.csv has no cpg_id column; it is built from chr/start as in 13_DMR_DSS_analysis.R.
  hdr  <- names(fread(opt$all_path, nrows = 0, sep = ","))
  ecol <- grep("^EPIC_(Blood|Fibro)[0-9]+$", hdr, value = TRUE)
  if (length(ecol) == 0) stop("No EPIC_Blood*/EPIC_Fibro* columns in --all_path")
  epic <- fread(opt$all_path, sep = ",", na.strings = "NA", select = c("chr", "start", ecol))
  epic[, cpg_id := paste0(chr, ":", start)]
  keep <- Reduce(`&`, lapply(c("Blood", "Fibro"), function(ti) {
    cols <- grep(paste0("^EPIC_", ti), ecol, value = TRUE)
    rowSums(!is.na(as.matrix(epic[, ..cols]))) >= opt$min_samples
  }))
  sp$EPIC_all <- spacing_stats(ids_to_gr(epic$cpg_id[keep]), "EPIC (all QC-passing probes)")
  rm(epic); invisible(gc())
}
spacing_summary <- rbindlist(lapply(sp, `[[`, "summary"))
fwrite(spacing_summary, file.path(opt$datadir, "TableS11A_spacing.tsv"), sep = "\t")
print(spacing_summary)

# ---- 3. DMCs in DMR-eligible clusters and in DMRs -----------------------------
cat("[3/5] DMC eligibility and DMR membership...\n")
elig_gr <- lapply(sp[c("Tier1", "Tier2")], function(s) {
  d <- s$per_cpg[eligible == TRUE]; GRanges(d$chr, IRanges(d$pos, width = 1))
})
dmc_tab <- rbindlist(lapply(c("Tier1", "Tier2"), function(tier) {
  rbindlist(lapply(names(dml[[tier]]), function(m) {
    t <- dml[[tier]][[m]]; s <- is_dmc(t)
    g <- GRanges(t$chr[s], IRanges(as.integer(t$pos[s]), width = 1))
    data.table(Tier = tier, Method = m, DMCs = length(g),
               pct_DMC_eligible = 100 * mean(overlapsAny(g, elig_gr[[tier]])),
               pct_DMC_in_DMR   = 100 * mean(overlapsAny(g, dmr[[tier]][[m]])),
               DMRs = length(dmr[[tier]][[m]]),
               median_CpGs_per_DMR = median(dmr[[tier]][[m]]$no.cpgs),
               pct_DMR_at_min_cpgs = 100 * mean(dmr[[tier]][[m]]$no.cpgs <= MIN_CPGS))
  }))
}))
fwrite(dmc_tab, file.path(opt$datadir, "TableS11B_DMC_eligibility.tsv"), sep = "\t")
print(dmc_tab)

# ---- 4. Density-only loss (same platform, same data, Tier 1 -> Tier 2) --------
cat("[4/5] Density-only loss of Tier 1 DMRs...\n")
loss_tab <- rbindlist(lapply(SEQ_METHODS, function(m) {
  d1 <- dmr$Tier1[[m]]; d2 <- dmr$Tier2[[m]]
  n2 <- countOverlaps(d1, t2_gr)
  rec <- n2 >= MIN_CPGS
  data.table(Method = m, Tier1_DMRs = length(d1),
             pct_0_Tier2_CpGs        = 100 * mean(n2 == 0),
             pct_1_to_min_Tier2_CpGs = 100 * mean(n2 > 0 & n2 < MIN_CPGS),
             pct_structurally_recoverable = 100 * mean(rec),
             pct_recoverable_recovered    = 100 * mean(overlapsAny(d1[rec], d2)),
             pct_Tier1_DMR_bp_recovered   = 100 * sum(width(GenomicRanges::intersect(reduce(d1), reduce(d2)))) /
                                            sum(width(reduce(d1))))
}))
fwrite(loss_tab, file.path(opt$datadir, "TableS11C_density_loss.tsv"), sep = "\t")
print(loss_tab)

# ---- 5. DMRcate parameter sensitivity on Tier 2 --------------------------------
cat("[5/5] DMRcate parameter sensitivity on Tier 2...\n")
grid_dmr <- lapply(seq_len(nrow(PARAM_GRID)), function(i) {
  p <- PARAM_GRID[i, ]
  if (p$setting == "default") return(dmr$Tier2)   # reuse the primary results
  setNames(lapply(ALL_METHODS, function(m) run_dmrcate_setting(dml$Tier2[[m]], p, paste(p$setting, m))),
           ALL_METHODS)
})
names(grid_dmr) <- PARAM_GRID$setting
saveRDS(grid_dmr, file.path(opt$datadir, "Tier2_DMR_parameter_grid.rds"))

grid_counts <- rbindlist(lapply(names(grid_dmr), function(s) rbindlist(lapply(ALL_METHODS, function(m) {
  g <- grid_dmr[[s]][[m]]
  data.table(Setting = s, Method = m, DMRs = length(g),
             median_CpGs = if (length(g)) median(g$no.cpgs) else NA_real_,
             median_width_bp = if (length(g)) median(width(g)) else NA_real_,
             pct_hyper = if (length(g)) 100 * mean(g$meandiff > 0) else NA_real_)
}))))
grid_jacc <- rbindlist(lapply(names(grid_dmr), function(s) {
  rbindlist(lapply(combn(ALL_METHODS, 2, simplify = FALSE), function(pr) {
    a <- grid_dmr[[s]][[pr[1]]]; b <- grid_dmr[[s]][[pr[2]]]
    data.table(Setting = s, Pair = paste(pr[1], "vs.", pr[2]),
               Class = if ("EPIC" %in% pr) "EPIC vs. sequencing" else "sequencing vs. sequencing",
               Jaccard_bp = jaccard_bp(a, b), Jaccard_CpG = jaccard_cpg(a, b, t2_gr))
  }))
}))
fwrite(grid_counts, file.path(opt$datadir, "TableS11D_Tier2_grid_counts.tsv"), sep = "\t")
fwrite(grid_jacc,   file.path(opt$datadir, "TableS11E_Tier2_grid_jaccard.tsv"), sep = "\t")
print(grid_jacc[, .(min_bp = min(Jaccard_bp, na.rm = TRUE), max_bp = max(Jaccard_bp, na.rm = TRUE),
                    min_cpg = min(Jaccard_CpG, na.rm = TRUE), max_cpg = max(Jaccard_CpG, na.rm = TRUE)),
                by = .(Setting, Class)])

# ---- Supplementary Figure 20 ------------------------------------------------
th <- theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5))

nn_df <- rbindlist(lapply(names(sp), function(k) data.table(Set = sp[[k]]$summary$Set,
                                                           nn = sp[[k]]$per_cpg$nn_dist)))[is.finite(nn) & nn > 0]
pA <- ggplot(nn_df, aes(nn, colour = Set)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = c(SIGMA, LAMBDA), linetype = c("dotted", "dashed")) +
  scale_x_log10(labels = scales::comma) + scale_y_continuous(labels = scales::percent) +
  labs(title = "Distance to nearest tested CpG", x = "bp (log10)", y = "Cumulative fraction", colour = NULL) +
  th + theme(legend.position = "bottom", legend.direction = "vertical")

pB <- ggplot(spacing_summary, aes(Set, pct_DMR_eligible)) + geom_col(fill = "grey40", width = 0.6) +
  labs(title = sprintf("CpGs in clusters of >= %d CpGs (gaps < %g bp)", MIN_CPGS, LAMBDA),
       x = NULL, y = "% of tested CpGs") + th + theme(axis.text.x = element_text(angle = 20, hjust = 1))

pC <- dmc_tab |>
  pivot_longer(c(pct_DMC_eligible, pct_DMC_in_DMR), names_to = "metric", values_to = "pct") |>
  mutate(metric = recode(metric, pct_DMC_eligible = "in DMR-eligible cluster", pct_DMC_in_DMR = "in called DMR"),
         Method = factor(Method, levels = ALL_METHODS)) |>
  ggplot(aes(Method, pct, fill = metric)) + geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ Tier, scales = "free_x") + scale_fill_manual(values = c("grey70", "grey25")) +
  labs(title = "DMCs eligible for / contained in DMRs", x = NULL, y = "% of DMCs", fill = NULL) +
  th + theme(legend.position = "bottom")

pD <- grid_jacc |>
  mutate(Setting = factor(Setting, levels = PARAM_GRID$setting)) |>
  ggplot(aes(Setting, Jaccard_bp, colour = Class)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c("EPIC vs. sequencing" = unname(MCOL["EPIC"]),
                                 "sequencing vs. sequencing" = "grey30")) +
  labs(title = "Tier 2 DMR overlap across DMRcate settings", x = NULL,
       y = "Base-pair-weighted Jaccard", colour = NULL) +
  th + theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(opt$outdir, "SuppFig20_dmr_probe_density.png"),
       (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A"), width = 18, height = 14, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(opt$datadir, "sessionInfo.txt"))
cat("\nDone. Tables: ", normalizePath(opt$datadir), "\n", sep = "")
