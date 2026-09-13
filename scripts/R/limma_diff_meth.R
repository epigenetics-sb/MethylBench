#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Limma Differential Methylation Analysis
# =============================================================================
# Description:
#   Runs limma-based differential methylation analysis (Blood vs Fibroblast)
#   for each platform independently. Generates one output CSV per method in
#   the format expected by 06_differential_methylation.R for the UpSet plot.
#
#   Output files per method (in --outdir):
#     EPIC_Blood_vs_Fibroblast.csv
#     ONT_Blood_vs_Fibroblast.csv
#     TWIST_Blood_vs_Fibroblast.csv
#     RRBS_Blood_vs_Fibroblast.csv
#     WGEC_Blood_vs_Fibroblast.csv
#
#   Output columns (no header, positional):
#     V1 = CpG identifier (chr:start)
#     V2 = logFC (log2 fold change, Blood vs Fibro)
#     V3 = AveExpr (average expression/methylation)
#     V4 = t statistic
#     V5 = P.Value
#     V6 = adj.P.Val (BH-corrected FDR) 
#     V7 = B (log-odds)
#     V8 = delta_beta (mean_blood - mean_fibro, on beta scale)
#
#   Additionally writes:
#     all_methods_limma_combined.csv 
#
# CHANGELOG (post-review fix, Reviewer 1 / Major comment 2, pairing part
# ONLY -- see note below on scope):
#   FIXED missing subject-level pairing:
#     Blood and fibroblast samples come from the SAME 5 individuals (paired
#     design), but the design matrix used to be group-only (~group),
#     ignoring which blood sample belongs to which fibroblast sample.
#     Subject identity is now parsed directly from the sample suffix in the
#     column names actually used for each method (get_subject_id()) --
#     deriving it from the real column names rather than assuming
#     blood_cols[i]/fibro_cols[i] line up positionally, the same principle
#     used to fix the PCA label desync in 11_pca.R -- and the design is now
#     ~subject + group, so the Blood-vs-Fibro test is adjusted for each
#     individual's baseline methylation level.
#
#   NOT CHANGED HERE, ON PURPOSE -- beta- vs. M-value scale:
#     This script implements the EXPLORATORY differential methylation stage
#     (Section 2.5.1 in the manuscript: limma and Wilcoxon applied to the
#     full overlapping CpG set, n=146,704, on beta-values, for a like-for-
#     like comparison of the two statistical frameworks across all five
#     platforms). The Methods text's statement that "EPIC was analyzed
#     using limma on M-values" refers specifically to the PRIMARY analysis
#     (Section 2.5.2: DSS Beta-Binomial for sequencing methods on the
#     Tier1/Tier2 consensus set, with limma-on-M-values as the necessary
#     EPIC counterpart, since EPIC yields no read counts for DSS). Applying
#     an M-value transform to EPIC only in THIS script would (a) not match
#     what generated the already-reported Fig. 7A exploratory numbers, and
#     (b) break the intentional same-scale comparison between limma and
#     Wilcoxon in Section 2.5.1. The beta-vs-M-value scale question is
#     therefore intentionally left untouched here; see the response letter
#     for how the primary-analysis (Tier1/Tier2) EPIC M-value step is
#     addressed.
#
# Input:
#   --all_path    Path to ALL_with_EPIC.csv
#                 Column naming convention:
#                   EPIC_Blood1..5,  EPIC_Fibro1..5
#                   ONT_Blood1..5,   WGEC_Blood1..5,  TWIST_Blood1..5, RRBS_Blood1..5
#                   ONT_Fibro1..5,   WGEC_Fibro1..5,  TWIST_Fibro1..5, RRBS_Fibro1..5
#   --outdir      Output directory for per-method CSV files
#   --fdr_cutoff  BH-adjusted p-value threshold for significance [default: 0.05]
#   --delta_cutoff Absolute delta-beta threshold [default: 0.1]
#
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(limma)
  library(dplyr)
  library(tibble)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path",
    type    = "character",
    help    = "Path to ALL_with_EPIC.csv [required]",
    metavar = "FILE"
  ),
  make_option("--outdir",
    type    = "character",
    default = "data/diff_meth/",
    help    = "Output directory for limma result files [default: data/diff_meth/]",
    metavar = "DIR"
  ),
  make_option("--fdr_cutoff",
    type    = "double",
    default = 0.05,
    help    = "BH-adjusted FDR threshold [default: 0.05]",
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

if (is.null(opt$all_path)) stop("ERROR: --all_path is required")
if (!file.exists(opt$all_path)) stop(paste("File not found:", opt$all_path))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

FDR_CUTOFF   <- opt$fdr_cutoff
DELTA_CUTOFF <- opt$delta_cutoff
BLOOD_IDS    <- paste0("Blood", 1:5)
FIBRO_IDS    <- paste0("Fibro",  1:5)

#' Parse the subject/individual ID directly from a sample column name
#' (e.g. "EPIC_Blood3" -> "3", "WGEC_Fibro3" -> "3"). Blood/Fibro pairs
#' with the same trailing number are assumed to come from the same
#' individual (matched blood/fibroblast sampling per subject; see Methods,
#' "Matched blood and fibroblast samples from five human individuals").
#'
#' Deriving this from the actual column names -- rather than assuming
#' blood_cols[i] and fibro_cols[i] line up positionally -- means pairing
#' cannot silently break if columns are ever reordered or subsetted
#' upstream (the same principle as the PCA label-desync fix in 11_pca.R).
get_subject_id <- function(col_names) {
  m  <- regmatches(col_names,
                    regexec("(?:Blood|Fibro)([0-9]+)$", col_names))
  ok <- lengths(m) == 2
  if (!all(ok)) {
    stop("get_subject_id(): could not parse subject ID from column(s): ",
         paste(col_names[!ok], collapse = ", "))
  }
  vapply(m, `[[`, character(1), 2)
}

cat("[1/4] Loading ALL matrix...\n")

all <- fread(opt$all_path, header = TRUE, sep = ",", na.strings = "NA")
cat(sprintf("  %d CpGs x %d columns\n", nrow(all), ncol(all)))

get_method_cols <- function(prefix, blood_ids, fibro_ids, all_cols) {
  blood <- intersect(paste0(prefix, "_", blood_ids), all_cols)
  fibro <- intersect(paste0(prefix, "_", fibro_ids), all_cols)
  list(blood = blood, fibro = fibro)
}

all_cols <- colnames(all)

method_defs <- list(
  EPIC  = get_method_cols("EPIC",  BLOOD_IDS, FIBRO_IDS, all_cols),
  ONT   = get_method_cols("ONT",   BLOOD_IDS, FIBRO_IDS, all_cols),
  WGEC  = get_method_cols("WGEC",  BLOOD_IDS, FIBRO_IDS, all_cols),
  TWIST = get_method_cols("TWIST", BLOOD_IDS, FIBRO_IDS, all_cols),
  RRBS  = get_method_cols("RRBS",  BLOOD_IDS, FIBRO_IDS, all_cols)
)

for (m in names(method_defs)) {
  nd <- method_defs[[m]]
  n_blood <- length(nd$blood)
  n_fibro <- length(nd$fibro)
  if (n_blood == 0 || n_fibro == 0) {
    warning(sprintf("Method %s: %d blood cols, %d fibro cols found – check column names",
      m, n_blood, n_fibro))
  } else {
    cat(sprintf("  %-5s : %d blood + %d fibro columns\n", m, n_blood, n_fibro))
  }
}

cat("\n[2/4] Running limma per method...\n")

cpg_ids <- if ("coord" %in% all_cols) {
  all$coord
} else if (all(c("chr", "start") %in% all_cols)) {
  paste0(all$chr, ":", all$start)
} else if (all(c("Chr", "Pos") %in% all_cols)) {
  paste0(all$Chr, ":", all$Pos)
} else {
  paste0("CpG_", seq_len(nrow(all)))
}

results_list <- list()

for (method_label in names(method_defs)) {

  nd         <- method_defs[[method_label]]
  blood_cols <- nd$blood
  fibro_cols <- nd$fibro

  if (length(blood_cols) == 0 || length(fibro_cols) == 0) {
    warning(sprintf("Skipping %s – no columns found", method_label))
    next
  }

  cat(sprintf("\n  [%s]\n", method_label))

  beta_mat <- as.matrix(all[, c(blood_cols, fibro_cols), with = FALSE])
  rownames(beta_mat) <- cpg_ids

  # Remove rows with any NA (limma requires complete data)
  complete_rows <- complete.cases(beta_mat)
  beta_mat      <- beta_mat[complete_rows, , drop = FALSE]
  cat(sprintf("    CpGs after NA removal: %d\n", nrow(beta_mat)))

  if (nrow(beta_mat) == 0) {
    warning(sprintf("  %s: no complete CpGs – skipping", method_label))
    next
  }

  # Design matrix: Blood = 1, Fibro = 0, ADJUSTED for subject identity
  # (paired design: the same 5 individuals contribute both a blood and a
  # fibroblast sample -- see Methods). Subject IDs are parsed from the
  # actual column names (get_subject_id()), not assumed from column
  # position, so pairing stays correct even if blood_cols/fibro_cols were
  # ever reordered upstream.
  n_blood     <- length(blood_cols)
  n_fibro     <- length(fibro_cols)
  subject_ids <- get_subject_id(c(blood_cols, fibro_cols))
  subject     <- factor(subject_ids, levels = sort(unique(subject_ids)))

  group     <- factor(c(rep("Blood", n_blood), rep("Fibro", n_fibro)),
                      levels = c("Fibro", "Blood"))
  design    <- model.matrix(~ subject + group)
  colnames(design) <- c(
    "Intercept",
    paste0("Subject_", levels(subject)[-1]),
    "Blood_vs_Fibro"
  )

  fit  <- lmFit(beta_mat, design)
  fit2 <- eBayes(fit)

  top <- topTable(fit2,
    coef       = "Blood_vs_Fibro",
    number     = Inf,
    adjust.method = "BH",
    sort.by    = "P"
  )
  top <- as.data.frame(top)
  top <- rownames_to_column(top, var = "CpG")

  mean_blood  <- rowMeans(beta_mat[top$CpG, seq_len(n_blood),  drop = FALSE], na.rm = TRUE)
  mean_fibro  <- rowMeans(beta_mat[top$CpG, (n_blood + 1):(n_blood + n_fibro), drop = FALSE], na.rm = TRUE)
  top$delta_beta <- mean_blood - mean_fibro
  top$Method     <- method_label
  top$Significant <- (!is.na(top$adj.P.Val) &
                       top$adj.P.Val < FDR_CUTOFF &
                       abs(top$delta_beta) > DELTA_CUTOFF)

  n_sig <- sum(top$Significant, na.rm = TRUE)
  cat(sprintf("    Significant DMCs (FDR < %.2f, |Δβ| > %.2f): %d\n",
    FDR_CUTOFF, DELTA_CUTOFF, n_sig))

  results_list[[method_label]] <- top

  # ----- Write per-method file in format expected by 06_differential_methylation.R
  # Columns: CpG, logFC, AveExpr, t, P.Value, adj.P.Val (=FDR), B, delta_beta
  # No header (V1..V8 positional access in downstream script)
  out_cols <- c("CpG", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "delta_beta")
  out_cols <- intersect(out_cols, colnames(top))
  out_df   <- top[, out_cols, drop = FALSE]

  file_prefix <- method_label
  out_file    <- file.path(opt$outdir,
                   paste0(file_prefix, "_Blood_vs_Fibroblast.csv"))

  write.table(out_df,
    file      = out_file,
    row.names = FALSE,
    col.names = FALSE,
    sep       = ",",
    quote     = FALSE
  )
  cat(sprintf("    Written: %s\n", basename(out_file)))
}

cat("\n[3/4] Writing combined results...\n")

all_results <- bind_rows(results_list)

combined_file <- file.path(opt$outdir, "all_methods_limma_combined.csv")
fwrite(as.data.table(all_results),
  combined_file,
  sep   = ",",
  quote = FALSE
)

cat("\n[4/4] Summary of significant DMCs per method:\n")
cat(sprintf("  %-8s  %8s  %8s\n", "Method", "Total", "Sig DMCs"))
cat(sprintf("  %s\n", strrep("-", 28)))

for (m in names(results_list)) {
  df     <- results_list[[m]]
  n_tot  <- nrow(df)
  n_sig  <- sum(df$Significant, na.rm = TRUE)
  cat(sprintf("  %-8s  %8d  %8d\n", m, n_tot, n_sig))
}

cat(sprintf("\nDone. Results written to: %s\n", opt$outdir))
