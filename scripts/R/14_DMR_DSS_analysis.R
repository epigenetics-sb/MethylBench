#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – DMR Analysis using DSS
# =============================================================================
# Description:
#   Identifies differentially methylated regions (DMRs) between blood and
#   fibroblast samples using the DSS beta-binomial framework.
#
#   The analysis follows the existing MethylBench differential-methylation
#   structure, but separates regional DSS analysis from DMC detection and
#   annotation:
#
#     - Tier 1: CpGs passing the coverage criterion in all sequencing methods
#                (ONT, TWIST, WGEC, RRBS) in both tissues.
#     - Tier 2: Tier 1 CpGs additionally covered by EPIC in both tissues.
#                DSS DMRs are still called only for sequencing methods because
#                EPIC does not provide methylated/unmethylated read counts.
#     - Blood vs fibroblast comparison is performed with DSS::DMLtest().
#     - Regional DMRs are called with DSS::callDMR().
#     - DMR tables and pairwise DMR Jaccard matrices are exported.
#     - DMR count, width and Jaccard figures are generated.
#
# Input:
#   --all_path      Path to the combined methylation matrix
#                   (EPIC + sequencing methods).
#
#   Column naming convention:
#     EPIC_Blood1..n, EPIC_Fibro1..n
#     ONT_Blood1..n,  ONT_Fibro1..n
#     TWIST_Blood1..n, TWIST_Fibro1..n
#     WGEC_Blood1..n, WGEC_Fibro1..n
#     RRBS_Blood1..n, RRBS_Fibro1..n
#     ONT_cov_Blood1..n, ...
#     WGEC_cov_Blood1..n, ...
#
# Output:
#   --outdir        Directory for figures.
#   --datadir       Directory for DMR tables and intermediate RDS files.
#
# Usage:
#   Rscript scripts/R/14_DMR_DSS_analysis.R \
#     --all_path    data/matrices/ALL_with_EPIC.csv \
#     --outdir      results/figures/ \
#     --datadir     results/dmr_dss/
#
# Author: Lukas Laufer
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
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path",
    type = "character",
    help = "Path to combined methylation matrix [required]",
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
    default = "results/dmr_dss/",
    help    = "Output directory for DSS DMR results [default: results/dmr_dss/]",
    metavar = "DIR"
  ),
  make_option("--min_cov",
    type    = "double",
    default = 10,
    help    = "Minimum coverage per CpG [default: 10]",
    metavar = "FLOAT"
  ),
  make_option("--min_samples",
    type    = "integer",
    default = 4,
    help    = "Minimum number of samples passing coverage [default: 4]",
    metavar = "INT"
  ),
  make_option("--fdr_cutoff",
    type    = "double",
    default = 0.05,
    help    = "FDR threshold used for DML reporting [default: 0.05]",
    metavar = "FLOAT"
  ),
  make_option("--delta_cutoff",
    type    = "double",
    default = 0.1,
    help    = "Delta-beta threshold used for DMR calling [default: 0.1]",
    metavar = "FLOAT"
  ),
  make_option("--p_threshold",
    type    = "double",
    default = 1e-5,
    help    = "Raw p-value threshold used by DSS::callDMR [default: 1e-5]",
    metavar = "FLOAT"
  ),
  make_option("--smoothing_span",
    type    = "integer",
    default = 500,
    help    = "DSS smoothing span in bp [default: 500]",
    metavar = "INT"
  ),
  make_option("--minlen",
    type    = "integer",
    default = 50,
    help    = "Minimum DMR length in bp [default: 50]",
    metavar = "INT"
  ),
  make_option("--min_cpgs",
    type    = "integer",
    default = 3,
    help    = "Minimum number of CpGs per DMR [default: 3]",
    metavar = "INT"
  ),
  make_option("--dis_merge",
    type    = "integer",
    default = 100,
    help    = "Maximum distance for merging nearby DMRs [default: 100]",
    metavar = "INT"
  ),
  make_option("--pct_sig",
    type    = "double",
    default = 0.5,
    help    = "Minimum fraction of significant CpGs per DMR [default: 0.5]",
    metavar = "FLOAT"
  ),
  make_option("--genome",
    type    = "character",
    default = "hg38",
    help    = "Genome build used for output metadata [default: hg38]",
    metavar = "STRING"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$all_path)) {
  stop("ERROR: --all_path is required")
}
if (!file.exists(opt$all_path)) {
  stop(paste("File not found:", opt$all_path))
}

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

MIN_COV       <- opt$min_cov
MIN_SAMPLES   <- opt$min_samples
FDR_CUTOFF    <- opt$fdr_cutoff
DELTA_CUTOFF  <- opt$delta_cutoff
P_THRESHOLD   <- opt$p_threshold
SMOOTHING_SPAN <- opt$smoothing_span
MINLEN        <- opt$minlen
MIN_CPGS      <- opt$min_cpgs
DIS_MERGE     <- opt$dis_merge
PCT_SIG       <- opt$pct_sig
GENOME        <- opt$genome

METHOD_PREFIX <- c(
  ONT   = "ONT",
  TWIST = "TWIST",
  WGEC  = "WGEC",
  RRBS  = "RRBS"
)

METHOD_COLORS <- c(
  ONT   = "#E69F00",
  TWIST = "#009E73",
  WGEC  = "#5654E9",
  RRBS  = "#0072B2"
)

TISSUES <- c("Blood", "Fibro")

cat("[1/6] Loading methylation matrix...\n")
combined_df <- fread(
  opt$all_path,
  header = TRUE,
  sep = ",",
  na.strings = "NA"
)

required_coord_cols <- c("Chr", "Pos")
missing_coord_cols <- setdiff(required_coord_cols, colnames(combined_df))
if (length(missing_coord_cols) > 0) {
  stop(
    "Required coordinate columns not found: ",
    paste(missing_coord_cols, collapse = ", ")
  )
}

combined_df[, cpg_id := paste0(Chr, ":", Pos)]

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

get_meth_cols <- function(df, method, tissue) {
  prefix <- METHOD_PREFIX[[method]]
  grep(
    paste0("^", prefix, "_", tissue, "[0-9]+$"),
    colnames(df),
    value = TRUE
  )
}

get_cov_cols <- function(df, method, tissue) {
  prefix <- METHOD_PREFIX[[method]]
  grep(
    paste0("^", prefix, "_cov_", tissue, "[0-9]+$"),
    colnames(df),
    value = TRUE
  )
}

get_passing_cpgs_seq <- function(df,
                                 method,
                                 tissue,
                                 min_cov = MIN_COV,
                                 min_samples = MIN_SAMPLES) {

  cov_cols <- get_cov_cols(df, method, tissue)

  if (length(cov_cols) == 0) {
    stop(sprintf(
      "[%s | %s] No coverage columns found.",
      method, tissue
    ))
  }

  cov_mat <- as.matrix(df[, ..cov_cols])
  storage.mode(cov_mat) <- "numeric"

  pass <- rowSums(
    !is.na(cov_mat) & cov_mat >= min_cov
  ) >= min_samples

  df$cpg_id[pass]
}

get_passing_cpgs_epic <- function(df,
                                  tissue,
                                  min_samples = MIN_SAMPLES) {

  beta_cols <- grep(
    paste0("^EPIC_", tissue, "[0-9]+$"),
    colnames(df),
    value = TRUE
  )

  if (length(beta_cols) == 0) {
    stop(sprintf(
      "[EPIC | %s] No EPIC columns found.",
      tissue
    ))
  }

  beta_mat <- as.matrix(df[, ..beta_cols])

  pass <- rowSums(!is.na(beta_mat)) >= min_samples

  df$cpg_id[pass]
}

make_bsseq_from_combined <- function(df, method, tissue) {

  beta_cols <- get_meth_cols(df, method, tissue)
  cov_cols  <- get_cov_cols(df, method, tissue)

  if (length(beta_cols) == 0) {
    stop(sprintf(
      "[%s | %s] No methylation columns found.",
      method, tissue
    ))
  }

  if (length(cov_cols) == 0) {
    stop(sprintf(
      "[%s | %s] No coverage columns found.",
      method, tissue
    ))
  }

  if (length(beta_cols) != length(cov_cols)) {
    stop(sprintf(
      "[%s | %s] Number of methylation columns (%d) does not match \
number of coverage columns (%d).",
      method, tissue, length(beta_cols), length(cov_cols)
    ))
  }

  cov_mat <- as.matrix(df[, ..cov_cols])
  storage.mode(cov_mat) <- "numeric"

  beta_mat <- as.matrix(df[, ..beta_cols])
  storage.mode(beta_mat) <- "numeric"

  meth_mat <- round(beta_mat * cov_mat)

  na_mask <- is.na(cov_mat) | cov_mat <= 0

  cov_mat[na_mask] <- 0
  meth_mat[is.na(meth_mat) | na_mask] <- 0

  sample_names <- sub(
    paste0("^", METHOD_PREFIX[[method]], "_"),
    paste0(method, "_"),
    beta_cols
  )

  BSseq(
    chr         = as.character(df$Chr),
    pos         = as.integer(df$Pos),
    M           = unname(meth_mat),
    Cov         = unname(cov_mat),
    sampleNames = sample_names
  )
}

run_dss_tissue <- function(bs_blood,
                           bs_fibro,
                           smoothing_span = SMOOTHING_SPAN) {

  blood_samples <- sampleNames(bs_blood)
  fibro_samples <- sampleNames(bs_fibro)

  bs_combined <- BiocGenerics::combine(
    bs_blood,
    bs_fibro
  )

  DMLtest(
    bs_combined,
    group1         = blood_samples,
    group2         = fibro_samples,
    smoothing      = TRUE,
    smoothing.span = smoothing_span
  )
}

call_dss_dmr <- function(dml_result,
                         delta = DELTA_CUTOFF,
                         p.threshold = P_THRESHOLD,
                         minlen = MINLEN,
                         minCG = MIN_CPGS,
                         dis.merge = DIS_MERGE,
                         pct.sig = PCT_SIG) {

  callDMR(
    dml_result,
    delta     = delta,
    p.threshold = p.threshold,
    minlen    = minlen,
    minCG     = minCG,
    dis.merge = dis.merge,
    pct.sig   = pct.sig
  )
}

jaccard_dmr <- function(gr1, gr2) {

  if (is.null(gr1) || is.null(gr2) ||
      length(gr1) == 0 || length(gr2) == 0) {
    return(NA_real_)
  }

  gr1 <- GenomicRanges::reduce(gr1)
  gr2 <- GenomicRanges::reduce(gr2)

  inter <- GenomicRanges::intersect(gr1, gr2)
  union <- GenomicRanges::reduce(c(gr1, gr2))

  if (length(union) == 0) {
    return(NA_real_)
  }

  sum(width(inter)) / sum(width(union))
}

build_jaccard_mat <- function(dmr_list) {

  keys <- names(dmr_list)[
    vapply(dmr_list, function(x) {
      !is.null(x) && length(x) > 0
    }, logical(1))
  ]

  if (length(keys) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  mat <- outer(
    keys,
    keys,
    Vectorize(function(i, j) {
      jaccard_dmr(
        dmr_list[[i]],
        dmr_list[[j]]
      )
    })
  )

  rownames(mat) <- colnames(mat) <- keys

  mat
}

# -------------------------------------------------------------------------
# 2. Determine consensus CpG sets
# -------------------------------------------------------------------------

cat("[2/6] Determining consensus CpG sets...\n")

passing_cpgs <- list()

for (method in names(METHOD_PREFIX)) {

  blood_pass <- get_passing_cpgs_seq(
    combined_df,
    method,
    "Blood"
  )

  fibro_pass <- get_passing_cpgs_seq(
    combined_df,
    method,
    "Fibro"
  )

  both_pass <- intersect(
    blood_pass,
    fibro_pass
  )

  cat(sprintf(
    "  [%s] Blood: %d | Fibro: %d | Blood ∩ Fibro: %d\n",
    method,
    length(blood_pass),
    length(fibro_pass),
    length(both_pass)
  ))

  passing_cpgs[[method]] <- both_pass
}

epic_blood_pass <- get_passing_cpgs_epic(
  combined_df,
  "Blood"
)

epic_fibro_pass <- get_passing_cpgs_epic(
  combined_df,
  "Fibro"
)

epic_both_pass <- intersect(
  epic_blood_pass,
  epic_fibro_pass
)

cat(sprintf(
  "  [EPIC] Blood: %d | Fibro: %d | Blood ∩ Fibro: %d\n",
  length(epic_blood_pass),
  length(epic_fibro_pass),
  length(epic_both_pass)
))

passing_cpgs[["EPIC"]] <- epic_both_pass

consensus_tier1 <- Reduce(
  intersect,
  passing_cpgs[names(METHOD_PREFIX)]
)

consensus_tier2 <- intersect(
  consensus_tier1,
  passing_cpgs[["EPIC"]]
)

cat(sprintf(
  "  Tier 1 consensus: %d CpGs\n",
  length(consensus_tier1)
))

cat(sprintf(
  "  Tier 2 consensus: %d CpGs\n",
  length(consensus_tier2)
))

combined_tier1 <- combined_df[
  combined_df$cpg_id %in% consensus_tier1
]

combined_tier2 <- combined_df[
  combined_df$cpg_id %in% consensus_tier2
]

stopifnot(
  nrow(combined_tier1) == length(consensus_tier1),
  nrow(combined_tier2) == length(consensus_tier2)
)

fwrite(
  data.table(cpg_id = consensus_tier1),
  file.path(opt$datadir, "Tier1_consensus_CpGs.tsv"),
  sep = "\t"
)

fwrite(
  data.table(cpg_id = consensus_tier2),
  file.path(opt$datadir, "Tier2_consensus_CpGs.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# 3. Build BSseq objects
# -------------------------------------------------------------------------

cat("[3/6] Building BSseq objects...\n")

bsseq_t1 <- list()
bsseq_t2 <- list()

for (method in names(METHOD_PREFIX)) {
  for (tissue in TISSUES) {

    key <- paste0(method, "_", tissue)

    bsseq_t1[[key]] <- make_bsseq_from_combined(
      combined_tier1,
      method,
      tissue
    )

    bsseq_t2[[key]] <- make_bsseq_from_combined(
      combined_tier2,
      method,
      tissue
    )

    cat(sprintf(
      "  [%s] Tier1: %d CpGs | Tier2: %d CpGs | samples: %d\n",
      key,
      nrow(bsseq_t1[[key]]),
      nrow(bsseq_t2[[key]]),
      ncol(bsseq_t1[[key]])
    ))
  }
}

saveRDS(
  bsseq_t1,
  file.path(opt$datadir, "BSseq_Tier1.rds")
)

saveRDS(
  bsseq_t2,
  file.path(opt$datadir, "BSseq_Tier2.rds")
)

# -------------------------------------------------------------------------
# 4. DSS DML + DMR analysis
# -------------------------------------------------------------------------

cat("[4/6] Running DSS DML and DMR analysis...\n")

run_tier <- function(bsseq_list, tier_label) {

  dml_list <- list()
  dmr_list <- list()

  for (method in names(METHOD_PREFIX)) {

    key_blood <- paste0(method, "_Blood")
    key_fibro <- paste0(method, "_Fibro")

    cat(sprintf(
      "  [%s | %s] DMLtest...\n",
      tier_label,
      method
    ))

    dml <- run_dss_tissue(
      bsseq_list[[key_blood]],
      bsseq_list[[key_fibro]]
    )

    dml_list[[method]] <- dml

    dml_sig <- callDML(
      dml,
      p.threshold = FDR_CUTOFF,
      delta = DELTA_CUTOFF
    )

    fwrite(
      as.data.table(dml_sig),
      file.path(
        opt$datadir,
        paste0(tier_label, "_DML_significant_", method, ".tsv")
      ),
      sep = "\t"
    )

    cat(sprintf(
      "    [%s | %s] %d significant DMLs\n",
      tier_label,
      method,
      nrow(dml_sig)
    ))

    cat(sprintf(
      "  [%s | %s] callDMR...\n",
      tier_label,
      method
    ))

    dmr <- tryCatch(
      call_dss_dmr(dml),
      error = function(e) {
        cat(sprintf(
          "    [%s | %s] ERROR: %s\n",
          tier_label,
          method,
          conditionMessage(e)
        ))
        NULL
      }
    )

    dmr_list[[method]] <- dmr

    cat(sprintf(
      "    [%s | %s] %d DMRs\n",
      tier_label,
      method,
      ifelse(is.null(dmr), 0L, nrow(dmr))
    ))

    fwrite(
      as.data.table(dml),
      file.path(
        opt$datadir,
        paste0(tier_label, "_DML_", method, ".tsv")
      ),
      sep = "\t"
    )

    if (!is.null(dmr) && nrow(dmr) > 0) {
      fwrite(
        as.data.table(dmr),
        file.path(
          opt$datadir,
          paste0(tier_label, "_DMR_", method, ".tsv")
        ),
        sep = "\t"
      )
    }
  }

  saveRDS(
    dml_list,
    file.path(
      opt$datadir,
      paste0(tier_label, "_DML.rds")
    )
  )

  saveRDS(
    dmr_list,
    file.path(
      opt$datadir,
      paste0(tier_label, "_DMR.rds")
    )
  )

  list(
    dml = dml_list,
    dmr = dmr_list
  )
}

dss_t1 <- run_tier(
  bsseq_t1,
  "Tier1"
)

dss_t2 <- run_tier(
  bsseq_t2,
  "Tier2"
)

# -------------------------------------------------------------------------
# 5. DMR summaries
# -------------------------------------------------------------------------

cat("[5/6] Summarizing DMR results...\n")

make_dmr_summary <- function(dmr_list, tier_label) {

  bind_rows(lapply(names(dmr_list), function(method) {

    dmr <- dmr_list[[method]]

    if (is.null(dmr) || nrow(dmr) == 0) {
      return(data.frame(
        Tier       = tier_label,
        Method     = method,
        DMRs       = 0L,
        MedianWidth = NA_real_,
        MedianCpGs  = NA_real_
      ))
    }

    data.frame(
      Tier        = tier_label,
      Method      = method,
      DMRs        = nrow(dmr),
      MedianWidth = median(dmr$length, na.rm = TRUE),
      MedianCpGs  = median(dmr$nCG, na.rm = TRUE)
    )
  }))
}

summary_t1 <- make_dmr_summary(
  dss_t1$dmr,
  "Tier1"
)

summary_t2 <- make_dmr_summary(
  dss_t2$dmr,
  "Tier2"
)

dmr_summary <- bind_rows(
  summary_t1,
  summary_t2
)

fwrite(
  as.data.table(dmr_summary),
  file.path(opt$datadir, "DSS_DMR_summary.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# DMR count plot
# -------------------------------------------------------------------------

plot_dmr_counts <- function(dmr_list,
                            tier_label,
                            method_col) {

  count_df <- bind_rows(lapply(names(dmr_list), function(method) {

    dmr <- dmr_list[[method]]

    if (is.null(dmr) ||
        nrow(dmr) == 0 ||
        !method %in% names(method_col)) {
      return(NULL)
    }

    data.frame(
      method = method,
      n_hyper = sum(
        dmr$diff.Methy > 0,
        na.rm = TRUE
      ),
      n_hypo = sum(
        dmr$diff.Methy < 0,
        na.rm = TRUE
      )
    )
  })) |>
    pivot_longer(
      c(n_hyper, n_hypo),
      names_to = "direction",
      values_to = "count"
    ) |>
    mutate(
      direction = factor(
        direction,
        levels = c("n_hyper", "n_hypo"),
        labels = c(
          "Hypermethylated (Blood > Fibro)",
          "Hypomethylated (Blood < Fibro)"
        )
      ),
      method = factor(
        method,
        levels = names(method_col)
      )
    )

  ggplot(
    count_df,
    aes(
      x = method,
      y = count,
      fill = direction
    )
  ) +
    geom_col(width = 0.6) +
    scale_fill_manual(
      values = c(
        "Hypermethylated (Blood > Fibro)" = "#d73027",
        "Hypomethylated (Blood < Fibro)" = "#4575b4"
      )
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste(
        "DSS DMR Counts per Method –",
        tier_label
      ),
      x = NULL,
      y = "#DMRs",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        hjust = 0.5,
        size = 22
      ),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 22),
      text = element_text(size = 22),
      axis.text.x = element_text(
        colour = method_col[
          levels(count_df$method)
        ]
      )
    )
}

ggsave(
  file.path(
    opt$outdir,
    "dmr_DSS_counts_Tier1.png"
  ),
  plot_dmr_counts(
    dss_t1$dmr,
    "Tier1",
    METHOD_COLORS
  ),
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    opt$outdir,
    "dmr_DSS_counts_Tier2.png"
  ),
  plot_dmr_counts(
    dss_t2$dmr,
    "Tier2",
    METHOD_COLORS
  ),
  width = 9,
  height = 7,
  dpi = 300
)

# -------------------------------------------------------------------------
# DMR width plots
# -------------------------------------------------------------------------

plot_dmr_width <- function(dmr_list,
                           tier_label,
                           method_col) {

  width_df <- bind_rows(lapply(names(dmr_list), function(method) {

    dmr <- dmr_list[[method]]

    if (is.null(dmr) ||
        nrow(dmr) == 0 ||
        !method %in% names(method_col)) {
      return(NULL)
    }

    data.frame(
      method = method,
      width_bp = dmr$length
    )
  })) |>
    mutate(
      method = factor(
        method,
        levels = names(method_col)
      )
    )

  ggplot(
    width_df,
    aes(
      x = method,
      y = width_bp,
      fill = method
    )
  ) +
    geom_violin(
      alpha = 0.8,
      trim = TRUE,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      fill = "white",
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = method_col
    ) +
    scale_y_log10(
      labels = scales::comma
    ) +
    labs(
      title = paste(
        "DSS DMR Width Distribution –",
        tier_label
      ),
      x = NULL,
      y = "DMR width (bp, log10)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 22
      ),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 22),
      text = element_text(size = 22),
      axis.text.x = element_text(
        colour = method_col[
          levels(width_df$method)
        ]
      )
    )
}

ggsave(
  file.path(
    opt$outdir,
    "dmr_DSS_width_Tier1.png"
  ),
  plot_dmr_width(
    dss_t1$dmr,
    "Tier1",
    METHOD_COLORS
  ),
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    opt$outdir,
    "dmr_DSS_width_Tier2.png"
  ),
  plot_dmr_width(
    dss_t2$dmr,
    "Tier2",
    METHOD_COLORS
  ),
  width = 9,
  height = 7,
  dpi = 300
)

# -------------------------------------------------------------------------
# Pairwise DMR Jaccard
# -------------------------------------------------------------------------

jaccard_dmr_t1 <- build_jaccard_mat(
  dss_t1$dmr
)

jaccard_dmr_t2 <- build_jaccard_mat(
  dss_t2$dmr
)

write_jaccard <- function(mat, path) {

  if (length(mat) == 0) {
    return(invisible(NULL))
  }

  fwrite(
    as.data.table(
      mat,
      keep.rownames = "Method"
    ),
    path,
    sep = "\t"
  )
}

write_jaccard(
  jaccard_dmr_t1,
  file.path(
    opt$datadir,
    "DSS_DMR_Jaccard_Tier1.tsv"
  )
)

write_jaccard(
  jaccard_dmr_t2,
  file.path(
    opt$datadir,
    "DSS_DMR_Jaccard_Tier2.tsv"
  )
)

plot_dmr_jaccard <- function(jmat,
                             tier_label,
                             method_col) {

  if (length(jmat) == 0) {
    return(NULL)
  }

  keys <- intersect(
    rownames(jmat),
    names(method_col)
  )

  jmat <- jmat[keys, keys, drop = FALSE]

  col_fun <- colorRamp2(
    c(0, 0.5, 1),
    c("#f7fbff", "#6baed6", "#08306b")
  )

  label_colors <- method_col[keys]

  top_ann <- HeatmapAnnotation(
    Method = keys,
    col = list(
      Method = setNames(
        method_col[keys],
        keys
      )
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  left_ann <- rowAnnotation(
    Method = keys,
    col = list(
      Method = setNames(
        method_col[keys],
        keys
      )
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  ht <- Heatmap(
    jmat,
    name = "Jaccard\nIndex",
    col = col_fun,
    na_col = "grey80",
    cell_fun = function(
      j, i, x, y, width, height, fill
    ) {

      val <- jmat[i, j]

      if (!is.na(val)) {
        grid.text(
          sprintf("%.3f", val),
          x,
          y,
          gp = gpar(
            fontsize = 11,
            col = ifelse(
              val > 0.5,
              "white",
              "black"
            )
          )
        )
      }
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      col = label_colors
    ),
    column_names_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      col = label_colors
    ),
    column_title = paste(
      "Pairwise DSS DMR Jaccard Index –",
      tier_label
    ),
    column_title_gp = gpar(
      fontsize = 13,
      fontface = "bold"
    ),
    left_annotation = left_ann,
    top_annotation = top_ann
  )

  png(
    file.path(
      opt$outdir,
      paste0(
        "jaccard_DMR_DSS_",
        tier_label,
        ".png"
      )
    ),
    units = "in",
    width = 9,
    height = 7,
    res = 500
  )

  draw(ht)

  dev.off()
}

plot_dmr_jaccard(
  jaccard_dmr_t1,
  "Tier1",
  METHOD_COLORS
)

plot_dmr_jaccard(
  jaccard_dmr_t2,
  "Tier2",
  METHOD_COLORS
)

# -------------------------------------------------------------------------
# Save session information
# -------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    opt$datadir,
    "sessionInfo.txt"
  )
)

cat("\n[6/6] Done.\n")
cat("  DMR tables : ", normalizePath(opt$datadir), "\n", sep = "")
cat("  Figures    : ", normalizePath(opt$outdir), "\n", sep = "")
