#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – DMR Analysis using DMRcate
# =============================================================================
# Description:
#   Identifies differentially methylated regions (DMRs) between blood and
#   fibroblast samples using DMRcate.
#
#   The analysis follows the existing MethylBench differential-methylation
#   structure and is kept separate from the DSS regional analysis:
#
#     - Tier 1: CpGs passing the coverage criterion in all sequencing methods
#                (ONT, TWIST, WGEC, RRBS) in both tissues.
#     - Tier 2: Tier 1 CpGs additionally covered by EPIC in both tissues.
#     - Blood vs fibroblast DMLs are obtained with DSS::DMLtest().
#     - Regional DMRs are called with DMRcate.
#     - EPIC is handled separately because it does not provide methylated /
#       unmethylated read counts required for the DSS input.
#     - DMR tables and pairwise DMR Jaccard matrices are exported.
#     - DMR count, width and Jaccard figures are generated.
#
# Input:
#   --all_path      Path to the combined methylation matrix
#                   (EPIC + sequencing methods).
#
#   Column naming convention:
#     EPIC_Blood1..n, EPIC_Fibro1..n
#     ONT_Blood1..n, ONT_Fibro1..n
#     TWIST_Blood1..n, TWIST_Fibro1..n
#     WGEC_Blood1..n, WGEC_Fibro1..n
#     RRBS_Blood1..n, RRBS_Fibro1..n
#     ONT_cov_Blood1..n, ...
#     WGEC_cov_Blood1..n, ...
#
# Output:
#   --outdir        Directory for figures.
#   --datadir       Directory for DMRcate tables and intermediate RDS files.
#
# Usage:
#   Rscript scripts/R/15_DMR_DMRcate_analysis.R \
#     --all_path    data/matrices/ALL_with_EPIC.csv \
#     --outdir      results/figures/ \
#     --datadir     results/dmr_dmrcate/
#
# Author: Lukas Laufer
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(DMRcate)
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
    type = "character",
    default = "results/figures/",
    help = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"
  ),
  make_option("--datadir",
    type = "character",
    default = "results/dmr_dmrcate/",
    help = "Output directory for DMRcate results [default: results/dmr_dmrcate/]",
    metavar = "DIR"
  ),
  make_option("--min_cov",
    type = "double",
    default = 10,
    help = "Minimum coverage per CpG [default: 10]",
    metavar = "FLOAT"
  ),
  make_option("--min_samples",
    type = "integer",
    default = 4,
    help = "Minimum number of samples passing coverage [default: 4]",
    metavar = "INT"
  ),
  make_option("--fdr_cutoff",
    type = "double",
    default = 0.05,
    help = "FDR threshold for significant DSS DMLs [default: 0.05]",
    metavar = "FLOAT"
  ),
  make_option("--delta_cutoff",
    type = "double",
    default = 0.1,
    help = "Delta-beta threshold for significant DSS DMLs [default: 0.1]",
    metavar = "FLOAT"
  ),
  make_option("--lambda",
    type = "double",
    default = 1000,
    help = "DMRcate smoothing bandwidth lambda [default: 1000]",
    metavar = "FLOAT"
  ),
  make_option("--C",
    type = "double",
    default = 2,
    help = "DMRcate scaling parameter C [default: 2]",
    metavar = "FLOAT"
  ),
  make_option("--min_cpgs",
    type = "integer",
    default = 3,
    help = "Minimum CpGs per DMR [default: 3]",
    metavar = "INT"
  ),
  make_option("--genome",
    type = "character",
    default = "hg38",
    help = "Genome build used for extractRanges [default: hg38]",
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

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

MIN_COV        <- opt$min_cov
MIN_SAMPLES    <- opt$min_samples
FDR_CUTOFF     <- opt$fdr_cutoff
DELTA_CUTOFF   <- opt$delta_cutoff
LAMBDA         <- opt$lambda
C_PARAM        <- opt$C
MIN_CPGS       <- opt$min_cpgs
GENOME         <- opt$genome

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

run_dss_tissue <- function(dml_df,
                           bs_blood,
                           bs_fibro) {

  blood_samples <- sampleNames(bs_blood)
  fibro_samples <- sampleNames(bs_fibro)

  bs_combined <- BiocGenerics::combine(
    bs_blood,
    bs_fibro
  )

  DMLtest(
    bs_combined,
    group1 = blood_samples,
    group2 = fibro_samples,
    smoothing = TRUE,
    smoothing.span = 500
  )
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

  bsseq::BSseq(
    chr = as.character(df$Chr),
    pos = as.integer(df$Pos),
    M = unname(meth_mat),
    Cov = unname(cov_mat),
    sampleNames = sample_names
  )
}

dss_to_cpgannotated <- function(dml_result) {

  df <- as.data.frame(dml_result)

  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(
      start = df$pos,
      width = 1
    ),
    stat = df$stat,
    rawpval = df$pval,
    diff = df$diff,
    ind.fdr = df$fdr,
    is.sig = df$fdr < FDR_CUTOFF &
      abs(df$diff) >= DELTA_CUTOFF
  )

  names(gr) <- paste0(df$chr, ":", df$pos)

  new(
    "CpGannotated",
    ranges = gr
  )
}

run_dmrcate <- function(dml_result,
                        method_label,
                        lambda = LAMBDA,
                        C = C_PARAM,
                        min_cpgs = MIN_CPGS) {

  cat(sprintf(
    "  [%s] DMRcate...\n",
    method_label
  ))

  cpg_ann <- dss_to_cpgannotated(dml_result)

  dmr <- tryCatch(
    dmrcate(
      cpg_ann,
      lambda = lambda,
      C = C,
      min.cpgs = min_cpgs
    ),
    error = function(e) {
      cat(sprintf(
        "    [%s] ERROR: %s\n",
        method_label,
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(dmr)) {
    return(NULL)
  }

  dmr_gr <- extractRanges(
    dmr,
    genome = GENOME
  )

  cat(sprintf(
    "    [%s] %d DMRs\n",
    method_label,
    length(dmr_gr)
  ))

  dmr_gr
}

jaccard_dmr <- function(gr1, gr2) {

  if (is.null(gr1) ||
      is.null(gr2) ||
      length(gr1) == 0 ||
      length(gr2) == 0) {
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
    vapply(
      dmr_list,
      function(x) !is.null(x) && length(x) > 0,
      logical(1)
    )
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
  file.path(
    opt$datadir,
    "Tier1_consensus_CpGs.tsv"
  ),
  sep = "\t"
)

fwrite(
  data.table(cpg_id = consensus_tier2),
  file.path(
    opt$datadir,
    "Tier2_consensus_CpGs.tsv"
  ),
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
  file.path(
    opt$datadir,
    "BSseq_Tier1.rds"
  )
)

saveRDS(
  bsseq_t2,
  file.path(
    opt$datadir,
    "BSseq_Tier2.rds"
  )
)

# -------------------------------------------------------------------------
# 4. DSS DML + DMRcate analysis
# -------------------------------------------------------------------------

cat("[4/6] Running DSS DML and DMRcate analysis...\n")

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
      bs_blood = bsseq_list[[key_blood]],
      bs_fibro = bsseq_list[[key_fibro]]
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
        paste0(
          tier_label,
          "_DML_significant_",
          method,
          ".tsv"
        )
      ),
      sep = "\t"
    )

    fwrite(
      as.data.table(dml),
      file.path(
        opt$datadir,
        paste0(
          tier_label,
          "_DML_",
          method,
          ".tsv"
        )
      ),
      sep = "\t"
    )

    dmr <- run_dmrcate(
      dml,
      method
    )

    dmr_list[[method]] <- dmr

    if (!is.null(dmr) && length(dmr) > 0) {

      dmr_df <- as.data.frame(dmr)

      fwrite(
        as.data.table(dmr_df),
        file.path(
          opt$datadir,
          paste0(
            tier_label,
            "_DMR_",
            method,
            ".tsv"
          )
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

dmrcate_t1 <- run_tier(
  bsseq_t1,
  "Tier1"
)

dmrcate_t2 <- run_tier(
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

    if (is.null(dmr) || length(dmr) == 0) {
      return(data.frame(
        Tier = tier_label,
        Method = method,
        DMRs = 0L,
        MedianWidth = NA_real_,
        MedianCpGs = NA_real_
      ))
    }

    dmr_df <- as.data.frame(dmr)

    cpg_col <- intersect(
      c("no.cpgs", "nCG", "n_cpgs"),
      colnames(dmr_df)
    )

    data.frame(
      Tier = tier_label,
      Method = method,
      DMRs = nrow(dmr_df),
      MedianWidth = median(
        width(dmr),
        na.rm = TRUE
      ),
      MedianCpGs = if (length(cpg_col) > 0) {
        median(
          dmr_df[[cpg_col[1]]],
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
    )
  }))
}

summary_t1 <- make_dmr_summary(
  dmrcate_t1$dmr,
  "Tier1"
)

summary_t2 <- make_dmr_summary(
  dmrcate_t2$dmr,
  "Tier2"
)

dmr_summary <- bind_rows(
  summary_t1,
  summary_t2
)

fwrite(
  as.data.table(dmr_summary),
  file.path(
    opt$datadir,
    "DMRcate_DMR_summary.tsv"
  ),
  sep = "\t"
)

plot_dmr_counts <- function(dmr_list,
                            tier_label,
                            method_col) {

  count_df <- bind_rows(lapply(
    names(dmr_list),
    function(method) {

      dmr <- dmr_list[[method]]

      if (is.null(dmr) ||
          length(dmr) == 0 ||
          !method %in% names(method_col)) {
        return(NULL)
      }

      dmr_df <- as.data.frame(dmr)

      diff_col <- intersect(
        c("diff.Methy", "meanMethyDiff"),
        colnames(dmr_df)
      )

      if (length(diff_col) == 0) {
        return(data.frame(
          method = method,
          n_hyper = NA_integer_,
          n_hypo = NA_integer_
        ))
      }

      diff_values <- dmr_df[[diff_col[1]]]

      data.frame(
        method = method,
        n_hyper = sum(
          diff_values > 0,
          na.rm = TRUE
        ),
        n_hypo = sum(
          diff_values < 0,
          na.rm = TRUE
        )
      )
    }
  )) |>
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
        "DMRcate DMR Counts per Method –",
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
    "dmr_DMRcate_counts_Tier1.png"
  ),
  plot_dmr_counts(
    dmrcate_t1$dmr,
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
    "dmr_DMRcate_counts_Tier2.png"
  ),
  plot_dmr_counts(
    dmrcate_t2$dmr,
    "Tier2",
    METHOD_COLORS
  ),
  width = 9,
  height = 7,
  dpi = 300
)

plot_dmr_width <- function(dmr_list,
                           tier_label,
                           method_col) {

  width_df <- bind_rows(lapply(
    names(dmr_list),
    function(method) {

      dmr <- dmr_list[[method]]

      if (is.null(dmr) ||
          length(dmr) == 0 ||
          !method %in% names(method_col)) {
        return(NULL)
      }

      data.frame(
        method = method,
        width_bp = width(dmr)
      )
    }
  )) |>
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
        "DMRcate DMR Width Distribution –",
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
    "dmr_DMRcate_width_Tier1.png"
  ),
  plot_dmr_width(
    dmrcate_t1$dmr,
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
    "dmr_DMRcate_width_Tier2.png"
  ),
  plot_dmr_width(
    dmrcate_t2$dmr,
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
  dmrcate_t1$dmr
)

jaccard_dmr_t2 <- build_jaccard_mat(
  dmrcate_t2$dmr
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
    "DMRcate_DMR_Jaccard_Tier1.tsv"
  )
)

write_jaccard(
  jaccard_dmr_t2,
  file.path(
    opt$datadir,
    "DMRcate_DMR_Jaccard_Tier2.tsv"
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

  jmat <- jmat[
    keys,
    keys,
    drop = FALSE
  ]

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
      "Pairwise DMRcate DMR Jaccard Index –",
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
        "jaccard_DMR_DMRcate_",
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
cat(
  "  DMR tables : ",
  normalizePath(opt$datadir),
  "\n",
  sep = ""
)
cat(
  "  Figures    : ",
  normalizePath(opt$outdir),
  "\n",
  sep = ""
)
