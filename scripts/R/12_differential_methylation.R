#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Differential Methylation Analysis
# =============================================================================
# Description:
#   Computes differentially methylated CpGs (DMCs) between blood and fibroblast
#   samples using a Wilcoxon rank-sum test per method. Generates:
#     - Delta-beta density plot per method
#     - Heatmaps of top 500 / 1000 / 5000 variable CpGs
#     - UpSet plot (limma-based DMC overlap, from pre-computed per-method files)
#     - UpSet plot (Wilcoxon-based DMC overlap, EPIC vs TWIST)
#     - Cross-platform EPIC vs TWIST delta-beta scatter
#     - Variance boxplot per method (Blood and Fibro)
#     - Coverage boxplot per method (Blood and Fibro)
#
# Input:
#   --all_path      Path to ALL_with_EPIC.csv (EPIC + sequencing matrix)
#                   Column naming convention in this matrix:
#                     EPIC_Blood1..5, EPIC_Fibro1..5
#                     ONT_Blood1..5,  WGEC_Blood1..5, TWIST_Blood1..5, RRBS_Blood1..5
#                     ONT_Fibro1..5,  WGEC_Fibro1..5, TWIST_Fibro1..5, RRBS_Fibro1..5
#                     ONT_cov_*,      WGEC_cov_*,     TWIST_cov_*,     RRBS_cov_*
#   --blood_path    Path to Blood_without_EPIC.csv
#   --fibro_path    Path to Fibro_without_EPIC.csv
#   --limma_dir     Directory with per-method limma DMC files:
#                     EPIC_Blood_vs_Fibroblast.csv, ONT_Blood_vs_Fibroblast.csv,
#                     TWIST_Blood_vs_Fibroblast.csv, RRBS_Blood_vs_Fibroblast.csv,
#                     WGEC_Blood_vs_Fibroblast.csv
#                   Expected columns: V1=CpG, V6=FDR
#   --outdir        Output directory for figures
#   --datadir       Output directory for intermediate data files
#   --delta_cutoff  Absolute delta-beta threshold for significance [default: 0.1]
#
# Usage:
#   Rscript 06_differential_methylation.R \
#     --all_path    data/matrices/ALL_with_EPIC.csv \
#     --blood_path  data/matrices/Blood_without_EPIC.csv \
#     --fibro_path  data/matrices/Fibro_without_EPIC.csv \
#     --limma_dir   data/diff_meth/ \
#     --outdir      results/figures/ \
#     --datadir     results/diff_meth/
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(UpSetR)
  library(circlize)
  library(grid)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--all_path",
    type = "character", help = "Path to ALL_with_EPIC.csv [required]",
    metavar = "FILE"),
  make_option("--blood_path",
    type = "character", help = "Path to Blood_without_EPIC.csv [required]",
    metavar = "FILE"),
  make_option("--fibro_path",
    type = "character", help = "Path to Fibro_without_EPIC.csv [required]",
    metavar = "FILE"),
  make_option("--limma_dir",
    type = "character", help = "Directory with per-method limma DMC CSV files [required]",
    metavar = "DIR"),
  make_option("--outdir",
    type    = "character", default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"),
  make_option("--datadir",
    type    = "character", default = "results/diff_meth/",
    help    = "Output directory for data files [default: results/diff_meth/]",
    metavar = "DIR"),
  make_option("--delta_cutoff",
    type    = "double", default = 0.1,
    help    = "Absolute delta-beta threshold for significance [default: 0.1]",
    metavar = "FLOAT")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("all_path", "blood_path", "fibro_path", "limma_dir")) {
  val <- opt[[arg]]
  if (is.null(val)) stop(paste("ERROR: --", arg, " is required", sep = ""))
}

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

DELTA_CUTOFF <- opt$delta_cutoff
col.vec      <- get_colors()
BLOOD_IDS    <- paste0("Blood", 1:5)
FIBRO_IDS    <- paste0("Fibro",  1:5)

# Platform annotation for heatmap and upset
PLATFORMS <- c(
  "EPIC"  = "Array",
  "TWIST" = "ShortRead",
  "WGEC"  = "ShortRead",
  "RRBS"  = "ShortRead",
  "ONT"   = "LongRead"
)

PLATFORM_COLORS <- c(
  "Array"      = "#D5E8D4",
  "ShortRead"  = "#E1D5E7",
  "LongRead"   = "#FFE6CC"
)

cat("[1/6] Loading matrices...\n")

all   <- fread(opt$all_path,   header = TRUE, sep = ",", na.strings = "NA")
blood <- fread(opt$blood_path, header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(opt$fibro_path, header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  ALL  : %d CpGs (genome-wide, before common-platform restriction)\n", nrow(all)))

# ---- Diagnostics: isolate WHERE the cross-platform intersection collapses -
cat("\n[Diagnostics] Per-platform coverage (>=1 non-NA Blood/Fibro sample):\n")
platform_masks <- list()
for (p in c("EPIC", "ONT", "WGEC", "TWIST", "RRBS")) {
  m <- platformsCoveredMask(all, p)
  platform_masks[[p]] <- m
  cat(sprintf("  %-6s: %d CpGs\n", p, sum(m)))
}

cat("\n[Diagnostics] Pairwise platform-coverage overlap:\n")
pnames <- names(platform_masks)
for (i in seq_along(pnames)) {
  for (j in seq_along(pnames)) {
    if (j <= i) next
    n_ov <- sum(platform_masks[[pnames[i]]] & platform_masks[[pnames[j]]])
    cat(sprintf("  %-6s vs %-6s: %d\n", pnames[i], pnames[j], n_ov))
  }
}

seq4_mask <- platform_masks[["ONT"]] & platform_masks[["WGEC"]] & platform_masks[["TWIST"]] & platform_masks[["RRBS"]]
cat(sprintf("\n[Diagnostics] ONT & WGEC & TWIST & RRBS (sequencing only, no EPIC): %d CpGs\n", sum(seq4_mask)))
cat(sprintf("[Diagnostics] ...of those, also covered by EPIC: %d CpGs\n", sum(seq4_mask & platform_masks[["EPIC"]])))

# Per the paper's stated methodology, this exploratory analysis is
# restricted to the common set of CpGs covered by every platform (at least
# one non-NA sample per platform) -- NOT the full genome-wide matrix. Genome-
# wide input matrices are mostly NA for any single platform pair, let alone
# all five simultaneously, so skipping this step leaves the heatmap/UpSet
# steps downstream with only a handful of literal all-5-complete rows.
common_platform_mask <- platformsCoveredMask(all, c("EPIC", "ONT", "WGEC", "TWIST", "RRBS"))
all <- all[common_platform_mask]
cat(sprintf("\n  ALL  : %d CpGs after restricting to those covered by all 5 platforms\n", nrow(all)))

cat("[2/6] Building per-method matrices...\n")

# Helper: select columns by method prefix and sampleset
select_meth_cols <- function(dt, prefix, sample_ids) {
  cols <- paste0(prefix, "_", sample_ids)
  cols <- intersect(cols, colnames(dt))
  dt[, ..cols]
}

methods <- list(
  EPIC  = cbind(
    select_meth_cols(all, "EPIC",  BLOOD_IDS),
    select_meth_cols(all, "EPIC",  FIBRO_IDS)
  ),
  ONT   = cbind(
    select_meth_cols(all, "ONT",   BLOOD_IDS),
    select_meth_cols(all, "ONT",   FIBRO_IDS)
  ),
  WGEC  = cbind(
    select_meth_cols(all, "WGEC",  BLOOD_IDS),
    select_meth_cols(all, "WGEC",  FIBRO_IDS)
  ),
  TWIST = cbind(
    select_meth_cols(all, "TWIST", BLOOD_IDS),
    select_meth_cols(all, "TWIST", FIBRO_IDS)
  ),
  RRBS  = cbind(
    select_meth_cols(all, "RRBS",  BLOOD_IDS),
    select_meth_cols(all, "RRBS",  FIBRO_IDS)
  )
)

cat("[3/6] Running Wilcoxon tests...\n")

results_list <- list()

for (method in names(methods)) {

  df <- as.data.table(methods[[method]])

  if (is.null(rownames(df))) {
    rownames(df) <- paste0("CpG_", seq_len(nrow(df)))
  }

  blood_cols <- grep(paste0("(?i)^", method, "_blood"),
                     colnames(df), value = TRUE, perl = TRUE)
  fibro_cols <- grep(paste0("(?i)^", method, "_fibro"),
                     colnames(df), value = TRUE, perl = TRUE)

  blood_mat <- as.matrix(df[, ..blood_cols, drop = FALSE])
  fibro_mat <- as.matrix(df[, ..fibro_cols, drop = FALSE])

  mean_blood  <- rowMeans(blood_mat, na.rm = TRUE)
  mean_fibro  <- rowMeans(fibro_mat, na.rm = TRUE)
  delta_beta  <- mean_blood - mean_fibro

  pvals <- vapply(
    seq_len(nrow(blood_mat)),
    FUN = function(i) {
      bt <- as.numeric(blood_mat[i, ])
      ft <- as.numeric(fibro_mat[i, ])
      bt <- bt[!is.na(bt)]
      ft <- ft[!is.na(ft)]
      if (length(bt) < 1 || length(ft) < 1) return(NA_real_)
      tryCatch(
        wilcox.test(bt, ft, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    },
    FUN.VALUE = numeric(1)
  )

  stopifnot(
    length(pvals)      == nrow(df),
    length(delta_beta) == nrow(df)
  )

  fdr <- p.adjust(pvals, method = "BH")

  results_list[[method]] <- data.frame(
    CpG         = rownames(df),
    delta_beta  = delta_beta,
    pvalue      = pvals,
    FDR         = fdr,
    Significant = (!is.na(fdr) & fdr < 0.05 & abs(delta_beta) > DELTA_CUTOFF),
    Method      = method,
    stringsAsFactors = FALSE
  )
}

all_results <- bind_rows(results_list)

fwrite(as.data.table(all_results),
  file.path(opt$datadir, "Wilcoxon_results.csv"),
  sep = ",", quote = FALSE
)
cat(sprintf("  Wilcoxon results: %d rows\n", nrow(all_results)))

cat("[4/6] Generating figures...\n")

# ---- 5.1 Delta-beta density -------------------------------------------------

p_db <- ggplot(all_results, aes(x = delta_beta, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 1.5) +
  scale_color_manual(values = col.vec) +
  labs(title = expression(Delta*beta ~ "distribution per method (Blood - Fibro)")) +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 22),
    axis.title = element_text(size = 22),
    text       = element_text(size = 22)
  )

ggsave(p_db,
  filename = file.path(opt$outdir, "Delta_Beta_Diff.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.2 Heatmaps (top 500 / 1000 / 5000 variable CpGs) --------------------

# ---- Diagnostics: cross-method CpG ID overlap ------------------------------
# The exploratory all_results object is built independently per method (see
# section 3), so a coordinate/ID mismatch between any pair of methods would
# silently collapse the cross-method intersection without an obvious error
# until the complete-case filter below. Surface this explicitly.
cat("\n[Diagnostics] CpGs tested per method (exploratory, before any filtering):\n")
print(table(all_results$Method))

cat("\n[Diagnostics] Sample CpG IDs per method (first 3 each):\n")
for (m in unique(all_results$Method)) {
  ids <- head(all_results$CpG[all_results$Method == m], 3)
  cat(sprintf("  %-6s: %s\n", m, paste(ids, collapse = " | ")))
}

cat("\n[Diagnostics] Pairwise CpG ID overlap between methods:\n")
method_cpg_sets <- split(all_results$CpG, all_results$Method)
method_names_diag <- names(method_cpg_sets)
for (i in seq_along(method_names_diag)) {
  for (j in seq_along(method_names_diag)) {
    if (j <= i) next
    m1 <- method_names_diag[i]; m2 <- method_names_diag[j]
    n_overlap <- length(intersect(method_cpg_sets[[m1]], method_cpg_sets[[m2]]))
    cat(sprintf("  %-6s vs %-6s: %d shared CpGs (of %d / %d)\n",
                m1, m2, n_overlap, length(method_cpg_sets[[m1]]), length(method_cpg_sets[[m2]])))
  }
}
cat("\n")

mat_wide <- all_results %>%
  select(CpG, Method, delta_beta) %>%
  pivot_wider(names_from = Method, values_from = delta_beta) %>%
  column_to_rownames("CpG") %>%
  as.matrix()

# Heatmap clustering (hclust) cannot handle NA/NaN/Inf. With genome-wide
# (not pre-intersected) input matrices, most CpGs are NOT covered by every
# method, so mat_wide is largely sparse/NA. Restrict to CpGs with complete
# data across all methods before ranking by variance -- a cross-method
# variance/heatmap comparison isn't meaningful for a CpG some methods never
# tested anyway.
n_before_complete <- nrow(mat_wide)
mat_wide <- mat_wide[complete.cases(mat_wide), , drop = FALSE]
cat(sprintf(
  "  Heatmap input: %d / %d CpGs have complete data across all methods (%d dropped, NA in at least one method)\n",
  nrow(mat_wide), n_before_complete, n_before_complete - nrow(mat_wide)
))
if (nrow(mat_wide) < 5000) {
  warning(sprintf(
    "Only %d complete-case CpGs available -- Heatmap_top5000.png (and possibly top1000/top500) will use fewer CpGs than its name suggests, or fail if fewer than requested.",
    nrow(mat_wide)
  ))
}

col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
vars    <- apply(mat_wide, 1, var, na.rm = TRUE)

ha <- HeatmapAnnotation(
  Platform = PLATFORMS[colnames(mat_wide)],
  col      = list(Platform = c(
    "Array"     = "#6D8764",
    "ShortRead" = "#76608A",
    "LongRead"  = "#F0A30A"
  )),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 16, fontface = "bold"),
  simple_anno_size     = unit(6, "mm"),
  annotation_legend_param = list(
    title      = "Platform",
    title_gp   = gpar(fontsize = 16, fontface = "bold"),
    labels_gp  = gpar(fontsize = 16),
    grid_height = unit(5, "mm"),
    grid_width  = unit(5, "mm")
  ),
  gp = gpar(col = "black")
)

make_heatmap <- function(n_cpgs, filename) {
  n_available <- nrow(mat_wide)
  if (n_cpgs > n_available) {
    warning(sprintf(
      "%s: requested top %d CpGs but only %d complete-case CpGs are available -- using all %d instead.",
      filename, n_cpgs, n_available, n_available
    ))
    n_cpgs <- n_available
  }
  top_cpgs <- names(sort(vars, decreasing = TRUE))[seq_len(n_cpgs)]
  mat_top  <- mat_wide[top_cpgs, , drop = FALSE]

  ht <- Heatmap(
    as.matrix(mat_top),
    name              = "Delta Beta",
    col               = col_fun,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    show_row_names    = FALSE,
    column_title      = "Delta Beta (Blood - Fibro) across Methods",
    row_title         = sprintf("Top %s variable CpGs",
                                  formatC(n_cpgs, big.mark = ",")),
    heatmap_legend_param = list(
      title           = "Delta Beta (Blood - Fibro)",
      title_gp        = gpar(fontsize = 17, fontface = "bold"),
      labels_gp       = gpar(fontsize = 16),
      legend_width    = unit(7, "cm"),
      legend_direction = "horizontal",
      title_position  = "topcenter"
    ),
    column_title_gp   = gpar(fontsize = 18, fontface = "bold"),
    row_title_gp      = gpar(fontsize = 18, fontface = "bold"),
    column_names_gp   = gpar(fontsize = 18, fontface = "bold"),
    top_annotation    = ha
  )

  png(file.path(opt$outdir, filename),
    width = 14, height = 12, units = "in", res = 300)
  draw(ht,
    merge_legend            = TRUE,
    heatmap_legend_side     = "bottom",
    annotation_legend_side  = "bottom",
    padding                 = unit(c(5, 10, 5, 10), "mm")
  )
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

make_heatmap(500,  "Heatmap_top500.png")
make_heatmap(1000, "Heatmap_top1000.png")
make_heatmap(5000, "Heatmap_top5000.png")

# ---- 5.3 UpSet plot – limma DMCs --------------------------------------------

limma_files <- list(
  EPIC  = file.path(opt$limma_dir, "EPIC_Blood_vs_Fibroblast.csv"),
  ONT   = file.path(opt$limma_dir, "ONT_Blood_vs_Fibroblast.csv"),
  TWIST = file.path(opt$limma_dir, "TWIST_Blood_vs_Fibroblast.csv"),
  RRBS  = file.path(opt$limma_dir, "RRBS_Blood_vs_Fibroblast.csv"),
  WGEC  = file.path(opt$limma_dir, "WGEC_Blood_vs_Fibroblast.csv")
)

missing_limma <- names(limma_files)[!sapply(limma_files, file.exists)]
if (length(missing_limma) > 0) {
  warning(paste("Limma files not found – skipping UpSet limma:",
                paste(missing_limma, collapse = ", ")))
} else {
  limma_data <- lapply(names(limma_files), function(m) {
    dt <- fread(limma_files[[m]], header = FALSE, sep = ",")
    dt[dt$V6 < 0.05, ]$V1
  })
  names(limma_data) <- names(limma_files)

  cat("\n[Diagnostics] Significant (p<0.05) DMCs per method (for UpSet_limma):\n")
  print(sapply(limma_data, length))

  upset_data <- UpSetR::fromList(limma_data)

  cat("\n[Diagnostics] upset_data structure (columns must match intersect = c(...) below):\n")
  cat("  dim: ", paste(dim(upset_data), collapse = " x "), "\n")
  cat("  colnames: ", paste(colnames(upset_data), collapse = ", "), "\n\n")

  # UpSetR::fromList() silently DROPS any method with zero elements (zero
  # significant DMCs at p<0.05) from the resulting columns. Hardcoding all
  # 5 platform names in intersect=/queries= then crashes with a cryptic
  # "undefined columns selected" if any platform's column is missing.
  # Restrict to platforms actually present, and say so explicitly.
  requested_platforms <- c("ONT", "WGEC", "RRBS", "TWIST", "EPIC")
  available_platforms <- intersect(requested_platforms, colnames(upset_data))
  dropped_platforms   <- setdiff(requested_platforms, available_platforms)
  if (length(dropped_platforms) > 0) {
    warning(sprintf(
      "UpSet_limma.png: platform(s) with zero significant (p<0.05) DMCs were dropped from the plot: %s",
      paste(dropped_platforms, collapse = ", ")
    ))
  }
  if (length(available_platforms) < 2) {
    warning("UpSet_limma.png: fewer than 2 platforms have any significant DMCs -- skipping this plot entirely.")
  } else {

  stripe_df <- data.frame(
    set    = names(PLATFORMS)[names(PLATFORMS) %in% available_platforms],
    labeli = PLATFORMS[names(PLATFORMS) %in% available_platforms]
  )

  set_queries       <- lapply(available_platforms, function(p) upset_query(set = p, fill = METHOD_COLORS[[p]]))
  intersect_queries <- lapply(available_platforms, function(p) upset_query(intersect = p, color = METHOD_COLORS[[p]], fill = METHOD_COLORS[[p]]))

  png(file.path(opt$outdir, "UpSet_limma.png"),
    width = 16, height = 10, units = "in", res = 400)

  print(ComplexUpset::upset(
    upset_data,
    intersect = available_platforms,
    queries   = c(set_queries, intersect_queries),
    set_sizes = upset_set_size(geom = geom_bar(width = 0.8)) +
      ylab("DMCs") +
      scale_y_continuous(
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      theme(text = element_text(size = 25)),
    base_annotations = list(
      "Intersection size" = intersection_size(width = 0.8, counts = FALSE) +
        theme(
          plot.background = element_rect(fill = "lightgray"),
          text = element_text(size = 25)
        ) +
        ylab("Overlapping DMCs")
    ),
    stripes  = upset_stripes(
      mapping = aes(color = labeli),
      colors  = PLATFORM_COLORS,
      data    = stripe_df
    ),
    name     = "Method specific overlap",
    min_size = 100
  ) +
    guides(color = guide_legend(title = NULL)) +
    theme(
      text         = element_text(size = 25),
      axis.text    = element_text(size = 18),
      axis.title   = element_text(size = 20),
      strip.text   = element_text(size = 18),
      legend.text  = element_text(size = 18),
      legend.title = element_text(size = 20)
    ))

  dev.off()
  cat("  Saved: UpSet_limma.png\n")
  }  # end: length(available_platforms) >= 2
}

# ---- 5.4 UpSet plot – Wilcoxon DMCs (EPIC vs TWIST) ------------------------

df_sig <- all_results %>%
  filter(FDR < 0.05) %>%
  mutate(Significant = TRUE) %>%
  select(CpG, Significant, Method) %>%
  pivot_wider(
    names_from  = Method,
    values_from = Significant,
    values_fill = FALSE
  )

stripe_df_wil <- data.frame(
  set    = c("TWIST", "EPIC"),
  labeli = c("ShortRead", "Array")
)

png(file.path(opt$outdir, "UpSet_wilcoxon.png"),
  width = 16, height = 10, units = "in", res = 400)

print(ComplexUpset::upset(
  df_sig,
  intersect = c("TWIST", "EPIC"),
  queries   = list(
    upset_query(set = "EPIC",  fill = "#009E73"),
    upset_query(set = "TWIST", fill = "#DC79A7"),
    upset_query(intersect = "EPIC",  color = "#009E73", fill = "#009E73"),
    upset_query(intersect = "TWIST", color = "#DC79A7", fill = "#DC79A7")
  ),
  set_sizes = upset_set_size(geom = geom_bar(width = 0.8)) +
    ylab("DMCs") +
    scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    theme(text = element_text(size = 25)),
  base_annotations = list(
    "Intersection size" = intersection_size(width = 0.8, counts = FALSE) +
      theme(
        plot.background = element_rect(fill = "lightgray"),
        text = element_text(size = 25)
      ) +
      ylab("Overlapping DMCs")
  ),
  stripes  = upset_stripes(
    mapping = aes(color = labeli),
    colors  = PLATFORM_COLORS,
    data    = stripe_df_wil
  ),
  name     = "Method specific overlap",
  min_size = 100
) +
  guides(color = guide_legend(title = NULL)) +
  theme(
    text         = element_text(size = 25),
    axis.text    = element_text(size = 18),
    axis.title   = element_text(size = 20),
    strip.text   = element_text(size = 18),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 20)
  ))

dev.off()
cat("  Saved: UpSet_wilcoxon.png\n")

# ---- 5.5 Cross-platform EPIC vs TWIST delta-beta scatter --------------------

epic_df  <- all_results %>%
  filter(Method == "EPIC") %>%
  select(CpG,
    delta_beta_epic = delta_beta, pvalue_epic = pvalue,
    FDR_epic = FDR, sig_epic = Significant)

twist_df <- all_results %>%
  filter(Method == "TWIST") %>%
  select(CpG,
    delta_beta_twist = delta_beta, pvalue_twist = pvalue,
    FDR_twist = FDR, sig_twist = Significant)

merged_scatter <- inner_join(epic_df, twist_df, by = "CpG") %>%
  mutate(
    concordance = case_when(
      sig_epic & sig_twist & sign(delta_beta_epic) == sign(delta_beta_twist) ~ "Concordant",
      sig_epic  & !sig_twist ~ "EPIC only",
      sig_twist & !sig_epic  ~ "TWIST only",
      TRUE ~ "Not significant"
    )
  )

r_val <- round(cor(merged_scatter$delta_beta_epic,
                    merged_scatter$delta_beta_twist,
                    use = "pairwise.complete.obs"), 2)

p_scatter <- ggplot(
  merged_scatter,
  aes(x = delta_beta_epic, y = delta_beta_twist, color = concordance)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1.2,
              color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "red",
              linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(values = c(
    "Concordant"      = "green",
    "EPIC only"       = "#009E73",
    "TWIST only"      = "#DC79A7",
    "Not significant" = "grey80"
  )) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title    = "Cross-platform comparison of methylation changes",
    subtitle = expression(Delta*beta ~ "concordance between EPIC and TWIST (Blood – Fibro)"),
    x        = expression(Delta*beta ~ "(EPIC)"),
    y        = expression(Delta*beta ~ "(TWIST)"),
    color    = "Concordance"
  ) +
  annotate("text", x = -0.9, y = 0.9,
           label = paste0("r = ", r_val),
           size = 9, hjust = 0, fontface = "bold") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text   = element_text(size = 25),
    axis.title  = element_text(size = 25),
    text        = element_text(size = 25)
  )

ggsave(p_scatter,
  filename = file.path(opt$outdir, "Cross_Platform_EPIC_TWIST.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.6 Variance per method (Blood and Fibro, no GIAB) ---------------------

cat("[5/6] Computing variance and coverage...\n")

value_cols_all <- colnames(all)[
  !colnames(all) %in% c("chr", "start", "coord") &
  !grepl("_cov_", colnames(all))
]

meta <- tibble(col = value_cols_all) %>%
  mutate(
    Method    = str_extract(col, "EPIC|ONT|TWIST|RRBS|WGEC"),
    SampleSet = case_when(
      str_detect(col, "Blood") ~ "Blood",
      str_detect(col, "Fibro") ~ "Fibro",
      str_detect(col, "GIAB")  ~ "GIAB"
    )
  ) %>%
  filter(!is.na(Method), !is.na(SampleSet), Method != "PacBio")

variance_df <- meta %>%
  group_by(Method, SampleSet) %>%
  summarise(
    Variance = list(apply(all[, ..col], 1, var, na.rm = TRUE)),
    .groups  = "drop"
  ) %>%
  unnest(Variance)

fwrite(as.data.table(variance_df),
  file.path(opt$datadir, "Variances.csv"),
  sep = ",", quote = FALSE
)

p_var <- ggplot(
  variance_df[variance_df$SampleSet != "GIAB", ],
  aes(x = Method, y = Variance, fill = Method)
) +
  geom_boxplot() +
  scale_fill_manual(values = col.vec) +
  facet_wrap(~SampleSet) +
  labs(title = "Variance of Blood/Fibro DMC-set across Methods") +
  theme_bw() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    axis.text   = element_text(size = 26),
    axis.title  = element_text(size = 26),
    text        = element_text(size = 26),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(p_var,
  filename = file.path(opt$outdir, "Variance_DMCs.png"),
  height = 12, width = 14, dpi = 300
)

# ---- 5.7 Coverage per method (Blood and Fibro) ------------------------------

merged_cov <- merge(all, blood, by = c("chr", "start"))
merged_cov <- merge(merged_cov, fibro, by = c("chr", "start"))

cov_col_names <- grep("_cov_", colnames(merged_cov), value = TRUE)
covs          <- as.data.table(merged_cov)[, ..cov_col_names]

coverage_df <- covs %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "col",
    values_to = "Coverage"
  ) %>%
  mutate(
    Method    = str_extract(col, "ONT|TWIST|WGEC|RRBS"),
    SampleSet = case_when(
      str_detect(col, "Blood") ~ "Blood",
      str_detect(col, "Fibro") ~ "Fibro"
    )
  ) %>%
  filter(!is.na(Method), !is.na(SampleSet)) %>%
  select(SampleSet, Method, Coverage)

p_cov <- ggplot(coverage_df, aes(x = Method, y = Coverage, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = col.vec) +
  scale_y_log10() +
  facet_wrap(~SampleSet) +
  labs(
    title = "Coverage of Blood/Fibro DMC-set across Methods",
    y     = "Coverage (log10 scale)"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    axis.text   = element_text(size = 26),
    axis.title  = element_text(size = 26),
    text        = element_text(size = 26),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(p_cov,
  filename = file.path(opt$outdir, "Coverage_DMCs.png"),
  height = 12, width = 14, dpi = 300
)

cat(sprintf("\n[6/6] Done. Figures written to: %s\n", opt$outdir))
