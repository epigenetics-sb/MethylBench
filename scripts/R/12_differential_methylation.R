#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - measurement properties on the common five-platform CpG set
# =============================================================================
# Restricts all five platforms to the CpGs covered by every platform and
# describes them without any statistical testing: blood-fibroblast delta-beta
# per platform and its pairwise concordance, per-CpG coverage and within-group
# variance. Differential methylation testing is done in 13_DMR_DSS_analysis.R
# (DSS for sequencing data, limma on M-values for EPIC).
# Figure 7; Suppl. Figure 7.
#
# Usage:
#   Rscript scripts/R/12_differential_methylation.R --all_path ALL.csv \
#     --blood_path Blood_without_EPIC.csv --fibro_path Fibro_without_EPIC.csv \
#     --outdir results/figures/ --datadir results/diff_meth/
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
  library(circlize)
  library(grid)
  library(patchwork)
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
  make_option("--outdir",
    type    = "character", default = "results/figures/",
    help    = "Output directory for figures [default: results/figures/]",
    metavar = "DIR"),
  make_option("--datadir",
    type    = "character", default = "results/diff_meth/",
    help    = "Output directory for data files [default: results/diff_meth/]",
    metavar = "DIR")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("all_path", "blood_path", "fibro_path")) {
  val <- opt[[arg]]
  if (is.null(val)) stop(paste("ERROR: --", arg, " is required", sep = ""))
}

dir.create(opt$outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(opt$datadir, recursive = TRUE, showWarnings = FALSE)

col.vec      <- get_colors()
BLOOD_IDS    <- paste0("Blood", 1:5)
FIBRO_IDS    <- paste0("Fibro",  1:5)

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

METHOD_PREFIX <- c(
  ONT   = "ONT",
  TWIST = "TWIST",
  WGEC  = "WGEC",
  RRBS  = "RRBS"
)

# ---- 1. Load matrices and restrict to the common CpG set --------------------

cat("[1/6] Loading matrices...\n")

all   <- fread(opt$all_path,   header = TRUE, sep = ",", na.strings = "NA")
blood <- fread(opt$blood_path, header = TRUE, sep = ",", na.strings = "NA")
fibro <- fread(opt$fibro_path, header = TRUE, sep = ",", na.strings = "NA")

cat(sprintf("  ALL  : %d CpGs (genome-wide, before common-platform restriction)\n", nrow(all)))

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

# Common set: CpGs with at least one non-NA sample on every platform.
common_platform_mask <- platformsCoveredMask(all, c("EPIC", "ONT", "WGEC", "TWIST", "RRBS"))
all <- all[common_platform_mask]
cat(sprintf("\n  ALL  : %d CpGs after restricting to those covered by all 5 platforms\n", nrow(all)))

cat("[2/6] Building per-method matrices...\n")

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

# ---- 2. Delta-beta per platform (no testing) --------------------------------

cat("[3/6] Computing delta-beta per platform...\n")

results_list <- list()

for (method in names(methods)) {

  df <- as.data.table(methods[[method]])

  blood_cols <- grep(paste0("(?i)^", method, "_blood"),
                     colnames(df), value = TRUE, perl = TRUE)
  fibro_cols <- grep(paste0("(?i)^", method, "_fibro"),
                     colnames(df), value = TRUE, perl = TRUE)

  mean_blood <- rowMeans(as.matrix(df[, ..blood_cols, drop = FALSE]), na.rm = TRUE)
  mean_fibro <- rowMeans(as.matrix(df[, ..fibro_cols, drop = FALSE]), na.rm = TRUE)

  # Rows of `all` are identical for all methods, so the row index is the CpG key.
  results_list[[method]] <- data.frame(
    CpG        = paste0(all$chr, ":", all$start),
    delta_beta = mean_blood - mean_fibro,
    Method     = method,
    stringsAsFactors = FALSE
  )
}

all_results <- bind_rows(results_list)
all_results$delta_beta[is.nan(all_results$delta_beta)] <- NA_real_

fwrite(as.data.table(all_results),
  file.path(opt$datadir, "DeltaBeta_common_set.csv"),
  sep = ",", quote = FALSE
)
cat(sprintf("  Delta-beta table: %d rows\n", nrow(all_results)))

# ---- 3. Figures -------------------------------------------------------------

cat("[4/6] Generating figures...\n")

# ---- Delta-beta density (Suppl. Figure 7B) ----------------------------------

p_db <- ggplot(all_results, aes(x = delta_beta, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 1.5) +
  scale_color_manual(values = col.vec) +
  labs(title = expression(Delta*beta ~ "distribution per method (Blood - Fibro)"),
       x     = expression(Delta*beta ~ ""),
       y     = "Density") +
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

# ---- Heatmaps of the most variable CpGs (Suppl. Figure 7C-E) ----------------
cat("\n[Diagnostics] CpGs per method on the common set:\n")
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

mat_wide_all <- mat_wide
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

# ---- Pairwise delta-beta concordance (Figure 7A, 7B) ------------------------

PLAT_ORDER <- c("EPIC", "TWIST", "WGEC", "RRBS", "ONT")
db_wide    <- mat_wide_all[, intersect(PLAT_ORDER, colnames(mat_wide_all)), drop = FALSE]
pairs      <- combn(colnames(db_wide), 2, simplify = FALSE)

cor_tab <- rbindlist(lapply(pairs, function(p) {
  ok <- complete.cases(db_wide[, p])
  data.table(
    Platform1 = p[1], Platform2 = p[2], n_CpGs = sum(ok),
    pearson_r  = cor(db_wide[ok, p[1]], db_wide[ok, p[2]], method = "pearson"),
    spearman_rho = cor(db_wide[ok, p[1]], db_wide[ok, p[2]], method = "spearman"),
    MAD        = mean(abs(db_wide[ok, p[1]] - db_wide[ok, p[2]])),
    Comparison = ifelse("EPIC" %in% p, "EPIC-sequencing", "sequencing-sequencing")
  )
}))
fwrite(cor_tab, file.path(opt$datadir, "DeltaBeta_pairwise_correlation.csv"))
cat("\n[Pairwise delta-beta concordance, common set]\n")
print(cor_tab)
cat(sprintf("  Pearson r, all pairs: %.2f-%.2f; EPIC-sequencing: %.2f-%.2f; sequencing-sequencing: %.2f-%.2f\n",
  min(cor_tab$pearson_r), max(cor_tab$pearson_r),
  min(cor_tab[Comparison == "EPIC-sequencing", pearson_r]),
  max(cor_tab[Comparison == "EPIC-sequencing", pearson_r]),
  min(cor_tab[Comparison == "sequencing-sequencing", pearson_r]),
  max(cor_tab[Comparison == "sequencing-sequencing", pearson_r])))

# Figure 7A: lower-triangle correlation matrix.
cor_plot <- rbind(
  cor_tab[, .(x = Platform1, y = Platform2, r = pearson_r)],
  data.table(x = colnames(db_wide), y = colnames(db_wide), r = NA_real_)
)
cor_plot[, `:=`(x = factor(x, levels = colnames(db_wide)),
                y = factor(y, levels = rev(colnames(db_wide))))]

p_cor <- ggplot(cor_plot, aes(x = x, y = y)) +
  geom_tile(data = cor_plot[!is.na(r)], aes(fill = r), color = "white", linewidth = 1) +
  geom_text(data = cor_plot[!is.na(r)], aes(label = sprintf("%.2f", r)), size = 9) +
  geom_text(data = cor_plot[is.na(r)], aes(label = x), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "#F1EFE8", high = "#0F6E56",
                      limits = c(floor(min(cor_tab$pearson_r) * 20) / 20, 1),
                      name = expression("Pearson r of" ~ Delta*beta),
                      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_fixed() +
  labs(title = expression(Delta*beta ~ "concordance (Blood - Fibro), common CpG set"),
       subtitle = sprintf("n = %s CpGs; no statistical testing",
                          formatC(nrow(all), big.mark = ","))) +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), text = element_text(size = 22),
        legend.position = "bottom", legend.key.width = unit(2, "cm"))

ggsave(p_cor, filename = file.path(opt$outdir, "DeltaBeta_correlation_matrix.png"),
       height = 12, width = 12, dpi = 300)

# Pairwise hexbin scatter in the style of Figure 8A.
delta_hex <- function(m1, m2, text_size = 12) {
  ok <- complete.cases(db_wide[, c(m1, m2)])
  d  <- data.frame(x = db_wide[ok, m1], y = db_wide[ok, m2])
  r  <- cor(d$x, d$y)
  ggplot(d, aes(x, y)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "CpGs (log10)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = -0.9, y = 0.9, label = sprintf("r = %.2f", r),
             hjust = 0, size = text_size / 3, fontface = "bold", color = "black") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(x = bquote(Delta*beta ~ .(m1)), y = bquote(Delta*beta ~ .(m2))) +
    theme_bw() +
    theme(text = element_text(size = text_size), axis.text = element_text(size = text_size))
}

# Figure 7B: EPIC vs. TWIST, without significance colouring.
p_scatter <- delta_hex("EPIC", "TWIST", text_size = 26) +
  labs(title = expression(Delta*beta ~ "concordance between EPIC and TWIST (Blood vs. Fibro)"))
ggsave(p_scatter,
  filename = file.path(opt$outdir, "Cross_Platform_EPIC_TWIST.png"),
  height = 12, width = 14, dpi = 300
)

# ---- Suppl. Figure 7A: pairwise delta-beta concordance, all ten pairs --------
# One facet_wrap, so all panels share one fill scale.
db_long <- dcast(as.data.table(all_results), CpG ~ Method, value.var = "delta_beta")
plats   <- intersect(PLAT_ORDER, names(db_long))
pair_df <- rbindlist(lapply(combn(plats, 2, simplify = FALSE), function(p)
  data.table(Pair = paste(p[1], "vs.", p[2]), x = db_long[[p[1]]], y = db_long[[p[2]]])))
pair_df <- pair_df[!is.na(x) & !is.na(y)]
pair_df[, Pair := factor(Pair, levels = unique(Pair))]
r_df <- pair_df[, .(r = cor(x, y)), by = Pair]

p_pairs <- ggplot(pair_df, aes(x, y)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "CpGs (log10)",
                       limits = c(1, NA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_text(data = r_df, aes(x = -0.95, y = 0.9, label = sprintf("r = %.2f", r)),
            hjust = 0, size = 5, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ Pair, ncol = 5) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(title = expression("Pairwise" ~ Delta*beta ~ "concordance on the common five-platform CpG set"),
       subtitle = sprintf("n = %s CpGs | red dashed = identity line | r = Pearson | no statistical testing",
                          formatC(nrow(all), big.mark = ",")),
       x = expression(Delta*beta ~ "(first platform)"),
       y = expression(Delta*beta ~ "(second platform)")) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave(p_pairs, filename = file.path(opt$outdir, "SuppFig7A_delta_beta_pairwise.png"),
       width = 20, height = 9, dpi = 300, bg = "white")

# ---- Within-group variance per method (Figure 7D) ---------------------------

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
  labs(title = "Variance of Blood/Fibro CpG-set across Methods") +
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

# ---- Per-CpG coverage on the common CpG set (Figure 7C) ---------------------

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
    title = "Coverage of Blood/Fibro CpG-set across Methods",
    y     = "Coverage (log10 scale)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
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
