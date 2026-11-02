#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - cpg_stats_methylbench.tab
# =============================================================================
# Per sample: merges the per-CpG calls of all platforms on genomic coordinates
# and counts the CpGs shared by all platforms at increasing coverage thresholds.
# Output columns: Sample, Sampleset, Cov_0, Cov_5, ..., Cov_40.
#
# Usage:
#   Rscript scripts/R/07_generate_cpg_stats.R --samplesheet samplesheet.csv \
#     (--calls_sheet calls.csv | --datadir data/) --outdir <stats>/ [--coverages 0,5,10,...]
#
# Input files: --calls_sheet (columns Sample, Method, File; as for
# 07b_build_matrices.R) or, under --datadir, <METHOD>/<sample><suffix>:
#   ONT/<sample>_modkit_pileup.bed
#   RRBS/<sample>.tab (lab-specific table) or <sample>.bismark.cov.gz
#   WGEC/, TWIST/<sample>.bismark.cov.gz
#   PacBio/<sample>.GRCh38.pbmm2.combined.bed   (GIAB only)
# All files are read with readPlatformCalls() (utils/helpers.R), as in 07b and
# 17; positions are 0-based.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--samplesheet",
    type    = "character",
    help    = "Sample sheet, comma- or tab-separated (columns: Sample, Sampleset) [required]",
    metavar = "FILE"
  ),
  make_option("--calls_sheet",
    type    = "character",
    default = NULL,
    help    = "Table with columns Sample, Method, File (alternative to --datadir), as for 07b_build_matrices.R",
    metavar = "FILE"
  ),
  make_option("--datadir",
    type    = "character",
    default = "data/",
    help    = "Root data directory containing per-method subdirectories [default: data/]",
    metavar = "DIR"
  ),
  make_option("--outdir",
    type    = "character",
    default = "results/qc/",
    help    = "Output directory [default: results/qc/]",
    metavar = "DIR"
  ),
  make_option("--coverages",
    type    = "character",
    default = "0,5,10,15,20,25,30,35,40",
    help    = "Comma-separated coverage thresholds [default: 0,5,10,15,20,25,30,35,40]",
    metavar = "STRING"
  ),
  make_option("--ont_suffix",
    type    = "character",
    default = "_modkit_pileup.bed",
    help    = "ONT file suffix [default: _modkit_pileup.bed]"
  ),
  make_option("--bismark_suffix",
    type    = "character",
    default = ".bismark.cov.gz",
    help    = "Bismark file suffix for RRBS/WGEC/TWIST [default: .bismark.cov.gz]"
  ),
  make_option("--rrbs_suffix",
              type    = "character",
              default = ".tab",
              help    = "RRBS file suffix (lab-specific .tab format) [default: .tab]"
  ),
  make_option("--pacbio_suffix",
    type    = "character",
    default = ".GRCh38.pbmm2.combined.bed",
    help    = "PacBio file suffix (GIAB only) [default: .GRCh38.pbmm2.combined.bed]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$samplesheet)) stop("ERROR: --samplesheet is required")
if (!file.exists(opt$samplesheet)) stop(paste("ERROR: File not found:", opt$samplesheet))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

COV_THRESHOLDS <- as.integer(strsplit(opt$coverages, ",")[[1]])

# Per-platform calls with the shared reader (utils/helpers.R): one row per CpG,
# 0-based position of the C on the + strand, both strands combined, cov > 0.
suffixes <- list(ONT = opt$ont_suffix, BISMARK = opt$bismark_suffix,
                 RRBS = opt$rrbs_suffix, PacBio = opt$pacbio_suffix)
calls_sheet <- if (is.null(opt$calls_sheet)) NULL else readCallsSheet(opt$calls_sheet)

load_method <- function(method, sample_id, datadir, opt) {
  d <- readPlatformCalls(datadir, method, sample_id, suffixes,
                         file = callsSheetFile(calls_sheet, method, sample_id))
  if (is.null(d)) {
    warning(sprintf("%s file not found for %s", method, sample_id))
    return(NULL)
  }
  d[, .(coord = paste0(chr, ":", start), coverage = cov, methylation = beta)]
}

# Number of CpGs at which all platforms reach each coverage threshold.
count_cpgs_at_thresholds <- function(merged_dt, methods_present, thresholds) {

  cov_cols <- paste0("cov_", methods_present)

  counts <- sapply(thresholds, function(thr) {
    if (thr == 0) {
      return(nrow(merged_dt))
    }
    mask <- Reduce(`&`, lapply(cov_cols, function(col) {
      merged_dt[[col]] >= thr
    }))
    return(sum(mask, na.rm = TRUE))
  })

  names(counts) <- paste0("Cov_", thresholds)
  return(counts)
}

ss <- as.data.table(readSampleSheet(opt$samplesheet))

required_cols <- c("Sample", "Sampleset")
missing_cols  <- setdiff(required_cols, colnames(ss))
if (length(missing_cols) > 0) {
  stop(paste("Samplesheet missing columns:", paste(missing_cols, collapse = ", ")))
}

results <- vector("list", nrow(ss))

for (i in seq_len(nrow(ss))) {

  sample_id  <- ss$Sample[i]
  sampleset  <- ss$Sampleset[i]
  is_giab    <- sampleset == "GIAB"

  methods <- if (is_giab) {
    c("ONT", "RRBS", "WGEC", "TWIST", "PacBio")
  } else {
    c("ONT", "RRBS", "WGEC", "TWIST")
  }

  method_dts <- lapply(methods, function(m) {
    dt <- load_method(m, sample_id, opt$datadir, opt)
    if (is.null(dt)) return(NULL)
    setnames(dt,
      old = c("coverage", "methylation"),
      new = c(paste0("cov_", m), paste0("meth_", m))
    )
    return(dt)
  })
  names(method_dts) <- methods

  method_dts <- Filter(Negate(is.null), method_dts)

  if (length(method_dts) == 0) {
    warning(sprintf("  No data loaded for sample %s.", sample_id))
    next
  }

  merged <- Reduce(
    function(a, b) merge(a, b, by = "coord", all = FALSE),
    method_dts
  )

  methods_present <- names(method_dts)
  counts <- count_cpgs_at_thresholds(merged, methods_present, COV_THRESHOLDS)

  row <- data.table(
    Sample    = sample_id,
    Sampleset = sampleset
  )
  for (col_name in names(counts)) {
    row[, (col_name) := counts[[col_name]]]
  }

  results[[i]] <- row
}

cpg_stats <- rbindlist(Filter(Negate(is.null), results), fill = TRUE)

col_order <- c("Sample", "Sampleset", paste0("Cov_", COV_THRESHOLDS))
col_order <- intersect(col_order, colnames(cpg_stats))
cpg_stats <- cpg_stats[, ..col_order]

cpg_stats[, sort_order := fcase(
  Sampleset == "Blood",      1L,
  Sampleset == "Fibroblast", 2L,
  Sampleset == "GIAB",       3L,
  default = 4L
)]
setorder(cpg_stats, sort_order, Sample)
cpg_stats[, sort_order := NULL]

out_file <- file.path(opt$outdir, "cpg_stats_methylbench.tab")
fwrite(cpg_stats, out_file, sep = "\t", quote = FALSE, na = "NA")
