#!/usr/bin/env Rscript
# MethylBench – Generate cpg_stats_methylbench.tab

# Description:
#   For each sample, reads raw per-CpG methylation files from all platforms,
#   corrects 0-based coordinates to 1-based (+1 to start), merges all methods
#   on genomic coordinates (chr:start), and counts the number of overlapping
#   CpGs at increasing coverage thresholds.
#
#   Output mirrors the format of cpg_stats_methylbench.tab:
#     Sample | Sampleset | Cov_0 | Cov_5 | Cov_10 | Cov_15 | Cov_20 |
#     Cov_25 | Cov_30 | Cov_35 | Cov_40
#
# Usage:
#   Rscript 07_data_summary.R \
#     --samplesheet  config/samples.tsv \
#     --datadir      data/ \
#     --outdir       results/qc/ \
#     --coverages    "0,5,10,15,20,25,30,35,40"
#
# Expected directory structure under --datadir:
#   data/
#   ├── ONT/   <sample_id>_modkit_pileup.bed
#   ├── RRBS/  <sample_id>.bismark.cov.gz
#   ├── WGEC/  <sample_id>.bismark.cov.gz
#   ├── TWIST/ <sample_id>.bismark.cov.gz
#   └── PacBio/<sample_id>.GRCh38.pbmm2.combined.bed   (GIAB only)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# Load shared helpers
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "scripts/R"
)
source(file.path(script_dir, "helpers.R"))

### Argument Parsing

option_list <- list(
  make_option("--samplesheet",
    type    = "character",
    help    = "Path to samplesheet TSV (columns: Sample, Sampleset) [required]",
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
  # Optional: override file path patterns if naming differs
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

### Resolve file path for a sample + method

get_path <- function(datadir, method, sample_id, suffix) {
  file.path(datadir, method, paste0(sample_id, suffix))
}

### Load and standardize one method's data for one sample

load_method <- function(method, sample_id, datadir, opt) {

  if (method == "ONT") {
    path <- get_path(datadir, "ONT", sample_id, opt$ont_suffix)
    if (!file.exists(path)) {
      warning(sprintf("ONT file not found for %s: %s", sample_id, path))
      return(NULL)
    }
    dt <- readBedMethyl(path)
    dt <- dt[modbase == "m"]
    dt <- dt[, .(
      coord      = paste0(chr, ":", start),
      coverage   = Nvalid_cov,
      methylation = fraction_mod
    )]

  } else if (method %in% c("RRBS", "WGEC", "TWIST")) {
    path <- get_path(datadir, method, sample_id, opt$bismark_suffix)
    if (!file.exists(path)) {
      warning(sprintf("%s file not found for %s: %s", method, sample_id, path))
      return(NULL)
    }
    dt <- readBismarkMeth(path, correct_coords = TRUE)
    dt[, start := start + 1L]
    dt <- dt[, .(
      coord       = paste0(chr, ":", start),
      coverage    = coverage,
      methylation = percentage / 100
    )]

  } else if (method == "PacBio") {
    path <- get_path(datadir, "PacBio", sample_id, opt$pacbio_suffix)
    if (!file.exists(path)) {
      warning(sprintf("PacBio file not found for %s: %s", sample_id, path))
      return(NULL)
    }
    dt <- readPacBio(path)
    dt <- dt[, .(
      coord       = paste0(chr, ":", start),
      coverage    = coverage,
      methylation = percentage / 100
    )]

  } else {
    stop(paste("Unknown method:", method))
  }

  # Remove duplicate coordinates
  dt <- dt[!duplicated(coord)]
  return(dt)
}

### Count overlapping CpGs at each coverage threshold

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

### Process each sample

ss <- fread(opt$samplesheet, header = TRUE, sep = "\t")

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
    # Rename columns before merging
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

### Assemble and Write Output

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