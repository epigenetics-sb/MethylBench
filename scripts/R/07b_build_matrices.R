#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - merged methylation matrices (input of scripts 08-16 and 18)
# =============================================================================
# Builds the five CpG x (platform, sample) matrices described in the README
# ("Input data") from the per-platform methylation calls:
#   Blood_without_EPIC.csv, Fibro_without_EPIC.csv, GIAB_without_EPIC.csv,
#   ALL_without_EPIC.csv (all sequencing platforms, all samples) and
#   ALL.csv (ALL_without_EPIC.csv joined with the EPIC beta-values).
#
# Steps
#   1. Per sample and platform, the calls are read with readPlatformCalls()
#      (utils/helpers.R; the same reader as in script 17): one row per CpG,
#      position of the C on the + strand, both strands combined
#      (WGEC, TWIST: Bismark coverage2cytosine --merge_CpG or merging of
#      per-strand rows; RRBS table: merging of per-strand rows;
#      ONT: modkit --combine-strands, 5mC; PacBio: pb-CpG-tools, type "Total",
#      beta = mod_score / 100), CpGs with zero coverage removed, restricted to
#      chr1-22, chrX and chrY (--chromosomes canonical).
#   2. Platforms are joined within each sample  (--platform_join),
#      samples within each sample set            (--sample_join),
#      and the three sample sets into ALL         (--set_join);
#      each join is "inner" (CpGs present in all parts) or "outer" (union).
#   3. The EPIC matrices are joined with each other and with ALL (--epic_join;
#      "inner": CpGs with an EPIC probe in every EPIC matrix, i.e. the common
#      five-platform CpG set; "left": all CpGs of ALL_without_EPIC.csv).
# The defaults (inner joins throughout, canonical chromosomes, 0-based
# positions) reproduce the matrices used for the manuscript.
#
# Positions: column 'start' is the 0-based position of the C on the + strand
#   (BED convention); --one_based writes 1-based positions instead. EPIC
#   positions are joined as given, i.e. they must use the same convention.
# A summary of the number of CpGs per platform, sample and file is written to
# <outdir>/matrix_build_summary.tsv.
#
# Usage:
#   Rscript scripts/R/07b_build_matrices.R \
#     --samplesheet samplesheet.csv \
#     (--calls_dir /path/to/calls | --calls_sheet calls.csv) \
#     --epic_path EPIC_blood_fibro.csv,EPIC_GIAB.csv \
#     --outdir "$DATA"
#
# Input of the per-platform calls, either
#   --calls_dir   <calls_dir>/<METHOD>/<Sample><suffix>, as for script 17
#                 (README, "Inputs for script 17"), or
#   --calls_sheet comma- or tab-separated table with columns Sample, Method,
#                 File (one row per sample and platform; METHOD in ONT, WGEC,
#                 TWIST, RRBS, PacBio). RRBS files ending in --rrbs_suffix are
#                 read as the lab-specific table, all others as Bismark files.
# --epic_path: EPIC beta matrices (tab or comma) with a 'coord' column (chr:pos)
#   or 'chr' and 'start' columns and one column per sample named EPIC_<Sample>
#   (normalized beta-values, e.g. from RnBeads); probes at the same position
#   are averaged.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

source("scripts/R/utils/helpers.R")

option_list <- list(
  make_option("--samplesheet", type = "character",
              help = "Sample sheet (comma-separated; columns Sample, Sampleset) [required]"),
  make_option("--calls_dir", type = "character", default = NULL,
              help = "Directory with per-platform calls, <calls_dir>/<METHOD>/<Sample><suffix>"),
  make_option("--calls_sheet", type = "character", default = NULL,
              help = "Table with columns Sample, Method, File (alternative to --calls_dir)"),
  make_option("--epic_path", type = "character", default = NULL,
              help = "Comma-separated EPIC beta matrices; without it, ALL.csv is not written"),
  make_option("--outdir", type = "character", default = "data/matrices/"),
  make_option("--platform_join", type = "character", default = "inner",
              help = "Join of the platforms within a sample: inner or outer [default %default]"),
  make_option("--sample_join", type = "character", default = "inner",
              help = "Join of the samples within a sample set: inner or outer [default %default]"),
  make_option("--set_join", type = "character", default = "inner",
              help = "Join of the sample sets into ALL: inner or outer [default %default]"),
  make_option("--epic_join", type = "character", default = "inner",
              help = "Join of the EPIC matrices with each other and with ALL: inner or left [default %default]"),
  make_option("--chromosomes", type = "character", default = "canonical",
              help = "canonical (chr1-22, chrX, chrY) or all [default %default]"),
  make_option("--pacbio_value", type = "character", default = "mod_score",
              help = "PacBio beta-value: mod_score or discretized [default %default]"),
  make_option("--ont_suffix",     type = "character", default = CALL_SUFFIXES$ONT),
  make_option("--bismark_suffix", type = "character", default = CALL_SUFFIXES$BISMARK),
  make_option("--rrbs_suffix",    type = "character", default = CALL_SUFFIXES$RRBS),
  make_option("--pacbio_suffix",  type = "character", default = CALL_SUFFIXES$PacBio),
  make_option("--one_based", action = "store_true", default = FALSE,
              help = "Write 1-based positions; default: 0-based position of the C on the + strand"),
  make_option("--digits", type = "integer", default = 6,
              help = "Decimal places of the beta-values in the output [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$samplesheet)) stop("--samplesheet is required")
if (is.null(opt$calls_dir) == is.null(opt$calls_sheet)) stop("give exactly one of --calls_dir and --calls_sheet")
for (a in c("platform_join", "sample_join", "set_join")) {
  if (!opt[[a]] %in% c("inner", "outer")) stop("--", a, " must be 'inner' or 'outer'")
}
if (!opt$epic_join %in% c("inner", "left")) stop("--epic_join must be 'inner' or 'left'")
if (!opt$chromosomes %in% c("canonical", "all")) stop("--chromosomes must be 'canonical' or 'all'")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# readPlatformCalls() returns 0-based positions (as stored in the matrices);
# --one_based adds 1.
OFFSET <- if (opt$one_based) 1L else 0L
CANONICAL <- c(paste0("chr", 1:22), "chrX", "chrY")

suffixes <- list(ONT = opt$ont_suffix, BISMARK = opt$bismark_suffix,
                 RRBS = opt$rrbs_suffix, PacBio = opt$pacbio_suffix)

SETS <- list(
  Blood      = list(file = "Blood_without_EPIC.csv", methods = c("ONT", "WGEC", "TWIST", "RRBS")),
  Fibroblast = list(file = "Fibro_without_EPIC.csv", methods = c("ONT", "WGEC", "TWIST", "RRBS")),
  GIAB       = list(file = "GIAB_without_EPIC.csv",  methods = c("ONT", "WGEC", "TWIST", "RRBS", "PacBio"))
)

ss <- as.data.table(readSampleSheet(opt$samplesheet))
if (!all(c("Sample", "Sampleset") %in% names(ss))) stop("Sample sheet needs columns Sample and Sampleset")
if (!all(ss$Sampleset %in% names(SETS))) {
  stop("Unknown Sampleset(s): ", paste(setdiff(ss$Sampleset, names(SETS)), collapse = ", "),
       " (expected Blood, Fibroblast, GIAB)")
}

calls_sheet <- if (is.null(opt$calls_sheet)) NULL else readCallsSheet(opt$calls_sheet)
call_file <- function(m, smp) callsSheetFile(calls_sheet, m, smp)

summary_rows <- list()
add_summary <- function(level, name, n) {
  summary_rows[[length(summary_rows) + 1]] <<- data.table(level = level, name = name, n_CpGs = n)
}

# Natural chromosome order (chr1..chr22, chrX, chrY, chrM, others).
order_cpgs <- function(dt) {
  lev <- c(CANONICAL, "chrM")
  dt[, chr_rank := match(chr, lev, nomatch = length(lev) + 1L)]
  setorder(dt, chr_rank, chr, start)
  dt[, chr_rank := NULL]
  dt
}

write_matrix <- function(dt, file) {
  out <- copy(dt)
  beta_cols <- setdiff(grep("^(ONT|WGEC|TWIST|RRBS|PacBio|EPIC)_", names(out), value = TRUE),
                       grep("_cov_", names(out), value = TRUE))
  for (cc in beta_cols) set(out, j = cc, value = round(out[[cc]], opt$digits))
  fwrite(order_cpgs(out), file.path(opt$outdir, file), sep = ",", na = "NA", quote = FALSE)
  cat(sprintf("  wrote %-26s %10d CpGs x %3d columns\n", file, nrow(out), ncol(out)))
  add_summary("file", file, nrow(out))
}

# ---- 1./2. Per sample, per sample set ----------------------------------------
set_dts <- list()
for (set_name in names(SETS)) {
  samples <- ss[Sampleset == set_name, Sample]
  if (length(samples) == 0) next
  methods <- SETS[[set_name]]$methods
  cat(sprintf("\n[%s] samples: %s\n", set_name, paste(samples, collapse = ", ")))

  sample_dts <- lapply(samples, function(smp) {
    platform_dts <- lapply(methods, function(m) {
      d <- readPlatformCalls(opt$calls_dir, m, smp, suffixes,
                             file = call_file(m, smp), pacbio_value = opt$pacbio_value)
      if (is.null(d)) {
        warning(sprintf("No %s calls found for %s", m, smp))
        return(NULL)
      }
      src <- basename(attr(d, "src"))
      if (opt$chromosomes == "canonical") d <- d[chr %in% CANONICAL]
      cat(sprintf("  %-7s %-8s %10d CpGs  [%s]\n", m, smp, nrow(d), src))
      add_summary("calls", paste(m, smp, sep = "_"), nrow(d))
      d <- d[, .(chr, start = start + OFFSET, beta, cov = as.integer(round(cov)))]
      setnames(d, c("beta", "cov"), c(paste0(m, "_", smp), paste0(m, "_cov_", smp)))
      d
    })
    joined <- joinOnCpG(platform_dts, opt$platform_join)
    if (!is.null(joined)) {
      cat(sprintf("  %-16s %10d CpGs after %s join of the platforms\n", smp, nrow(joined), opt$platform_join))
      add_summary("sample", smp, nrow(joined))
    }
    joined
  })

  set_dt <- joinOnCpG(sample_dts, opt$sample_join)
  cat(sprintf("  %-16s %10d CpGs after %s join of the samples\n", set_name, nrow(set_dt), opt$sample_join))
  write_matrix(set_dt, SETS[[set_name]]$file)
  set_dts[[set_name]] <- set_dt
}

# ---- 2. All sequencing platforms, all samples --------------------------------
cat("\n[ALL]\n")
all_dt <- joinOnCpG(set_dts, opt$set_join)
write_matrix(all_dt, "ALL_without_EPIC.csv")

# ---- 3. EPIC ------------------------------------------------------------------
if (!is.null(opt$epic_path)) {
  epic_list <- lapply(strsplit(opt$epic_path, ",")[[1]], function(f) {
    e <- readEPICMatrix(f)
    cat(sprintf("  EPIC %-30s %8d probe positions, %d samples\n", basename(f), nrow(e), ncol(e) - 2))
    e
  })
  epic <- joinOnCpG(epic_list, if (opt$epic_join == "inner") "inner" else "outer")
  keep <- c("chr", "start", intersect(paste0("EPIC_", ss$Sample), names(epic)))
  epic <- epic[, ..keep]
  all_epic <- merge(all_dt, epic, by = c("chr", "start"), all.x = (opt$epic_join == "left"))
  cat(sprintf("  ALL.csv: %s join with %d EPIC positions -> %d CpGs\n", opt$epic_join, nrow(epic), nrow(all_epic)))
  write_matrix(all_epic, "ALL.csv")
} else {
  cat("  --epic_path not given: ALL.csv not written\n")
}

summary_dt <- rbindlist(summary_rows)
summary_dt[, `:=`(platform_join = opt$platform_join, sample_join = opt$sample_join,
                  set_join = opt$set_join, epic_join = opt$epic_join, chromosomes = opt$chromosomes)]
fwrite(summary_dt, file.path(opt$outdir, "matrix_build_summary.tsv"), sep = "\t")
cat(sprintf("\nSummary: %s\n", file.path(opt$outdir, "matrix_build_summary.tsv")))
