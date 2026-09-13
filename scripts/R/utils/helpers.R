#!/usr/bin/env Rscript
# =============================================================================
# MethylBench – Utility & Helper Functions
# =============================================================================
# Description:
#   Shared helper functions used across all MethylBench analysis scripts.
#   Covers:
#     - File readers for all methylation platforms (ONT/modkit, PacBio,
#       Bismark, EPIC array, samplesheets)
#     - Vectorized long-format data construction for density/ridge plots
#       (replaces repetitive rbind-based create_max_* functions)
#     - Cross-platform Pearson correlation computation
#     - Shared color palettes
#
# NOTE on the original create_max_* functions:
#   The original implementations contained a systematic bug where the TWIST
#   coverage filter always used sample index 1 (e.g. TWIST_cov_Blood1) instead
#   of the correct per-sample index (TWIST_cov_Blood2, Blood3, ...) for all
#   non-first samples. This is fixed here by constructing the filter mask
#   programmatically per sample.
#
# Author:  MethylBench – Laufer et al.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

readBedMethyl <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = FALSE, sep = "\t")

  expected_cols <- 18L
  if (ncol(dataset) != expected_cols) {
    warning(sprintf(
      "readBedMethyl: expected %d columns, got %d in file: %s",
      expected_cols, ncol(dataset), path
    ))
  }

  colnames(dataset) <- c(
    "chr", "start", "end", "modbase", "score", "strand",
    "start1", "end1", "color", "Nvalid_cov", "fraction_mod",
    "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
    "Nfail", "Ndiff", "Nnocall"
  )[seq_len(ncol(dataset))]

  return(dataset)
}

readPacBio <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = FALSE, sep = "\t")

  expected_cols <- 9L
  if (ncol(dataset) != expected_cols) {
    warning(sprintf(
      "readPacBio: expected %d columns, got %d in file: %s",
      expected_cols, ncol(dataset), path
    ))
  }

  colnames(dataset) <- c(
    "chr", "start", "end", "score", "haplotype",
    "coverage", "N_modified", "N_unmodified", "percentage"
  )[seq_len(ncol(dataset))]

  return(dataset)
}

#' @return data.table with named columns.
readBismarkMeth <- function(path, correct_coords = FALSE) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = FALSE, sep = " ")

  colnames(dataset) <- c("chr", "start", "coverage", "methylated", "percentage")

  if (isTRUE(correct_coords)) {
    dataset[, start := start + 1L]
  }

  return(dataset)
}

readEPIC <- function(path, sample.name) {
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.character(sample.name), length(sample.name) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = TRUE, sep = "\t")

  if (!"ID" %in% colnames(dataset)) {
    stop("readEPIC: 'ID' column not found in file: ", path)
  }
  if (!sample.name %in% colnames(dataset)) {
    stop(sprintf("readEPIC: sample '%s' not found in file: %s", sample.name, path))
  }

  return(dataset[, c("ID", sample.name), with = FALSE])
}

readSampleSheet <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  ss <- read.table(path, header = TRUE, sep = ",")
  return(ss)
}

buildMethLong <- function(data,
                          samples,
                          methods,
                          coverages      = c(0, 10),
                          filter_methods = NULL) {

  stopifnot(is.data.table(data))
  stopifnot(is.character(samples),    length(samples) >= 1)
  stopifnot(is.character(methods),    length(methods) >= 1)
  stopifnot(is.numeric(coverages),    length(coverages) >= 1)

  if (is.null(filter_methods)) filter_methods <- names(methods)

  rows <- vector("list", length(coverages) * length(samples) * length(methods))
  idx  <- 1L

  for (cov_threshold in coverages) {
    for (smp in samples) {
      if (cov_threshold == 0) {
        mask <- rep(TRUE, nrow(data))
      } else {
        mask <- Reduce(
          `&`,
          lapply(filter_methods, function(m) {
            cov_col <- paste0(methods[m], "_cov_", smp)
            if (!cov_col %in% colnames(data)) {
              stop(sprintf(
                "buildMethLong: coverage column '%s' not found. ",
                cov_col,
                "Check 'methods' prefixes and 'samples' names."
              ))
            }
            data[[cov_col]] >= cov_threshold
          })
        )
      }

      for (method_label in names(methods)) {
        meth_col <- paste0(methods[method_label], "_", smp)

        if (!meth_col %in% colnames(data)) {
          warning(sprintf(
            "buildMethLong: methylation column '%s' not found – skipping.",
            meth_col
          ))
          next
        }

        meth_values <- data[mask, ][[meth_col]]

        rows[[idx]] <- data.frame(
          Methylation = meth_values,
          Coverage    = cov_threshold,
          Method      = method_label,
          Sample      = smp,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  result <- do.call(rbind, rows[seq_len(idx - 1L)])
  result$Method   <- factor(result$Method,   levels = names(methods))
  result$Coverage <- factor(result$Coverage, levels = sort(unique(coverages)))

  return(result)
}

computePairwiseCorr <- function(data,
                                methods,
                                sample,
                                coverage = 0,
                                existing = NULL) {

  stopifnot(is.data.table(data))
  stopifnot(is.character(methods), !is.null(names(methods)))
  stopifnot(is.character(sample), length(sample) == 1)

  method_labels <- names(methods)
  pairs <- combn(method_labels, 2, simplify = FALSE)

  rows <- lapply(pairs, function(pair) {
    m1_label <- pair[1]
    m2_label <- pair[2]
    col1     <- paste0(methods[m1_label], "_", sample)
    col2     <- paste0(methods[m2_label], "_", sample)

    if (!col1 %in% colnames(data)) {
      warning(sprintf("computePairwiseCorr: column '%s' not found – returning NA", col1))
      corr_val <- NA_real_
    } else if (!col2 %in% colnames(data)) {
      warning(sprintf("computePairwiseCorr: column '%s' not found – returning NA", col2))
      corr_val <- NA_real_
    } else {
      corr_val <- cor(data[[col1]], data[[col2]],
                      use    = "pairwise.complete.obs",
                      method = "pearson")
    }

    data.frame(
      Comparison  = paste(m1_label, "vs", m2_label),
      Method1     = m1_label,
      Method2     = m2_label,
      Coverage    = coverage,
      Correlation = corr_val,
      Sample      = sample,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)

  if (!is.null(existing)) {
    stopifnot(is.data.frame(existing))
    result <- rbind(existing, result)
  }

  return(result)
}

C25 <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
  "lightgrey",   "gold1",   "skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6",     "#FDBF6F", "gray70", "khaki2", "maroon",
  "orchid1",     "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1",      "yellow4", "yellow3", "darkorange4", "brown"
)

C12 <- c(
  "dodgerblue2", "#E31A1C", "green4",      "#6A3D9A",
  "#FF7F00",     "grey",    "gold1",        "khaki2",
  "brown",       "darkturquoise", "palegreen2", "orchid1"
)

METHOD_COLORS <- c(
  "ONT"    = "#D69F00",
  "PacBio" = "#E55E00",
  "EPIC"   = "#009E73",
  "TWIST"  = "#DC79A7",
  "WGBS"   = "purple",
  "RRBS"   = "#0072B2",
  "WGEC"   = "purple"
)

get_colors <- function() {
  return(METHOD_COLORS)
}

load_plotting_environment <- function() {
  packages_to_load <- c(
    "ggplot2", "reshape2", "data.table",
    "UpSetR", "ComplexUpset"
  )
  invisible(lapply(packages_to_load, require, character.only = TRUE))
}

load_environment <- function() {
  packages_to_load <- c(
    "ggplot2", "reshape2", "data.table",
    "tidyr",   "dplyr",    "ggridges"
  )
  invisible(lapply(packages_to_load, require, character.only = TRUE))
}

load_environment_diff_meth <- function() {
  packages_to_load <- c(
    "ggplot2", "reshape2",      "data.table",
    "tidyr",   "annotatr",      "dplyr",
    "ggridges", "limma",        "DMRcaller",
    "GenomicRanges",            "ComplexHeatmap"
  )
  invisible(lapply(packages_to_load, require, character.only = TRUE))
}

extractCovDf <- function(data, threshold, cov_cols) {
  stopifnot(is.data.table(data))
  stopifnot(is.numeric(threshold), length(threshold) == 1)

  if (is.numeric(cov_cols)) {
    cov_cols <- colnames(data)[cov_cols]
  }

  missing <- setdiff(cov_cols, colnames(data))
  if (length(missing) > 0) {
    stop(paste("extractCovDf: columns not found:", paste(missing, collapse = ", ")))
  }

  mask <- Reduce(`&`, lapply(cov_cols, function(col) data[[col]] >= threshold))
  return(data[mask])
}

buildMergedMatrix <- function(samplesheet,
                               sampleset,
                               datadir,
                               include_epic   = FALSE,
                               epic_path      = NULL,
                               ont_suffix     = "_modkit_pileup.bed",
                               bismark_suffix = ".bismark.cov.gz",
                               pacbio_suffix  = ".GRCh38.pbmm2.combined.bed") {

  stopifnot(is.data.table(samplesheet))
  stopifnot(sampleset %in% c("Blood", "Fibroblast", "GIAB"))
  if (include_epic && is.null(epic_path)) {
    stop("buildMergedMatrix: epic_path required when include_epic = TRUE")
  }

  samples <- samplesheet[Sampleset == sampleset, Sample]
  is_giab <- sampleset == "GIAB"
  methods <- if (is_giab) c("ONT","RRBS","WGEC","TWIST","PacBio") else
                           c("ONT","RRBS","WGEC","TWIST")

  cat(sprintf("\n[buildMergedMatrix] %s | Samples: %s\n",
    sampleset, paste(samples, collapse=", ")))

  all_sample_dts <- lapply(samples, function(smp) {

    method_dts <- lapply(methods, function(m) {

      if (m == "ONT") {
        path <- file.path(datadir, "ONT", paste0(smp, ont_suffix))
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
        dt <- readBedMethyl(path)
        dt <- dt[modbase == "m"]
        dt[, start := start + 1L]
        dt <- dt[, .(coord = paste0(chr,":",start),
                     cov   = Nvalid_cov,
                     meth  = fraction_mod)]

      } else if (m %in% c("RRBS","WGEC","TWIST")) {
        path <- file.path(datadir, m, paste0(smp, bismark_suffix))
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
        dt <- readBismarkMeth(path, correct_coords = TRUE)
        dt <- dt[, .(coord = paste0(chr,":",start),
                     cov   = coverage,
                     meth  = percentage / 100)]

      } else if (m == "PacBio") {
        path <- file.path(datadir, "PacBio", paste0(smp, pacbio_suffix))
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
        dt <- readPacBio(path)
        dt[, start := start + 1L]
        dt <- dt[, .(coord = paste0(chr,":",start),
                     cov   = coverage,
                     meth  = percentage / 100)]
      }

      dt <- dt[!duplicated(coord)]
      prefix <- ifelse(m == "WGEC", "WGBS", m)
      setnames(dt,
        c("cov", "meth"),
        c(paste0(prefix, "_cov_", smp), paste0(prefix, "_", smp))
      )
      return(dt)
    })

    method_dts <- Filter(Negate(is.null), method_dts)
    if (length(method_dts) == 0) return(NULL)
    Reduce(function(a,b) merge(a, b, by="coord", all=FALSE), method_dts)
  })

  all_sample_dts <- Filter(Negate(is.null), all_sample_dts)

  cat("  Merging across samples...\n")
  merged <- Reduce(function(a,b) merge(a, b, by="coord", all=TRUE), all_sample_dts)

  if (include_epic) {
    cat("  Joining EPIC data...\n")
    epic <- fread(epic_path, header=TRUE, sep="\t")
    if (!"coord" %in% colnames(epic)) {
      warning("EPIC matrix has no coord column – add chr:start column before joining.")
    } else {
      merged <- merge(merged, epic, by="coord", all.x=TRUE)
    }
  }

  out_dir  <- file.path(datadir, "matrices")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  suffix   <- ifelse(include_epic, "_with_EPIC.csv", "_without_EPIC.csv")
  out_file <- file.path(out_dir, paste0(sampleset, suffix))
  fwrite(merged, out_file, sep=",", quote=FALSE, na="NA")
  cat(sprintf("  Saved: %s (%d CpGs x %d columns)\n",
    out_file, nrow(merged), ncol(merged)))

  return(invisible(merged))
}

computeCorrAcrossCoverages <- function(data,
                                        samples,
                                        methods,
                                        coverages = c(0, 5, 10, 15)) {

  stopifnot(is.data.table(data))

  result <- data.frame()

  for (smp in samples) {
    for (cov in coverages) {

      if (cov == 0) {
        filtered <- data
      } else {
        cov_cols <- paste0(methods, "_cov_", smp)
        cov_cols <- cov_cols[cov_cols %in% colnames(data)]
        filtered <- extractCovDf(data, threshold = cov, cov_cols = cov_cols)
      }

      result <- computePairwiseCorr(
        data     = filtered,
        methods  = methods,
        sample   = smp,
        coverage = cov,
        existing = result
      )
    }
  }

  result$sample <- sub("^_", "", result$Sample)
  cov_labels     <- ifelse(coverages == 0, "None", paste0(coverages, "x"))
  result$Coverage2 <- ifelse(result$Coverage == 0, "None",
                              paste0(result$Coverage, "x"))
  result$Coverage2 <- factor(result$Coverage2, levels = cov_labels)

  return(result)
}
