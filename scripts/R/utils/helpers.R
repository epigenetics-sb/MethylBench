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

  full_colnames <- c(
    "chr", "start", "end", "modbase", "score", "strand",
    "start1", "end1", "color", "Nvalid_cov", "fraction_mod",
    "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
    "Nfail", "Ndiff", "Nnocall"
  )

  if (ncol(dataset) == 10L) {
    # Malformed export seen in some GIAB modkit files: the first 9 fields
    # are properly tab-separated, but the remaining 9 numeric fields got
    # joined by spaces into a single tab-field instead of also being
    # tab-separated. Split that last column further on whitespace.
    split_vals <- tstrsplit(dataset[[10]], "\\s+", perl = TRUE)
    if (length(split_vals) != 9L) {
      stop(sprintf(
        "readBedMethyl: column 10 in %s did not split into the expected 9 whitespace-separated values (got %d) -- inspect the file's exact format before proceeding.",
        path, length(split_vals)
      ))
    }
    extra_dt <- as.data.table(lapply(split_vals, as.numeric))
    dataset  <- cbind(dataset[, 1:9], extra_dt)
    cat(sprintf(
      "  [readBedMethyl] %s: repaired mixed tab/space-delimited format (9 tab fields + 1 space-joined blob -> 18 columns)\n",
      basename(path)
    ))
  }

  if (ncol(dataset) != 18L) {
    stop(sprintf(
      "readBedMethyl: expected 18 columns (or 9 tab fields + 1 space-joined blob of 9 values), got %d in file: %s -- refusing to guess the column layout.",
      ncol(dataset), path
    ))
  }

  colnames(dataset) <- full_colnames
  return(dataset)
}

readPacBio <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  # pb-cpg-tools output starts with several "##key=value" metadata lines
  # (VCF-style), e.g. "##pb-cpg-tools-version=3.0.0". Strip any leading
  # comment lines before parsing, regardless of how many there are, and
  # regardless of .bed vs .bed.gz.
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  cmd <- if (is_gz) {
    sprintf("zcat %s | grep -v '^#'", shQuote(path))
  } else {
    sprintf("grep -v '^#' %s", shQuote(path))
  }
  dataset <- fread(cmd = cmd, header = FALSE, sep = "\t")

  full_colnames <- c(
    "chr", "start", "end", "score", "haplotype",
    "coverage", "N_modified", "N_unmodified", "percentage"
  )
  if (ncol(dataset) != length(full_colnames)) {
    stop(sprintf(
      "readPacBio: expected %d columns after stripping '#' header lines, got %d in file: %s -- refusing to guess the column layout.",
      length(full_colnames), ncol(dataset), path
    ))
  }
  colnames(dataset) <- full_colnames

  if (!is.numeric(dataset$start)) {
    stop(sprintf(
      "readPacBio: 'start' column is non-numeric in %s even after stripping '#' lines -- inspect the file's exact format.",
      path
    ))
  }

  return(dataset)
}

#' Read a raw Bismark coverage file (.bismark.cov / .bismark.cov.gz), as
#' produced by coverage2cytosine WITHOUT --merge_CpG. Format is tab-separated,
#' 6 columns, 1-based single-base coordinates:
#'   chr  start  end  percentage(0-100)  count_methylated  count_unmethylated
#' NOTE: start/end are already 1-based -- do NOT apply any +1/-1 shift here.
#' Because --merge_CpG was not used, the +strand C of a CpG (position N) and
#' the -strand C (position N+1) appear as two SEPARATE rows; use
#' mergeCpGStrands() afterwards to combine them into one row per CpG.
#' @return data.table with columns chr, start, count_meth, count_unmeth.
readBismarkCovRaw <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = FALSE, sep = "\t")

  expected_cols <- 6L
  if (ncol(dataset) != expected_cols) {
    stop(sprintf(
      "readBismarkCovRaw: expected %d tab-separated columns (chr,start,end,percentage,count_meth,count_unmeth), got %d in file: %s",
      expected_cols, ncol(dataset), path
    ))
  }

  colnames(dataset) <- c(
    "chr", "start", "end", "percentage", "count_meth", "count_unmeth"
  )

  dataset[, .(chr, start, count_meth, count_unmeth)]
}

#' Read the lab's custom RRBS .tab format, produced from raw Bismark
#' .cov.gz via: awk '{print $1,$2,$6+$5,$6,$4/100}'
#' i.e. columns: chr, start, coverage, count_unmeth, fraction(0-1).
#' start is untouched from the raw .cov.gz, which is already 1-based -- do
#' NOT apply any coordinate shift. Like WGEC/TWIST, this is NOT strand-
#' merged (derived from the same un-merged raw Bismark output) -- follow
#' with mergeCpGStrands() after computing count_meth.
#' @return data.table with columns chr, start, count_meth, count_unmeth.
#' Restrict a wide methylation matrix to CpGs "covered by every platform",
#' matching the paper's stated exploratory-analysis restriction ("all
#' analyses were restricted to the common set of CpG sites covered across
#' platforms"). A CpG counts as covered by a platform if at least one
#' sample column for that platform is non-NA at that row -- this is
#' deliberately more lenient than the Tier1/Tier2 10x-coverage consensus
#' sets used elsewhere in the pipeline, which is appropriate for this
#' specific exploratory step (see Methods: limma/Wilcoxon on n=146,704
#' common CpGs, distinct from the later coverage-matched consensus sets).
#'
#' @param dt data.table containing the method's sample columns.
#' @param method_prefixes character vector of method name prefixes (e.g.
#'   c("EPIC","ONT","WGEC","TWIST","RRBS")); columns are matched by
#'   `^<prefix>_` (case-insensitive), excluding any `_cov_` columns.
#' @return logical vector, one per row of dt, TRUE where every prefix has
#'   at least one non-NA sample column.
platformsCoveredMask <- function(dt, method_prefixes) {
  covered <- matrix(TRUE, nrow = nrow(dt), ncol = length(method_prefixes))
  for (i in seq_along(method_prefixes)) {
    p <- method_prefixes[i]
    cols <- grep(paste0("(?i)^", p, "_(?!cov_)"), colnames(dt), value = TRUE, perl = TRUE)
    if (length(cols) == 0) {
      covered[, i] <- FALSE
      next
    }
    sub <- dt[, ..cols]
    covered[, i] <- rowSums(!is.na(sub)) > 0
  }
  apply(covered, 1, all)
}

readRRBSTab <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

  dataset <- fread(path, header = FALSE, sep = " ")

  expected_cols <- 5L
  if (ncol(dataset) != expected_cols) {
    stop(sprintf(
      "readRRBSTab: expected %d space-separated columns (chr,start,coverage,count_unmeth,fraction), got %d in file: %s",
      expected_cols, ncol(dataset), path
    ))
  }

  colnames(dataset) <- c("chr", "start", "coverage", "count_unmeth", "fraction")
  dataset[, count_meth := coverage - count_unmeth]

  dataset[, .(chr, start, count_meth, count_unmeth)]
}


#'
#' Bismark (without --merge_CpG) reports the +strand C of a CpG dinucleotide
#' at position N and the -strand C at position N+1 as two independent rows.
#' A naive "does start-1 exist" check is NOT sufficient to pair them
#' correctly when CpGs are directly adjacent (e.g. positions N, N+1, N+2 all
#' present, where N+2 is actually the NEXT CpG's own +strand, not a partner
#' of N+1) -- that would silently mismatch unrelated CpGs.
#'
#' This instead does a greedy left-to-right pairing within each maximal run
#' of consecutive integer positions (per chromosome): (1st,2nd), (3rd,4th),
#' etc. If a run has odd length, the trailing element is kept as an
#' unpaired singleton (rather than merged with an unrelated neighbor) --
#' this happens when one strand of a CpG had zero coverage and Bismark
#' therefore never emitted a row for it.
#'
#' Validated against a reference sequential implementation across pair/
#' singleton/run-of-3/run-of-4/run-of-5/multi-chromosome-boundary cases.
#'
#' @param dt data.table with columns chr, start, count_meth, count_unmeth
#'   (one row per Bismark-reported strand-specific cytosine call).
#' @return data.table with columns chr, start, count_meth, count_unmeth,
#'   one row per CpG, anchored at the lower (+strand) position.
mergeCpGStrands <- function(dt) {
  stopifnot(is.data.table(dt))
  dt <- dt[, .(chr, start, count_meth, count_unmeth)]
  setorder(dt, chr, start)

  # Row index within each chromosome (resets at every chr boundary).
  dt[, row_in_chr := seq_len(.N), by = chr]

  # While start increases by exactly 1 per row, (start - row_in_chr) stays
  # constant -- so this value changes exactly at the boundaries of maximal
  # runs of consecutive positions. rleid() also always bumps on chr change.
  dt[, run_id := rleid(chr, start - row_in_chr)]
  dt[, row_in_chr := NULL]

  # Within each run, alternate anchor (odd local_idx) / partner (even
  # local_idx): (1st,2nd) is a pair, (3rd,4th) is the next pair, etc.
  dt[, local_idx := seq_len(.N), by = run_id]
  dt[, is_partner := (local_idx %% 2L == 0L)]

  # For each anchor, look up the immediately following row's counts WITHIN
  # THE SAME RUN (shift() with by= returns NA at run boundaries, which
  # correctly identifies an anchor with no partner).
  dt[, partner_meth   := shift(count_meth,   -1L), by = run_id]
  dt[, partner_unmeth := shift(count_unmeth, -1L), by = run_id]

  anchors <- dt[is_partner == FALSE]
  anchors[, count_meth   := count_meth   + ifelse(is.na(partner_meth),   0, partner_meth)]
  anchors[, count_unmeth := count_unmeth + ifelse(is.na(partner_unmeth), 0, partner_unmeth)]

  n_pairs     <- sum(!is.na(anchors$partner_meth))
  n_singleton <- nrow(anchors) - n_pairs
  cat(sprintf(
    "    mergeCpGStrands: %d input rows -> %d CpG positions (%d strand-pairs merged, %d singleton/unpaired)\n",
    nrow(dt), nrow(anchors), n_pairs, n_singleton
  ))

  anchors[, .(chr, start, count_meth, count_unmeth)]
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

#' Build a long-format per-CpG coverage table across samples and methods,
#' analogous to buildMethLong() but reading raw <Method>_cov_<Sample> columns
#' instead of <Method>_<Sample> methylation fractions.
#'
#' @param data data.table containing <Method>_cov_<Sample> columns (e.g. the
#'   merged CpG-level matrix produced by buildMergedMatrix()).
#' @param samples character vector of sample labels (e.g. samplesheet$Sample).
#' @param methods named character vector, names = display labels, values =
#'   column-name prefixes, e.g. c(ONT="ONT", RRBS="RRBS", WGEC="WGEC",
#'   TWIST="TWIST", PacBio="PacBio"). A method/sample combination whose
#'   column doesn't exist (e.g. PacBio for non-GIAB samples) is silently
#'   skipped, since PacBio is only present for GIAB1/GIAB2.
#' @return data.table with columns Coverage (numeric, raw per-CpG coverage,
#'   zero/NA entries dropped), Method (factor, levels = names(methods)),
#'   Sample (character).
buildCovLong <- function(data, samples, methods) {

  stopifnot(is.data.table(data))
  stopifnot(is.character(samples), length(samples) >= 1)
  stopifnot(is.character(methods), !is.null(names(methods)))

  rows <- vector("list", length(samples) * length(methods))
  idx  <- 1L

  for (smp in samples) {
    for (method_label in names(methods)) {
      cov_col <- paste0(methods[method_label], "_cov_", smp)

      if (!cov_col %in% colnames(data)) {
        # e.g. PacBio_cov_<sample> only exists for GIAB1/GIAB2 -- skip
        # quietly rather than erroring, mirroring the sparse method
        # coverage across sample sets (see buildMergedMatrix()).
        next
      }

      cov_values <- data[[cov_col]]
      cov_values <- cov_values[!is.na(cov_values) & cov_values > 0]
      if (length(cov_values) == 0) next

      rows[[idx]] <- data.table(
        Coverage = cov_values,
        Method   = method_label,
        Sample   = smp
      )
      idx <- idx + 1L
    }
  }

  result <- rbindlist(rows[seq_len(idx - 1L)])
  result[, Method := factor(Method, levels = names(methods))]

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
  "black",   "gold1",   "skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6",     "#FDBF6F", "gray70", "khaki2", "maroon",
  "orchid1",     "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1",      "yellow4", "yellow3", "darkorange4", "brown"
)

C12 <- c(
  "dodgerblue2", "#E31A1C", "green4",      "#6A3D9A",
  "#FF7F00",     "black",    "gold1",        "khaki2",
  "brown",       "darkturquoise", "palegreen2", "orchid1"
)

METHOD_COLORS <- c(
  "ONT"    = "#D69F00",
  "PacBio" = "#E55E00",
  "EPIC"   = "#009E73",
  "TWIST"  = "#DC79A7",
  "WGEC"   = "purple",
  "RRBS"   = "#0072B2"
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
                               pacbio_suffix  = ".GRCh38.pbmm2.combined.bed",
                               pacbio_paths   = NULL) {
  # pacbio_paths: optional named character vector mapping this pipeline's
  # sample label (e.g. "GIAB1") to the EXACT raw PacBio file path (e.g.
  # ".../HG001.GRCh38.cpg_pileup.combined.bed"). Use this when the raw
  # PacBio files use different sample names than the rest of the pipeline
  # and/or mix compressed (.bed.gz) and uncompressed (.bed) files -- fread()
  # transparently handles either. If a sample isn't in pacbio_paths, falls
  # back to the datadir/PacBio/<sample><pacbio_suffix> convention.

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
        dt[, start := start + 1L]  # 0-based BED -> 1-based
        dt <- dt[, .(coord = paste0(chr,":",start),
                     cov   = Nvalid_cov,
                     meth  = fraction_mod / 100)]  # fraction_mod is 0-100

      } else if (m %in% c("RRBS","WGEC","TWIST")) {
        path <- file.path(datadir, m, paste0(smp, bismark_suffix))
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
        # Raw Bismark .cov.gz is already 1-based -- NO coordinate shift.
        # It also reports +/- strand of each CpG as separate rows, so
        # these must be merged before computing per-CpG methylation.
        dt_raw <- readBismarkCovRaw(path)
        dt_cpg <- mergeCpGStrands(dt_raw)
        cov_total <- dt_cpg$count_meth + dt_cpg$count_unmeth
        dt <- dt_cpg[, .(
          coord = paste0(chr, ":", start),
          cov   = cov_total,
          meth  = fifelse(cov_total > 0, count_meth / cov_total, NA_real_)
        )]

      } else if (m == "PacBio") {
        path <- if (!is.null(pacbio_paths) && smp %in% names(pacbio_paths)) {
          pacbio_paths[[smp]]
        } else {
          file.path(datadir, "PacBio", paste0(smp, pacbio_suffix))
        }
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
        dt <- readPacBio(path)
        dt[, start := start + 1L]  # 0-based BED -> 1-based
        dt <- dt[, .(coord = paste0(chr,":",start),
                     cov   = coverage,
                     meth  = percentage / 100)]
      }

      dt <- dt[!duplicated(coord)]
      prefix <- m
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
