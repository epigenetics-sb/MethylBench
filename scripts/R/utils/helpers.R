#!/usr/bin/env Rscript
# =============================================================================
# MethylBench - shared helper functions
# =============================================================================
# Readers for per-platform methylation calls, matrix construction, coverage
# filtering, correlation helpers and the shared color scheme. Sourced by all
# analysis scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- Readers ----------------------------------------------------------------

# modkit pileup (bedMethyl). Some modkit versions join the last nine fields
# with spaces instead of tabs; these are split back into 18 columns.
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

# pb-CpG-tools combined BED (header lines starting with '#' are skipped).
readPacBio <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))

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

# Bismark coverage2cytosine --merge_CpG output (*.merged_CpG_evidence.cov):
# Returns 0-based starts (position of the CpG)
readBismarkMergedCpG <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))
  
  dataset <- fread(path, header = FALSE, sep = "\t")
  
  full_colnames <- c(
    "chr", "start", "end", "percentage", "count_meth", "count_unmeth"
  )
  if (ncol(dataset) != length(full_colnames)) {
    stop(sprintf(
      "readBismarkMergedCpG: expected %d tab-separated columns (chr,start,end,percentage,count_meth,count_unmeth), got %d in file: %s",
      length(full_colnames), ncol(dataset), path
    ))
  }
  colnames(dataset) <- full_colnames
  
  width <- dataset$end - dataset$start
  if (all(width == 1L)) {
    dataset[, start := start - 1L]
  } else if (all(width == 2L)) {
    cat(sprintf(
      "  [readBismarkMergedCpG] %s: already 0-based (--zero_based), no shift applied\n",
      basename(path)
    ))
  } else {
    stop(sprintf(
      "readBismarkMergedCpG: unexpected interval widths in %s (end - start = 1: %d rows, = 2: %d rows, other: %d rows) -- is this a coverage2cytosine --merge_CpG file?",
      path, sum(width == 1L), sum(width == 2L), sum(!width %in% 1:2)
    ))
  }
  
  if (anyDuplicated(dataset, by = c("chr", "start"))) {
    stop(sprintf(
      "readBismarkMergedCpG: duplicated CpG positions in %s -- strands are not merged (run coverage2cytosine with --merge_CpG).",
      path
    ))
  }
  
  dataset[, coverage := count_meth + count_unmeth]
  return(dataset)
}

# Bismark coverage file (1-based, one row per strand).
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

# TRUE for rows where every platform has at least one non-NA sample.
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

# Space-separated RRBS table (chr, start, coverage, count_unmeth, fraction).
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

# Merges the +/- strand rows of each CpG (adjacent positions) into one row by
# summing methylated and unmethylated counts.
mergeCpGStrands <- function(dt) {
  stopifnot(is.data.table(dt))
  dt <- dt[, .(chr, start, count_meth, count_unmeth)]
  setorder(dt, chr, start)

  dt[, row_in_chr := seq_len(.N), by = chr]

  dt[, run_id := rleid(chr, start - row_in_chr)]
  dt[, row_in_chr := NULL]

  dt[, local_idx := seq_len(.N), by = run_id]
  dt[, is_partner := (local_idx %% 2L == 0L)]

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

# ---- Long-format tables -----------------------------------------------------

# Long-format methylation values per sample, method and coverage threshold. For
# thresholds > 0 a CpG is kept only if all `filter_methods` reach the threshold.
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

# Long-format per-CpG coverage (CpGs with coverage > 0).
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

# ---- Correlation ------------------------------------------------------------

# Pearson correlation of methylation values for all method pairs of one sample.
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

# ---- Colors ----------------------------------------------------------------
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
  "TWIST"  = "#CC79A7",
  "WGEC"   = "purple",
  "RRBS"   = "#0072B2"
)

get_colors <- function() {
  return(METHOD_COLORS)
}

# ---- Matrices ---------------------------------------------------------------

# Rows in which all `cov_cols` reach `threshold`.
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

# Merges the per-CpG calls of all platforms for one sample set into a matrix
# with <Method>_<Sample> (beta, 0-1) and <Method>_cov_<Sample> columns.
# Coordinates are converted to 1-based; platforms are inner-joined within a
# sample and samples are outer-joined.
buildMergedMatrix <- function(samplesheet,
                               sampleset,
                               datadir,
                               include_epic   = FALSE,
                               epic_path      = NULL,
                               ont_suffix     = "_modkit_pileup.bed",
                               bismark_suffix = ".bismark.cov.gz",
                               pacbio_suffix  = ".GRCh38.pbmm2.combined.bed",
                               pacbio_paths   = NULL) {

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
                     meth  = fraction_mod / 100)]

      } else if (m %in% c("RRBS","WGEC","TWIST")) {
        path <- file.path(datadir, m, paste0(smp, bismark_suffix))
        if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }
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

# Pairwise correlations per sample across coverage thresholds; at each threshold
# all methods of the sample must reach it.
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

# ---- UpSet helpers ----------------------------------------------------------

# Colors the platform-exclusive intersections of a ComplexUpset plot; only for
# intersections that are actually displayed (>= min_size).
singleton_queries <- function(df, sets, colors, min_size,
                              annotation = "Intersection size") {
  m <- as.matrix(df[, sets, drop = FALSE]) > 0
  shown <- vapply(sets, function(p) {
    others <- setdiff(sets, p)
    sum(m[, p] & rowSums(m[, others, drop = FALSE]) == 0) >= min_size
  }, logical(1))
  lapply(sets[shown], function(p)
    ComplexUpset::upset_query(
      intersect = p, color = colors[[p]], fill = colors[[p]],
      only_components = c("intersections_matrix", annotation)))
}
