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

  dt[, partner_meth   := data.table::shift(count_meth,   -1L), by = run_id]
  dt[, partner_unmeth := data.table::shift(count_unmeth, -1L), by = run_id]

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


# Sample sheet (comma- or tab-separated) with columns Sample and Sampleset.
readSampleSheet <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))
  ss <- fread(path, header = TRUE)
  if (!all(c("Sample", "Sampleset") %in% names(ss))) {
    stop("Sample sheet needs columns Sample and Sampleset: ", path)
  }
  as.data.frame(ss)
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
      Comparison  = pairLabel(m1_label, m2_label),
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

# Platform pairs: one canonical label per pair (fixed platform order, independent
# of the order in which methods are passed) and one fixed colour per pair, so
# that every comparison has the same label and colour in all correlation plots
# (Figure 4, Suppl. Figures 2 and 3).
PAIR_METHOD_ORDER <- c("ONT", "PacBio", "RRBS", "TWIST", "WGEC", "EPIC")

pairLabel <- function(m1, m2) {
  rk <- match(c(m1, m2), PAIR_METHOD_ORDER)
  if (anyNA(rk)) stop("pairLabel: unknown platform in '", m1, "' / '", m2, "'")
  ord <- c(m1, m2)[order(rk)]
  paste(ord[1], "vs", ord[2])
}

PAIR_COLORS <- c(
  "ONT vs PacBio"   = "dodgerblue2",
  "ONT vs RRBS"     = "#E31A1C",
  "ONT vs TWIST"    = "green4",
  "ONT vs WGEC"     = "#6A3D9A",
  "PacBio vs RRBS"  = "#FF7F00",
  "PacBio vs TWIST" = "black",
  "PacBio vs WGEC"  = "gold1",
  "RRBS vs TWIST"   = "skyblue2",
  "RRBS vs WGEC"    = "#FB9A99",
  "TWIST vs WGEC"   = "palegreen3"
)

# Colour scale for platform pairs; only the pairs present in the data are
# shown in the legend, always in the same order.
scale_color_pairs <- function(...) {
  scale_color_manual(values = PAIR_COLORS,
                     limits = function(x) intersect(names(PAIR_COLORS), x), ...)
}

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

# ---- Per-platform calls and merged matrices ---------------------------------

# Bismark coverage file of any of the two layouts: per-strand .cov
# (end == start, 1-based; strands are merged with mergeCpGStrands()) or
# coverage2cytosine --merge_CpG output (end - start == 1: 1-based;
# end - start == 2: 0-based, --zero_based). Returns chr, start (0-based
# position of the C on the + strand), count_meth, count_unmeth.
readBismarkAuto <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))
  d <- fread(path, header = FALSE, sep = "\t")
  if (ncol(d) != 6L) stop("readBismarkAuto: expected 6 columns, got ", ncol(d), " in ", path)
  setnames(d, c("chr", "start", "end", "percentage", "count_meth", "count_unmeth"))
  w <- unique(d$end - d$start)
  if (all(w == 0L)) {
    d <- mergeCpGStrands(d[, .(chr, start, count_meth, count_unmeth)])
    d[, start := as.integer(start) - 1L]
  } else if (all(w == 1L)) {
    d[, start := as.integer(start) - 1L]
  } else if (all(w == 2L)) {
    d[, start := as.integer(start)]
  } else {
    stop("readBismarkAuto: unexpected interval widths in ", path,
         " -- expected a per-strand .cov or a coverage2cytosine --merge_CpG file")
  }
  d[, .(chr, start, count_meth, count_unmeth)]
}

# Default file suffixes of the per-platform calls (<calls_dir>/<METHOD>/<Sample><suffix>).
CALL_SUFFIXES <- list(ONT = "_modkit_pileup.bed", BISMARK = ".bismark.cov.gz",
                      RRBS = ".methylation.tab", PacBio = ".combined.bed")

# Per-platform methylation calls of one sample, harmonized to one row per CpG:
# chr (UCSC style), start (0-based position of the C on the + strand, BED
# convention, as in the merged matrices; both strands combined), cov (reads)
# and beta (0-1). CpGs with zero coverage are dropped. Returns NULL
# if no file is found; attr(, "src") gives the file read.
#   ONT     modkit pileup bedMethyl (--cpg --combine-strands), 0-based, 5mC (code m) only
#   WGEC, TWIST
#           Bismark coverage, per-strand or coverage2cytosine --merge_CpG (readBismarkAuto)
#   RRBS    lab-specific table (chr pos coverage n_unmeth fraction; 1-based, per
#           strand; strands merged) if present, otherwise the Bismark file
#   PacBio  pb-CpG-tools combined BED, 0-based (strands combined by the tool);
#           rows of type "Total"; beta = mod_score / 100 (pacbio_value =
#           "mod_score", as for the matrices of the manuscript) or
#           discretized_mod_score / 100 ("discretized")
# With `file`, that file is read instead of <calls_dir>/<platform>/<sample><suffix>;
# an RRBS file ending in suffixes$RRBS is read as the lab-specific table.
readPlatformCalls <- function(calls_dir, platform, sample, suffixes = CALL_SUFFIXES,
                              file = NULL, pacbio_value = c("mod_score", "discretized")) {
  pacbio_value <- match.arg(pacbio_value)
  call_path <- function(sfx) {
    if (!is.null(file)) return(file)
    file.path(calls_dir, platform, paste0(sample, sfx))
  }
  is_rrbs_table <- function(p) endsWith(p, suffixes$RRBS)
  d <- NULL; src <- NA_character_
  counts_to_beta <- function(x) x[, .(chr, start, cov = count_meth + count_unmeth,
                                     beta = fifelse(count_meth + count_unmeth > 0,
                                                    count_meth / (count_meth + count_unmeth), NA_real_))]
  if (platform == "ONT") {
    p <- call_path(suffixes$ONT)
    if (file.exists(p)) {
      d <- readBedMethyl(p)[modbase == "m"]
      d <- d[, .(chr, start = as.integer(start), cov = Nvalid_cov, beta = fraction_mod / 100)]
      src <- p
    }
  } else if (platform == "RRBS" && file.exists(call_path(suffixes$RRBS)) &&
             is_rrbs_table(call_path(suffixes$RRBS))) {
    p <- call_path(suffixes$RRBS)
    d <- mergeCpGStrands(readRRBSTab(p))
    d[, start := as.integer(start) - 1L]   # 1-based table -> 0-based
    d <- counts_to_beta(d); src <- p
  } else if (platform %in% c("RRBS", "WGEC", "TWIST")) {
    p <- call_path(suffixes$BISMARK)
    if (file.exists(p)) { d <- counts_to_beta(readBismarkAuto(p)); src <- p }
  } else if (platform == "PacBio") {
    p <- call_path(suffixes$PacBio)
    if (file.exists(p)) {
      d <- readPacBio(p)
      if (any(d$haplotype == "Total")) d <- d[haplotype == "Total"]
      d <- d[, .(chr, start = as.integer(start), cov = coverage,
                 beta = if (pacbio_value == "mod_score") score / 100 else percentage / 100)]
      src <- p
    }
  } else {
    stop("readPlatformCalls: unknown platform ", platform)
  }
  if (is.null(d)) return(NULL)
  d[, chr := as.character(chr)]
  d[!startsWith(chr, "chr"), chr := paste0("chr", sub("^MT$", "M", chr))]
  d <- d[cov > 0]
  d <- d[!duplicated(d, by = c("chr", "start"))]
  setattr(d, "src", src)
  d
}

# EPIC beta matrix (tab- or comma-separated) with a 'coord' column (chr:pos) or
# 'chr' and 'start' columns and one column per sample named EPIC_<Sample>.
# Positions are used as given. Returns chr, start and the EPIC_<Sample> columns;
# probes at the same position are averaged.
readEPICMatrix <- function(path) {
  stopifnot(is.character(path), length(path) == 1)
  if (!file.exists(path)) stop(paste("File not found:", path))
  e <- fread(path)
  if ("coord" %in% names(e)) {
    e[, c("chr", "start") := tstrsplit(coord, ":", fixed = TRUE)]
    e[, start := as.integer(start)]
  } else if (!all(c("chr", "start") %in% names(e))) {
    stop("readEPICMatrix: needs a 'coord' column or 'chr' and 'start' columns: ", path)
  }
  cols <- grep("^EPIC_", names(e), value = TRUE)
  if (length(cols) == 0) stop("readEPICMatrix: no EPIC_<Sample> columns in ", path)
  e <- e[, c("chr", "start", cols), with = FALSE]
  e[, chr := as.character(chr)]
  e[!startsWith(chr, "chr"), chr := paste0("chr", chr)]
  e[, start := as.integer(start)]
  # Probes mapping to the same position (e.g. EPIC v2 replicates) are averaged.
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  e[, lapply(.SD, mean_or_na), by = .(chr, start), .SDcols = cols]
}

# Table of input files (comma- or tab-separated; columns Sample, Method, File),
# an alternative to the <calls_dir>/<METHOD>/<Sample><suffix> layout.
readCallsSheet <- function(path) {
  cs <- fread(path)
  if (!all(c("Sample", "Method", "File") %in% names(cs))) stop("Calls sheet needs columns Sample, Method, File: ", path)
  if (anyDuplicated(cs, by = c("Sample", "Method"))) stop("Calls sheet: duplicated Sample/Method in ", path)
  miss <- cs[!file.exists(File), File]
  if (length(miss) > 0) stop("Files not found:\n", paste(miss, collapse = "\n"))
  cs
}
# File for one sample and platform: NULL without a calls sheet (directory
# layout is used), "" if the sheet has no entry (treated as missing).
callsSheetFile <- function(cs, method, sample) {
  if (is.null(cs)) return(NULL)
  f <- cs[Sample == sample & Method == method, File]
  if (length(f) == 0) "" else f
}

# Joins a list of data.tables on chr/start; how = "inner" or "outer".
joinOnCpG <- function(dts, how = c("inner", "outer")) {
  how <- match.arg(how)
  dts <- Filter(Negate(is.null), dts)
  if (length(dts) == 0) return(NULL)
  Reduce(function(a, b) merge(a, b, by = c("chr", "start"), all = (how == "outer")), dts)
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
