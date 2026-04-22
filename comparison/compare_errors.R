#!/usr/bin/env Rscript
# compare_errors.R
#
# Run R's learnErrors() on a filtered FASTQ (the same file produced by the
# Rust filter-and-trim step) and compare the resulting error matrices against
# the Rust errors-from-sample output.
#
# Usage:
#   Rscript compare_errors.R <filtered.fastq.gz> <errors_rust.json> [out_prefix]
#
# Arguments:
#   filtered.fastq.gz — FASTQ produced by `dada2-rs filter-and-trim`
#                       (used as-is by R's derepFastq, guaranteeing identical
#                        input sequences to both pipelines)
#   errors_rust.json  — Rust `errors-from-sample` output
#   out_prefix        — optional output prefix for CSVs and R error JSON
#                       (default: same dir/stem as errors_rust.json)
#
# Filter parameters applied by run_rust_errors.sh (MiSeq SOP, forward reads):
#   truncLen=240, maxN=0, maxEE=2, truncQ=2, rm.phix=FALSE
#   (rm.phix omitted here so both pipelines see the same reads)

suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compare_errors.R <filtered.fastq.gz> <errors_rust.json> [out_prefix]\n")
  quit(status = 1)
}

fastq_path <- args[1]
rust_path  <- args[2]
out_prefix <- if (length(args) >= 3) args[3] else sub("\\.json$", "", rust_path)

# ---------------------------------------------------------------------------
# Step 1: Dereplicate the filtered FASTQ with R (informational only)
# ---------------------------------------------------------------------------
cat("==> R derepFastq:", fastq_path, "\n")
derep_r <- derepFastq(fastq_path)
cat(sprintf(
  "    %d unique sequences, %d total reads\n",
  length(derep_r$uniques),
  sum(derep_r$uniques)
))

# ---------------------------------------------------------------------------
# Step 2: Run R learnErrors() with parameters matching Rust / R defaults
#
# learnErrors() is the intended high-level API; it calls dada(selfConsist=TRUE)
# internally and returns getErrors(dds, detailed=TRUE) in the correct matrix
# format.  Calling dada() + getErrors() directly can return unexpected shapes
# for $trans and $err_in depending on dada2 version.
#
# Rust cli.rs (run_rust_errors.sh):    R dada_opts defaults:
#   omega_a   = 1e-40                    OMEGA_A    = 1e-40  (match)
#   omega_c   = 1e-40                    OMEGA_C    = 1e-40  (match)
#   omega_p   = 1e-4                     OMEGA_P    = 1e-4   (match)
#   min_fold  = 1.0                      MIN_FOLD   = 1      (match)
#   min_hamming = 1                      MIN_HAMMING = 1     (match)
#   min_abund = 1                        MIN_ABUNDANCE = 1   (match)
#   max_consist = 10                     MAX_CONSIST = 10    (match)
#   errfun = loess                       loessErrfun         (match)
#
# nbases: the filtered file is ~1.7 M bases, well below the 1e8 default,
# so all reads are used without subsampling.
# ---------------------------------------------------------------------------
cat("==> Running R learnErrors()...\n")
err_r <- learnErrors(
  fastq_path,
  errorEstimationFunction = loessErrfun,
  multithread             = FALSE,
  OMEGA_A                 = 1e-40,
  OMEGA_C                 = 1e-40,
  OMEGA_P                 = 1e-4,
  MIN_FOLD                = 1,
  MIN_HAMMING             = 1L,
  MIN_ABUNDANCE           = 1L,
  USE_QUALS               = TRUE,
  GREEDY                  = TRUE,
  MAX_CONSIST             = 10
)

# Inspect what learnErrors actually returned so we can coerce correctly.
cat("    learnErrors() returned fields:", paste(names(err_r), collapse=", "), "\n")
for (nm in names(err_r)) {
  x <- err_r[[nm]]
  cat(sprintf("      $%-10s  class=%-12s  dim=%s  length=%d\n",
              nm, class(x)[1],
              if (is.null(dim(x))) "NULL" else paste(dim(x), collapse="x"),
              length(x)))
}

# Coerce each field to a 16 × nq matrix regardless of how it was stored.
# $err_in comes back as a list of matrices (one per selfConsist round);
# take the last element, which is the error model fed into the final DADA run.
to_16_matrix <- function(x, nq) {
  if (is.null(x))      return(matrix(NA_real_, nrow = 16, ncol = nq))
  if (is.list(x))      x <- x[[length(x)]]   # last iteration for err_in
  if (is.null(dim(x))) return(matrix(x, nrow = 16))
  x
}

nq_r      <- ncol(err_r$err_out)
err_out_r <- to_16_matrix(err_r$err_out, nq_r)
err_in_r  <- to_16_matrix(err_r$err_in,  nq_r)
trans_r   <- to_16_matrix(err_r$trans,   nq_r)

cat(sprintf("    R: nq=%d\n", nq_r))

# ---------------------------------------------------------------------------
# Step 3: Load the Rust error model JSON
# ---------------------------------------------------------------------------
cat("==> Loading Rust error model JSON:", rust_path, "\n")
rj    <- fromJSON(rust_path)
nq_rs <- rj$nq
cat(sprintf("    Rust: nq=%d, converged=%s, iterations=%d\n",
            nq_rs, rj$converged, rj$iterations))

# Rust JSON: each matrix is a list of 16 row-vectors of length nq (row-major).
to_matrix <- function(lst, nq) {
  matrix(unlist(lst), nrow = 16, ncol = nq, byrow = TRUE)
}
trans_rs   <- to_matrix(rj$trans,   nq_rs)
err_in_rs  <- to_matrix(rj$err_in,  nq_rs)
err_out_rs <- to_matrix(rj$err_out, nq_rs)

# ---------------------------------------------------------------------------
# Align nq (take the overlap if they differ)
# ---------------------------------------------------------------------------
nq <- min(nq_r, nq_rs)
if (nq_r != nq_rs) {
  cat(sprintf(
    "  WARNING: nq mismatch (R=%d, Rust=%d); comparing first %d columns.\n",
    nq_r, nq_rs, nq
  ))
}
trans_r    <- trans_r[,    seq_len(nq), drop = FALSE]
err_in_r   <- err_in_r[,  seq_len(nq), drop = FALSE]
err_out_r  <- err_out_r[, seq_len(nq), drop = FALSE]
trans_rs   <- trans_rs[,   seq_len(nq), drop = FALSE]
err_in_rs  <- err_in_rs[,  seq_len(nq), drop = FALSE]
err_out_rs <- err_out_rs[, seq_len(nq), drop = FALSE]

# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------
nts        <- c("A", "C", "G", "T")
row_labels <- paste0(rep(nts, each = 4), "->", rep(nts, times = 4))

compare_mat <- function(A, B, label) {
  # Focus on cells where at least one value is non-trivial
  mask <- (A > 1e-10 | B > 1e-10)
  d    <- abs(A - B)[mask]
  cat(sprintf(
    "  %-10s  max|diff|=%.3e  mean|diff|=%.3e  median|diff|=%.3e\n",
    label,
    max(d,    na.rm = TRUE),
    mean(d,   na.rm = TRUE),
    median(d, na.rm = TRUE)
  ))
  invisible(abs(A - B))
}

compare_trans <- function(A, B) {
  total <- A + B
  mask  <- total > 0
  rel   <- abs(A - B)[mask] / (total[mask] / 2)
  cat(sprintf(
    "  %-10s  max_rel=%.3e  mean_rel=%.3e  max_abs=%.0f  total_R=%d  total_Rust=%d\n",
    "trans",
    max(rel,  na.rm = TRUE),
    mean(rel, na.rm = TRUE),
    max(abs(A - B)),
    sum(A),
    sum(B)
  ))
}

# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
cat("\n=== Matrix comparison (R vs Rust) ===\n")
cat(sprintf("  nq=%d, 16 transition rows\n\n", nq))

compare_trans(trans_r, trans_rs)
compare_mat(err_in_r,  err_in_rs,  "err_in")
compare_mat(err_out_r, err_out_rs, "err_out")

cat("\n--- err_out max|diff| per transition type ---\n")
for (r in seq_len(16)) {
  d <- max(abs(err_out_r[r, ] - err_out_rs[r, ]), na.rm = TRUE)
  cat(sprintf("  %s  %.3e\n", row_labels[r], d))
}

# ---------------------------------------------------------------------------
# Export comparison CSVs
# ---------------------------------------------------------------------------
export_comparison <- function(mat_r, mat_rs, name) {
  quality <- seq(0, nq - 1)
  rows <- lapply(seq_len(16), function(r) {
    a <- mat_r[r, ]
    b <- mat_rs[r, ]
    data.frame(
      transition = row_labels[r],
      quality    = quality,
      R          = a,
      Rust       = b,
      abs_diff   = abs(a - b),
      rel_diff   = ifelse((a + b) > 0, abs(a - b) / ((a + b) / 2), NA_real_)
    )
  })
  df   <- do.call(rbind, rows)
  path <- paste0(out_prefix, "_compare_", name, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat("  Written:", path, "\n")
  invisible(df)
}

cat("\n==> Exporting CSVs...\n")
export_comparison(trans_r,   trans_rs,   "trans")
export_comparison(err_in_r,  err_in_rs,  "err_in")
export_comparison(err_out_r, err_out_rs, "err_out")

# Export R's error model in Rust JSON format so plot_errors.R can compare them
r_out <- list(
  nq         = nq_r,
  converged  = TRUE,
  iterations = nrow(err_r$err_in),   # number of selfConsist rounds
  trans      = lapply(seq_len(16), function(r) as.integer(trans_r[r, ])),
  err_in     = lapply(seq_len(16), function(r) err_in_r[r, ]),
  err_out    = lapply(seq_len(16), function(r) err_out_r[r, ])
)
r_json_path <- paste0(out_prefix, "_errors_r.json")
write(toJSON(r_out, auto_unbox = TRUE, digits = 10), r_json_path)
cat("  Written:", r_json_path, "\n")
cat("  (Run plot_errors.R on the R and Rust JSON files to compare visually.)\n")

cat("\nDone.\n")
