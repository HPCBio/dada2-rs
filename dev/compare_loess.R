#!/usr/bin/env Rscript
# compare_loess.R
#
# Isolate the errfun fit difference between dada2-rs's native Rust
# implementation and R's stock stats::loess on a single, identical trans
# matrix.
#
# Reads one or more `learn-errors` JSON files. From each file pulls out the
# accumulated `trans` matrix and the embedded `err_out` (= whatever errfun
# produced this JSON), then re-fits the matching R reference (loessErrfun for
# `--errfun loess`, PacBioErrfun for `--errfun pacbio`) on the same trans and
# reports per-cell diffs.
#
# Surface (direct vs interpolate) is read from the JSON's
# `params.loess.surface` field (added in the loess-provenance commit).
# Older JSONs without that field fall back to "direct".
#
# Usage:
#   Rscript dev/compare_loess.R <learn_errors.json> [more.json ...]
#
# Interpretation:
#   - If the JSON was produced with native Rust loess (default learn-errors),
#     the diff reveals the per-cell discrepancy between Rust and R loess.
#   - For pacbio JSONs, q=93 is fit by a closed-form Laplace MLE in both
#     implementations and should match to ~machine epsilon; the loess columns
#     (q<93) are where any divergence lives.  The summary splits them so the
#     two estimators are not lumped together.
#   - If the JSON was produced with --errfun-cmd "Rscript loess_reference.R",
#     the diff should be ~ machine epsilon (sanity check on this script).
#   - If the JSON has no `trans` block (e.g. produced by
#     learnerrors_to_dada2rs.R), the file is skipped with a note.

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript compare_loess.R <learn_errors.json> [more.json ...]")
}

ROW_NAMES <- paste0(rep(c("A", "C", "G", "T"), each = 4L), "2",
                    rep(c("A", "C", "G", "T"), times = 4L))

# Verbatim port of dada2:::loessErrfun (same code as
# examples/external_errfun/loess_reference.R).
loessErrfun <- function(trans, surface = "interpolate") {
  qq  <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs  <- trans[paste0(nti, "2", ntj), ]
        tot   <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df    <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = tot, surface = surface)
        pred   <- predict(mod.lo, qq)
        # surface="direct" extrapolates the polynomial outside the data
        # range; mirror dada2-rs's None-then-flat-fill by NA-ing those.
        if (surface == "direct") {
          valid_q <- df$q[!is.na(df$rlogp) & df$tot > 0]
          pred[qq < min(valid_q) | qq > max(valid_q)] <- NA
        }

        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  # Post-fit clamp from dada2 errorModels.R:53-56 (R DADA2's `# HACKY` step).
  # Off-diagonals pinned to [1e-7, 0.25] before the diagonal is computed.
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  err <- rbind(1 - colSums(est[1:3, ]),  est[1:3, ],
               est[4, ],     1 - colSums(est[4:6, ]),  est[5:6, ],
               est[7:8, ],   1 - colSums(est[7:9, ]),  est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- ROW_NAMES
  colnames(err) <- colnames(trans)
  err
}

# Verbatim port of dada2::PacBioErrfun (R/errorModels.R:183-196), parameterized
# by `surface` so it can mirror either dada2-rs preset.  The q<93 columns are
# fit with loessErrfun and inherit its [1e-7, 0.25] clamp; the q=93 column is
# a Laplace MLE (count+1)/(total+4) and is deliberately *not* clamped — this
# matches both the R reference and pacbio_errfun in error_models.rs.
pacbioErrfun <- function(trans, surface = "interpolate") {
  if ("93" %in% colnames(trans)) {
    i.93 <- which(colnames(trans) %in% "93")
    if (i.93 != ncol(trans)) stop("Q93 must be the last column")
    err <- loessErrfun(trans[, 1:(i.93 - 1), drop = FALSE], surface)
    tot93 <- rep(c(sum(trans[1:4,  "93"]), sum(trans[5:8,  "93"]),
                   sum(trans[9:12, "93"]), sum(trans[13:16,"93"])), each = 4)
    err93 <- (trans[, "93"] + 1) / (tot93 + 4)
    err <- cbind(err, "93" = err93)
    rownames(err) <- ROW_NAMES
  } else {
    message("  PacBio JSON without q=93 column — falling back to plain loess")
    err <- loessErrfun(trans, surface)
  }
  err
}

list_to_mat <- function(x) {
  m <- do.call(rbind, lapply(x, as.numeric))
  rownames(m) <- ROW_NAMES
  colnames(m) <- as.character(seq_len(ncol(m)) - 1L)
  m
}

summarize_diff <- function(d, err_orig, label) {
  cat(sprintf("\n  per-cell |diff| (%s):\n", label))
  cat(sprintf("    max     = %.3e\n", max(d)))
  cat(sprintf("    mean    = %.3e\n", mean(d)))
  cat(sprintf("    median  = %.3e\n", median(d)))
  cat(sprintf("    p99     = %.3e\n", quantile(d, 0.99, names = FALSE)))
}

compare_one <- function(path) {
  cat(sprintf("\n=== %s ===\n", path))
  em <- tryCatch(fromJSON(path, simplifyVector = FALSE),
                 error = function(e) { cat("  read error: ", conditionMessage(e), "\n"); NULL })
  if (is.null(em)) return(invisible(NULL))

  if (is.null(em$trans) || length(em$trans) == 0L) {
    cat("  no `trans` block in JSON — skipping (likely converted from R RDS)\n")
    return(invisible(NULL))
  }

  trans    <- list_to_mat(em$trans)
  err_orig <- list_to_mat(em$err_out)

  errfun  <- if (!is.null(em$params$errfun))         em$params$errfun         else "loess"
  surface <- if (!is.null(em$params$loess$surface))  em$params$loess$surface  else "direct"

  cat(sprintf("  nq = %d, errfun (per JSON) = %s, surface = %s\n",
              ncol(trans), errfun, surface))
  if (!is.null(em$iterations)) cat(sprintf("  iterations = %d\n", em$iterations))
  if (!is.null(em$converged))  cat(sprintf("  converged = %s\n", em$converged))

  err_R <- switch(errfun,
    "loess"  = loessErrfun(trans, surface),
    "pacbio" = pacbioErrfun(trans, surface),
    {
      cat(sprintf("  errfun = %s not supported for comparison — skipping\n", errfun))
      return(invisible(NULL))
    }
  )

  d <- abs(err_R - err_orig)

  # For pacbio: q<93 is fit by loess (interesting), q=93 by Laplace MLE
  # (should match to machine epsilon).  Summarize separately so the loess
  # signal isn't masked by the MLE column (or vice versa).
  has_q93 <- errfun == "pacbio" && "93" %in% colnames(trans)
  if (has_q93) {
    q93_col <- which(colnames(trans) == "93")
    d_loess <- d[, -q93_col, drop = FALSE]
    d_mle   <- d[,  q93_col, drop = FALSE]
    summarize_diff(d_loess, err_orig[, -q93_col, drop = FALSE], "q<93, loess columns")
    summarize_diff(d_mle,   err_orig[,  q93_col, drop = FALSE], "q=93, Laplace MLE column")
  } else {
    summarize_diff(d, err_orig, "all cells")
  }

  cat(sprintf("\n  largest diffs (top 10 cells):\n"))
  ord <- order(d, decreasing = TRUE)[1:10]
  for (k in ord) {
    r <- ((k - 1L) %% 16L) + 1L
    c_ <- ((k - 1L) %/% 16L) + 1L
    rel <- if (err_orig[r, c_] > 0) 100 * d[r, c_] / err_orig[r, c_] else NA_real_
    tag <- if (has_q93 && c_ == q93_col) " [MLE]" else ""
    cat(sprintf("    %s  q=%2d%s:  orig=%.3e  R=%.3e  diff=%.3e  (%s%%)\n",
                ROW_NAMES[r], c_ - 1L, tag,
                err_orig[r, c_], err_R[r, c_], d[r, c_],
                if (is.na(rel)) "n/a" else sprintf("%.1f", rel)))
  }

  cat(sprintf("\n  per-transition max |diff|:\n"))
  for (i in seq_len(16)) {
    cat(sprintf("    %s : %.3e\n", ROW_NAMES[i], max(d[i, ])))
  }
  invisible(NULL)
}

for (path in args) compare_one(path)
