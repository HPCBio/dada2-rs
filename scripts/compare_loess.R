#!/usr/bin/env Rscript
# compare_loess.R
#
# Isolate the loess fit difference between dada2-rs's native Rust loess and
# R's stock stats::loess on a single, identical trans matrix.
#
# Reads one or more `learn-errors` JSON files. From each file pulls out the
# accumulated `trans` matrix and the embedded `err_out` (= whatever errfun
# produced this JSON), then re-fits R's stock loessErrfun on the same trans
# and reports per-cell diffs.
#
# Usage:
#   Rscript scripts/compare_loess.R <learn_errors.json> [more.json ...]
#
# Interpretation:
#   - If the JSON was produced with native Rust loess (default learn-errors),
#     the diff reveals the per-cell discrepancy between Rust and R loess.
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

list_to_mat <- function(x) {
  m <- do.call(rbind, lapply(x, as.numeric))
  rownames(m) <- ROW_NAMES
  colnames(m) <- as.character(seq_len(ncol(m)) - 1L)
  m
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

  errfun <- if (!is.null(em$params$errfun)) em$params$errfun else "<unknown>"
  cat(sprintf("  nq = %d, errfun (per JSON) = %s\n", ncol(trans), errfun))
  if (!is.null(em$iterations)) cat(sprintf("  iterations = %d\n", em$iterations))
  if (!is.null(em$converged)) cat(sprintf("  converged = %s\n", em$converged))

  err_R <- loessErrfun(trans, "direct")
  d <- abs(err_R - err_orig)

  cat(sprintf("\n  per-cell |diff|:\n"))
  cat(sprintf("    max     = %.3e\n", max(d)))
  cat(sprintf("    mean    = %.3e\n", mean(d)))
  cat(sprintf("    median  = %.3e\n", median(d)))
  cat(sprintf("    p99     = %.3e\n", quantile(d, 0.99, names = FALSE)))

  cat(sprintf("\n  largest diffs (top 10 cells):\n"))
  ord <- order(d, decreasing = TRUE)[1:10]
  for (k in ord) {
    r <- ((k - 1L) %% 16L) + 1L
    c_ <- ((k - 1L) %/% 16L) + 1L
    rel <- if (err_orig[r, c_] > 0) 100 * d[r, c_] / err_orig[r, c_] else NA_real_
    cat(sprintf("    %s  q=%2d:  orig=%.3e  R=%.3e  diff=%.3e  (%s%%)\n",
                ROW_NAMES[r], c_ - 1L,
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
