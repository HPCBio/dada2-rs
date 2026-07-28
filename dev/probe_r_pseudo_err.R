#!/usr/bin/env Rscript
# probe_r_pseudo_err.R — does R DADA2's pool="pseudo" change the error model
# between its two rounds? (dada2-rs issue #100)
#
# The table-level evidence already says yes: enabling our
# `--reestimate-err-between-rounds` moves us from 31 ASVs apart from R to 3, and
# 30 of R's 31 extra ASVs are exactly the ones the re-fit produces. That is
# inference from output. This observes the mechanism directly.
#
# Method: trace dada2:::dada_uniques, the Rcpp entry point every per-sample
# inference goes through, and record the `err` matrix (4th formal) handed to it on
# every call. Then compare what each call received against the matrix the caller
# supplied. No dada2 source changes; nothing but a trace.
#
# Expectation if R re-fits (dada.R:371-378, no selfConsist guard):
#   pool="pseudo"  -> 2 turns x n samples calls; turn-1 calls get the supplied
#                     err, turn-2 calls get something else.
#   pool=FALSE     -> 1 turn; every call gets the supplied err. (Control: shows
#                     the difference is pseudo-specific, not just "dada refits".)
#
# Deliberately small: this is a binary property of the code path, not of the data,
# so two shallow samples answer it exactly as well as a full run would.
#
# Usage:
#   Rscript dev/probe_r_pseudo_err.R [sample1.fastq.gz sample2.fastq.gz]
# Defaults to the dada2 package fixtures shipped in data/dada2/.

suppressMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  fns <- args[seq_len(2)]
} else {
  here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  fns <- file.path(here, "..", "data", "dada2", c("sam1F.fastq.gz", "sam2F.fastq.gz"))
}
stopifnot(all(file.exists(fns)))
cat("samples:\n"); cat(paste0("  ", fns, collapse = "\n"), "\n\n")

drps <- lapply(fns, derepFastq, verbose = FALSE)

# One error model, learned once, then supplied to every run below. selfConsist is
# off, so a caller would reasonably expect this exact matrix to be used throughout.
cat("learning a reference error model...\n")
err0 <- learnErrors(drps, multithread = FALSE, verbose = FALSE)
cat("done:", nrow(err0$err_out), "x", ncol(err0$err_out), "\n\n")
ERR <- err0$err_out

# ---- trace dada_uniques and capture the `err` it is handed ------------------
seen <- new.env(parent = emptyenv())
seen$mats <- list()
trace(
  dada2:::dada_uniques,
  tracer = quote({
    e <- get("err", envir = environment())
    assign("mats", c(get("mats", envir = seen), list(e)), envir = seen)
  }),
  print = FALSE
)
on.exit(try(untrace(dada2:::dada_uniques), silent = TRUE), add = TRUE)

run_probe <- function(label, pool) {
  seen$mats <- list()
  invisible(capture.output(suppressMessages(
    dada(drps, err = ERR, pool = pool, multithread = FALSE, verbose = FALSE)
  )))
  mats <- seen$mats
  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  dada_uniques calls: %d  (%d sample(s))\n", length(mats), length(drps)))
  for (i in seq_along(mats)) {
    m <- mats[[i]]
    # Compare against the SUPPLIED matrix, trimming to the common column count:
    # dada() extends an error matrix when observed Q exceeds its columns.
    nc <- min(ncol(m), ncol(ERR))
    d <- max(abs(m[, seq_len(nc)] - ERR[, seq_len(nc)]))
    cat(sprintf("  call %d: %s supplied err   (max |delta| = %.3e)\n",
                i, if (d == 0) "==" else "!=", d))
  }
  changed <- vapply(mats, function(m) {
    nc <- min(ncol(m), ncol(ERR))
    max(abs(m[, seq_len(nc)] - ERR[, seq_len(nc)])) > 0
  }, logical(1))
  cat(sprintf("  -> %d of %d calls used a DIFFERENT error model\n\n",
              sum(changed), length(changed)))
  invisible(changed)
}

cat("Does every inference use the error model that was supplied?\n\n")
ctrl <- run_probe('pool=FALSE  (control: one round)', FALSE)
pseu <- run_probe('pool="pseudo"  (two rounds)', "pseudo")

untrace(dada2:::dada_uniques)

cat("========================================================================\n")
if (any(pseu) && !any(ctrl)) {
  cat("CONFIRMED: pool=\"pseudo\" denoises part of its work with an error model\n")
  cat("that is NOT the one supplied, while pool=FALSE uses the supplied model\n")
  cat("throughout. The re-fit is pseudo-specific -- consistent with dada.R:371-378\n")
  cat("refitting err at the end of each turn with no selfConsist guard.\n")
} else if (!any(pseu)) {
  cat("NOT CONFIRMED: every pool=\"pseudo\" call used the supplied error model.\n")
  cat("The table-level result in issue #100 then needs another explanation.\n")
} else {
  cat("AMBIGUOUS: pool=FALSE also deviated from the supplied model, so the\n")
  cat("deviation is not specific to pseudo-pooling. Investigate before reporting.\n")
}
cat("========================================================================\n")
