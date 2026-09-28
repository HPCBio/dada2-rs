#!/usr/bin/env Rscript
# learnerrors_to_dada2rs.R
#
# Convert R DADA2 learnErrors() output (.rds) → dada2-rs JSON error model
# consumable by `dada2-rs dada` / `dada2-rs dada-pooled` via --error-model.
#
# Usage:
#   Rscript scripts/learnerrors_to_dada2rs.R <input.rds> <output.json>
#   Rscript scripts/learnerrors_to_dada2rs.R --help
#
# The input RDS may contain either:
#   (a) the list returned by learnErrors() — $err_out is used; or
#   (b) a 16-row error-rate matrix directly (e.g. saveRDS(getErrors(errR), ...)).
#
# The output JSON sets both `err_in` and `err_out` to the same matrix, so
# the value of dada2-rs's --use-err-in flag has no effect on downstream
# inference.  Row order must be A2A,A2C,A2G,A2T,C2A,...,T2T.
#
# Dependencies: jsonlite, optparse

suppressPackageStartupMessages({
  library(jsonlite)
  library(optparse)
})

parser <- OptionParser(
  usage = "usage: %prog <input.rds> <output.json>",
  description = paste0(
    "\nConvert an R DADA2 learnErrors() result (or a bare 16-row error\n",
    "matrix) saved as .rds into a dada2-rs --error-model JSON.\n")
)
args     <- parse_args(parser, positional_arguments = 2)$args
in_path  <- args[[1]]
out_path <- args[[2]]

obj <- readRDS(in_path)
err <- if (is.list(obj) && !is.null(obj$err_out)) {
  obj$err_out
} else if (is.matrix(obj)) {
  obj
} else {
  stop("Input RDS must be either a learnErrors() list (with $err_out) ",
       "or a 16-row matrix")
}

if (!is.matrix(err) || nrow(err) != 16L) {
  stop("Error matrix must have 16 rows; got ", nrow(err))
}
expected_rows <- paste0(rep(c("A", "C", "G", "T"), each = 4L), "2",
                        rep(c("A", "C", "G", "T"), times = 4L))
if (!is.null(rownames(err)) && !identical(rownames(err), expected_rows)) {
  stop("Row order must be ", paste(expected_rows, collapse = ","),
       "; got ", paste(rownames(err), collapse = ","))
}

nq <- ncol(err)
err_rows <- lapply(seq_len(16L), function(i) unname(as.numeric(err[i, ])))

output <- list(
  dada2_rs_command = "learn-errors",
  dada2_rs_version = "r-import",
  nq               = nq,
  err_in           = err_rows,
  err_out          = err_rows
)

writeLines(
  toJSON(output, auto_unbox = TRUE, digits = NA, pretty = TRUE),
  out_path
)
cat(sprintf("Wrote %s (16 x %d error matrix)\n", out_path, nq))
