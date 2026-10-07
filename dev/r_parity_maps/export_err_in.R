#!/usr/bin/env Rscript
# export_err_in.R — write each self-consistency pass's input model from a
# learnErrors object as a dada2-rs --error-model JSON (#277).
#
# learnErrors keeps $err_in as a list: element k is the model pass k ran under
# (dada.R stores it from pass 1 on; the max-rate initial pass is not kept).
# Feeding element k to both implementations replays pass k with identical input.
#
# Usage:
#   Rscript export_err_in.R <err.rds> <out-dir>
# Writes <out-dir>/err_in_<k>.json for every k, via
# scripts/learnerrors_to_dada2rs.R (17 significant digits, exact round trip).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: export_err_in.R <err.rds> <out-dir>")
e <- readRDS(args[1]); out_dir <- args[2]
if (!is.list(e$err_in)) stop(args[1], ": $err_in is not a per-pass list (not a selfConsist learnErrors object?)")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
conv <- file.path(script_dir, "..", "..", "scripts", "learnerrors_to_dada2rs.R")
for (k in seq_along(e$err_in)) {
  rds <- file.path(out_dir, sprintf("err_in_%d.rds", k))
  saveRDS(e$err_in[[k]], rds)
  status <- system2("Rscript", c(conv, rds, file.path(out_dir, sprintf("err_in_%d.json", k))))
  if (status != 0) stop("conversion failed for pass ", k)
}
cat(sprintf("%d pass(es); the last equals $err_out: %s\n", length(e$err_in),
            identical(unname(e$err_in[[length(e$err_in)]]), unname(e$err_out))))
