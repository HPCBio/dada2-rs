#!/usr/bin/env Rscript
# lambda_distribution.R
#
# For each cluster in a trace file, plot the distribution of member
# log10(λ) values vs hamming distance to the center.  Reveals which
# members sit near the acceptance boundary and which are deep in the
# cluster's basin.
#
# Usage:
#   Rscript lambda_distribution.R out.pdf trace.json [trace.json ...]

library(jsonlite)
library(ggplot2)

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: lambda_distribution.R <out.pdf> <trace.json>...")
outpdf <- args[1]
files  <- args[-1]

read_members <- function(path) {
  j <- jsonlite::read_json(path, simplifyVector = TRUE)
  if (j$dada2_rs_command != "cluster-trace") stop(sprintf("%s: not a cluster trace", path))
  if (isTRUE(j$trace_no_members)) {
    warning(sprintf("%s was written with --trace-no-members; skipping", path))
    return(NULL)
  }
  rows <- list()
  for (i in seq_along(j$clusters$id)) {
    cid <- j$clusters$id[i]
    m   <- j$clusters$members[[i]]
    if (is.null(m) || nrow(m) == 0) next
    rows[[length(rows) + 1]] <- data.frame(
      file       = basename(path),
      cluster_id = cid,
      hamming    = m$hamming,
      lambda     = m$lambda,
      abundance  = m$abundance,
      pval       = m$pval
    )
  }
  do.call(rbind, rows)
}

df <- do.call(rbind, lapply(files, read_members))
if (is.null(df) || nrow(df) == 0) stop("no member rows found in any file")
df$log10_lambda <- log10(pmax(df$lambda, 1e-300))

p <- ggplot(df, aes(x = hamming, y = log10_lambda, size = log10(abundance))) +
  geom_jitter(alpha = 0.4, width = 0.15) +
  facet_wrap(~ file) +
  labs(x = "Hamming distance to center",
       y = "log10(λ)",
       size = "log10(abundance)",
       title = "Member λ vs hamming, by cluster") +
  theme_bw()

ggsave(outpdf, p, width = 10, height = 7)
cat(sprintf("Wrote %s (%d member rows)\n", outpdf, nrow(df)))
