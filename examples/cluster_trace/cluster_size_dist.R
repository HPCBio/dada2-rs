#!/usr/bin/env Rscript
# cluster_size_dist.R
#
# Plot the distribution of cluster sizes for one or more cluster-trace
# JSON files (output of `--cluster-trace-dir` or `--cluster-trace`).
#
# Usage:
#   Rscript cluster_size_dist.R out.pdf trace_dir/cluster_iter_*.json
#   Rscript cluster_size_dist.R out.pdf clusters.json
#
# Produces a multi-panel PDF: one panel per file with cluster abundance on
# log-y, ordered by rank, colored by birth type. Useful for spotting
# convergence (panels stop changing) and outlier clusters at low abundance.

library(jsonlite)
library(ggplot2)

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: cluster_size_dist.R <out.pdf> <trace.json>...")
outpdf <- args[1]
files  <- args[-1]

read_one <- function(path) {
  j <- jsonlite::read_json(path, simplifyVector = TRUE)
  if (j$dada2_rs_command != "cluster-trace") {
    stop(sprintf("%s: not a cluster-trace file (tag=%s)", path, j$dada2_rs_command))
  }
  data.frame(
    file        = basename(path),
    sample      = j$sample,
    iteration   = if (is.null(j$iteration)) NA_integer_ else j$iteration,
    cluster_id  = j$clusters$id,
    abundance   = j$clusters$abundance,
    n_members   = j$clusters$n_members,
    birth_type  = j$clusters$birth_type,
    birth_hamming = j$clusters$birth_hamming
  )
}

df <- do.call(rbind, lapply(files, read_one))
df <- df[order(df$file, -df$abundance), ]
df$rank <- ave(df$abundance, df$file, FUN = function(x) seq_along(x))

p <- ggplot(df, aes(x = rank, y = abundance, color = birth_type, shape = birth_type)) +
  geom_point(alpha = 0.8) +
  scale_y_log10() +
  facet_wrap(~ file, scales = "free_x") +
  labs(x = "Cluster rank", y = "Reads (log10)",
       title = "Cluster size distribution",
       color = "Birth type", shape = "Birth type") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(outpdf, p, width = 10, height = 7)
cat(sprintf("Wrote %s (%d clusters across %d file(s))\n",
            outpdf, nrow(df), length(files)))
