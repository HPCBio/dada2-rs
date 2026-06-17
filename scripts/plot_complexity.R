#!/usr/bin/env Rscript
#
# Reproduce DADA2's plotComplexity() figure from one or more `summary` JSON
# outputs produced by `dada2-rs summary --complexity`. Mirrors the method and
# styling of dada2::plotComplexity: a histogram, per file, of each read's
# "effective oligonucleotide number" (exp(Shannon entropy) over its k-mer
# counts, range [1, 4^kmerSize]).
#
# The complexity statistic is a direct port of seqComplexity()/plotComplexity()
# from the DADA2 R package by Benjamin Callahan (original author); see
# dada2/R/filter.R and dada2/R/plot-methods.R. The per-read effective k-mer
# count computed by `summary --complexity` is bit-identical to R's
# dada2:::seqComplexity() on the fixtures. All credit for the method is his.
#
# dada2-rs emits a pre-binned histogram (it streams all reads rather than
# subsampling n, as R's FastqSampler does), so this script renders the bins as
# columns rather than recomputing geom_histogram() over raw values.
#
# Usage:
#   Rscript plot_complexity.R [--aggregate] [--out=plot.pdf] \
#                             [--width=8] [--height=5] \
#                             summary1.json [summary2.json ...]
#
# Defaults to writing complexity_profile.pdf in the current directory. Each
# input must have been produced with `summary --complexity`.

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

# ---- Argument parsing ---------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("usage: plot_complexity.R [--aggregate] [--out=FILE] [--width=N] [--height=N] summary.json [...]")
}

aggregate <- FALSE
out_file <- "complexity_profile.pdf"
width <- 8
height <- 5
files <- character(0)

for (a in args) {
  if (identical(a, "--aggregate")) {
    aggregate <- TRUE
  } else if (startsWith(a, "--out=")) {
    out_file <- sub("^--out=", "", a)
  } else if (startsWith(a, "--width=")) {
    width <- as.numeric(sub("^--width=", "", a))
  } else if (startsWith(a, "--height=")) {
    height <- as.numeric(sub("^--height=", "", a))
  } else if (startsWith(a, "--")) {
    stop(sprintf("unknown option: %s", a))
  } else {
    files <- c(files, a)
  }
}

if (length(files) == 0) stop("at least one summary JSON file is required")

# ---- Load complexity histograms -----------------------------------------

read_complexity <- function(path) {
  doc <- fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  if (!is.null(doc$data)) doc <- doc$data
  cx <- doc$complexity
  if (is.null(cx)) {
    stop(sprintf(
      "%s: no 'complexity' field — re-run `dada2-rs summary --complexity`",
      path
    ))
  }
  kmer_size <- as.integer(cx$kmer_size)
  bins      <- as.integer(cx$bins)
  counts    <- as.numeric(cx$histogram)
  max_c     <- 4^kmer_size
  # Bin b (0-based) covers [b, b+1) * max_c / bins; plot at its midpoint.
  midpoints <- (seq_len(bins) - 0.5) * (max_c / bins)
  data.frame(
    complexity = midpoints,
    Count      = counts,
    kmer_size  = kmer_size,
    max_c      = max_c,
    file       = basename(path)
  )
}

dfs <- lapply(files, read_complexity)
df  <- do.call(rbind, dfs)

kmer_size <- df$kmer_size[1]
max_c     <- df$max_c[1]
if (any(df$kmer_size != kmer_size)) {
  stop("inputs use different --complexity-kmer-size values; cannot combine on one axis")
}

if (aggregate) {
  agg <- aggregate(Count ~ complexity, df, sum)
  agg$file <- sprintf("%d files (aggregated)", length(files))
  df <- agg
}

# ---- Build plot ---------------------------------------------------------
# geom_col over bin midpoints reproduces plotComplexity's histogram shape from
# the pre-binned counts (binwidth = max_c / bins).

bins_n   <- nrow(dfs[[1]])
binwidth <- max_c / bins_n

p <- ggplot(df, aes(x = complexity, y = Count)) +
  geom_col(width = binwidth, fill = "grey35", na.rm = TRUE) +
  ylab("Count") + xlab("Effective Oligonucleotide Number") +
  theme_bw() +
  facet_wrap(~ file) +
  scale_x_continuous(
    limits = c(0, max_c),
    breaks = seq(0, max_c, max_c / 4)
  )

ggsave(out_file, plot = p, width = width, height = height)
message(sprintf("wrote %s", out_file))
