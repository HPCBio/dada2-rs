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
#   Rscript plot_complexity.R --help
#
# Defaults to writing complexity_profile.pdf in the current directory. Each
# input must have been produced with `summary --complexity`.
#
# Dependencies: jsonlite, ggplot2, optparse

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(optparse)
})

# ---- Argument parsing ---------------------------------------------------

opt_list <- list(
  make_option("--aggregate", action = "store_true", default = FALSE,
              help = "Pool all inputs into one panel instead of one per file"),
  make_option("--out", type = "character", default = "complexity_profile.pdf",
              metavar = "FILE", help = "Output path [default %default]"),
  make_option("--width", type = "double", default = 8, metavar = "N",
              help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 5, metavar = "N",
              help = "Figure height in inches [default %default]")
)

parser <- OptionParser(
  usage = "usage: %prog [options] summary1.json [summary2.json ...]",
  option_list = opt_list,
  description = paste0(
    "\nPlot read complexity from `dada2-rs summary --complexity` JSON, in the\n",
    "style of DADA2's plotComplexity().\n")
)
argv <- parse_args(parser, positional_arguments = c(1, Inf))
opt  <- argv$options

# optparse warns on a non-numeric double but passes the string through.
for (flag in c("width", "height")) {
  v <- suppressWarnings(as.numeric(opt[[flag]]))
  if (is.na(v) || v <= 0) stop(sprintf("--%s must be a positive number", flag))
  opt[[flag]] <- v
}

aggregate <- opt$aggregate
out_file  <- opt$out
width     <- opt$width
height    <- opt$height
files     <- argv$args

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
