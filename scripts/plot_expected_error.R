#!/usr/bin/env Rscript
#
# Plot per-position cumulative expected error (EE) from one or more `summary`
# JSON outputs produced by `dada2-rs summary --expected-error`. EE at a position
# is the read's cumulative Σ 10^(-Q/10) up to that base — the quantity
# `filter-and-trim` thresholds against `maxEE`. This shows how EE accumulates
# along the read so you can judge truncation length and maxEE choices: where the
# curves cross the maxEE reference lines (2/3/5/7), reads start being discarded.
#
# Inspired by Remi Maglione's Qual_vs_MaxEE plot
# (https://github.com/RemiMaglione/r-scripts), but driven by the *true* per-read
# cumulative-EE distribution that dada2-rs aggregates (exact mean/min/max plus
# quartiles), rather than EE derived from the per-position mean-quality curve.
#
# Usage:
#   Rscript plot_expected_error.R [--out=plot.pdf] [--linear] \
#                                 [--width=8] [--height=5] \
#                                 summary1.json [summary2.json ...]
#   Rscript plot_expected_error.R --help
#
# Defaults to writing expected_error.pdf with a log10 y-axis (use --linear for a
# linear y-axis). Each input must have been produced with `--expected-error`.
#
# Dependencies: jsonlite, ggplot2, optparse

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(optparse)
})

opt_list <- list(
  make_option("--out", type = "character", default = "expected_error.pdf",
              metavar = "FILE", help = "Output path [default %default]"),
  make_option("--linear", action = "store_true", default = FALSE,
              help = "Use a linear y-axis instead of log10"),
  make_option("--width", type = "double", default = 8, metavar = "N",
              help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 5, metavar = "N",
              help = "Figure height in inches [default %default]")
)

parser <- OptionParser(
  usage = "usage: %prog [options] summary1.json [summary2.json ...]",
  option_list = opt_list,
  description = paste0(
    "\nPlot per-position cumulative expected error from\n",
    "`dada2-rs summary --expected-error` JSON.\n")
)
argv <- parse_args(parser, positional_arguments = c(1, Inf))
opt  <- argv$options

# optparse warns on a non-numeric double but passes the string through.
for (flag in c("width", "height")) {
  v <- suppressWarnings(as.numeric(opt[[flag]]))
  if (is.na(v) || v <= 0) stop(sprintf("--%s must be a positive number", flag))
  opt[[flag]] <- v
}

out_file <- opt$out
log_y    <- !opt$linear
width    <- opt$width
height   <- opt$height
files    <- argv$args

read_ee <- function(path) {
  doc <- fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  if (!is.null(doc$data)) doc <- doc$data
  ee <- doc$expected_error
  if (is.null(ee)) {
    stop(sprintf(
      "%s: no 'expected_error' field — re-run `dada2-rs summary --expected-error`",
      path
    ))
  }
  data.frame(
    Position = seq_along(ee$mean),
    Mean     = as.numeric(ee$mean),
    Min      = as.numeric(ee$min),
    Q25      = as.numeric(ee$q25),
    Median   = as.numeric(ee$median),
    Q75      = as.numeric(ee$q75),
    Max      = as.numeric(ee$max),
    file     = basename(path)
  )
}

df <- do.call(rbind, lapply(files, read_ee))

# maxEE reference cutoffs, as used by filter-and-trim / DADA2 filterAndTrim.
cutoffs <- c(2, 3, 5, 7)

p <- ggplot(df, aes(x = Position)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = "#FC8D62", alpha = 0.25) +
  geom_line(aes(y = Median), color = "#FC8D62") +
  geom_line(aes(y = Mean), color = "#66C2A5") +
  geom_line(aes(y = Max), color = "grey50", linewidth = 0.25, linetype = "dashed") +
  geom_hline(yintercept = cutoffs, color = "red", linewidth = 0.25) +
  ylab("Cumulative expected error") + xlab("Position (cycle)") +
  theme_bw() +
  facet_wrap(~ file)

if (log_y) {
  p <- p + scale_y_log10()
}

ggsave(out_file, plot = p, width = width, height = height)
message(sprintf("wrote %s", out_file))
