#!/usr/bin/env Rscript
#
# plot_indel_errors.R — plot observed indel error rates (issue #291)
#
# The indel counterpart of plot_errors.R, from the read-weighted tables that
# `dada2-rs dada --indel-opportunities` writes (errors between each unique and
# its final cluster centre, `correct` members only). Three pages:
#
#   1. Insertions: rate per base by the inserted base's quality, one panel for
#      bases matching neither neighbour (novel) and one for bases that lengthen
#      a homopolymer run.
#   2. Deletions: rate per centre base by the lower quality of the read bases
#      either side, one panel per centre run length (1, 2, 3, 4, 5+).
#   3. Indel lengths: events by length, deletions negative and insertions
#      positive, as in variant-calling length spectra.
#
# Point area scales with the opportunities at that quality (log10), as
# plot_errors.R scales by observation mass, so points backed by little data are
# visibly small. The dashed red line is the Phred rate 10^(-Q/10), for scale
# only: indel rates have no nominal Phred value.
#
# Usage:
#   Rscript plot_indel_errors.R [options] <opportunities.tsv> [output.pdf]
#   Rscript plot_indel_errors.R --help
#
# Dependencies: ggplot2, scales, optparse

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

opt_list <- list(
  make_option("--out", type = "character", default = NULL, metavar = "FILE",
              help = "Output path [default: <input_stem>_indel_errors.pdf]"),
  make_option("--sample", type = "character", default = NULL, metavar = "NAME",
              help = "Plot one sample [default: all rows summed]"),
  make_option("--title", type = "character", default = NULL, metavar = "TEXT",
              help = "Title prefix [default: the input file name]"),
  make_option("--min-opp", type = "double", default = 100, metavar = "N",
              dest = "min_opp",
              help = "Drop points backed by fewer than N opportunities [default %default]"),
  make_option("--max-len", type = "integer", default = 20, metavar = "N",
              dest = "max_len",
              help = "Pool indel lengths above N into the last bar [default %default]"),
  make_option("--width", type = "double", default = 10, metavar = "N",
              help = "Page width in inches [default %default]"),
  make_option("--height", type = "double", default = 6, metavar = "N",
              help = "Page height in inches [default %default]"),
  make_option("--no-theoretical", action = "store_false", default = TRUE,
              dest = "theoretical", help = "Omit the Phred reference line")
)
parser <- OptionParser(
  usage = "usage: %prog [options] <opportunities.tsv> [output.pdf]",
  option_list = opt_list,
  description = paste0(
    "\nPlot observed indel error rates from `dada2-rs dada --indel-opportunities`:\n",
    "insertions by quality, deletions by run length and flanking quality, and\n",
    "the indel length distribution.\n")
)
argv <- parse_args(parser, positional_arguments = c(1, 2))
opt <- argv$options
tsv <- argv$args[1]
for (flag in c("width", "height", "min_opp")) {
  v <- suppressWarnings(as.numeric(opt[[flag]]))
  if (is.na(v) || v < 0) stop(sprintf("--%s must be a non-negative number", gsub("_", "-", flag)))
  opt[[flag]] <- v
}
out <- opt$out
if (is.null(out)) {
  out <- if (length(argv$args) >= 2) argv$args[2] else
    file.path(getwd(), paste0(sub("\\.tsv$", "", basename(tsv)), "_indel_errors.pdf"))
}
title <- if (is.null(opt$title)) basename(tsv) else opt$title

# ---------------------------------------------------------------------------
# Load: long table of (sample, table, key, opportunities, events)
# ---------------------------------------------------------------------------
d <- read.delim(tsv, comment.char = "#", stringsAsFactors = FALSE)
need <- c("sample", "table", "key", "opportunities", "events")
if (!all(need %in% names(d))) stop(tsv, ": not an --indel-opportunities table")
if (!is.null(opt$sample)) {
  d <- d[d$sample == opt$sample, ]
  if (nrow(d) == 0) stop("no rows for sample ", opt$sample)
}
d <- aggregate(cbind(opportunities, events) ~ table + key, data = d, FUN = sum)
have <- unique(d$table)
for (t in c("ins_q_novel", "ins_q_run", paste0("del_fq_L", 1:5), "ins_len", "del_len")) {
  if (!t %in% have) stop(tsv, ": no `", t, "` table; rerun dada with this dada2-rs version")
}
rates <- function(tables, labels) {
  x <- d[d$table %in% tables, ]
  x$panel <- factor(labels[match(x$table, tables)], levels = labels)
  x <- x[x$opportunities >= opt$min_opp, ]
  x$rate <- x$events / x$opportunities
  x
}
phred <- function(q) 10^(-q / 10)

rate_plot <- function(x, xlab, subtitle) {
  pos <- x[x$events > 0, ]
  zero <- x[x$events == 0, ]
  p <- ggplot(pos, aes(key, rate)) +
    geom_point(aes(size = opportunities), alpha = 0.6) +
    facet_wrap(~panel, nrow = 1) +
    scale_y_log10(labels = scales::label_math(10^.x, format = log10)) +
    scale_size_area(trans = "log10", max_size = 4, name = "opportunities") +
    labs(x = xlab, y = "error rate per base", title = title, subtitle = subtitle) +
    theme_bw()
  if (nrow(zero) > 0) {
    # A zero has no place on a log axis; mark the quality where none were seen.
    p <- p + geom_rug(data = zero, aes(x = key), inherit.aes = FALSE, sides = "b", alpha = 0.4)
  }
  if (opt$theoretical) {
    qs <- seq(min(x$key), max(x$key))
    p <- p + geom_line(data = data.frame(key = qs, rate = phred(qs)), aes(key, rate),
                       inherit.aes = FALSE, colour = "red", linetype = "dashed")
  }
  p
}

ins <- rates(c("ins_q_novel", "ins_q_run"), c("novel base", "lengthens a run"))
del <- rates(paste0("del_fq_L", 1:5), c("run 1", "run 2", "run 3", "run 4", "run 5+"))

lens <- d[d$table %in% c("ins_len", "del_len") & d$events > 0, ]
lens$length <- ifelse(lens$table == "del_len", -lens$key, lens$key)
lens$length <- pmax(pmin(lens$length, opt$max_len), -opt$max_len)
lens <- aggregate(events ~ length + table, data = lens, FUN = sum)
lens$type <- ifelse(lens$table == "del_len", "deletion", "insertion")

pdf(out, width = opt$width, height = opt$height)
print(rate_plot(ins, "quality of the inserted base",
                "Insertions; ticks on the axis mark qualities with no insertions"))
print(rate_plot(del, "lower quality of the read bases either side",
                "Deletions, by centre homopolymer run length"))
print(ggplot(lens, aes(length, events, fill = type)) +
        geom_col(width = 0.8) +
        scale_y_log10() +
        scale_fill_manual(values = c(deletion = "#b2182b", insertion = "#2166ac")) +
        labs(x = sprintf("indel length (deletions < 0; lengths beyond %d pooled)", opt$max_len),
             y = "read-weighted events", title = title, subtitle = "Indel length distribution") +
        theme_bw())
invisible(dev.off())
cat("wrote", out, "\n")
