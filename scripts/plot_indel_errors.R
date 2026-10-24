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
# --fit overlays fitted rates (`learn-errors --indel-rates-out`, or an
# `--indel-fit-cmd` output) as lines, one colour per fit, on pages 1 and 2.
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
  make_option("--overlay", type = "character", default = NULL, metavar = "A,B",
              help = paste("Overlay several samples, coloured by sample, e.g.",
                           "iter0,iter3 from a learn-errors --indel-history table")),
  make_option("--fit", type = "character", default = NULL, metavar = "NAME=TSV,...",
              help = paste("Overlay fitted rates as lines, e.g.",
                           "pooled=rates.tsv,gam=gam_rates.tsv")),
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
overlay <- !is.null(opt$overlay)
if (overlay) {
  keep <- trimws(strsplit(opt$overlay, ",", fixed = TRUE)[[1]])
  missing <- setdiff(keep, unique(d$sample))
  if (length(missing) > 0) stop("no rows for sample(s): ", paste(missing, collapse = ", "))
  d <- d[d$sample %in% keep, ]
  d$sample <- factor(d$sample, levels = keep)
} else {
  d$sample <- "all"
}
d <- aggregate(cbind(opportunities, events) ~ sample + table + key, data = d, FUN = sum)
fits <- NULL
if (!is.null(opt$fit)) {
  if (overlay) stop("--fit and --overlay both use colour; pick one")
  for (spec in trimws(strsplit(opt$fit, ",", fixed = TRUE)[[1]])) {
    kv <- strsplit(spec, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("--fit expects NAME=TSV, got `", spec, "`")
    f <- read.delim(kv[2], comment.char = "#", stringsAsFactors = FALSE)
    if (!identical(names(f), c("table", "key", "rate"))) stop(kv[2], ": not a fitted-rates table")
    f$fit <- kv[1]
    fits <- rbind(fits, f)
  }
  fits$fit <- factor(fits$fit, levels = unique(fits$fit))
}
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

rate_plot <- function(x, xlab, subtitle, tables, labels) {
  pos <- x[x$events > 0, ]
  zero <- x[x$events == 0, ]
  p <- if (is.null(fits)) {
    ggplot(pos, aes(key, rate, colour = sample)) +
      geom_point(aes(size = opportunities), alpha = 0.6) +
      scale_colour_brewer(palette = "Dark2", guide = if (overlay) "legend" else "none")
  } else {
    # Fitted lines span the qualities shown for each panel.
    fl <- fits[fits$table %in% tables, ]
    fl$panel <- factor(labels[match(fl$table, tables)], levels = labels)
    rng <- aggregate(key ~ panel, data = x, FUN = range)
    lo <- rng$key[, 1][match(fl$panel, rng$panel)]
    hi <- rng$key[, 2][match(fl$panel, rng$panel)]
    fl <- fl[!is.na(lo) & fl$key >= lo & fl$key <= hi & fl$rate > 0, ]
    ggplot(pos, aes(key, rate)) +
      geom_point(aes(size = opportunities), colour = "grey40", alpha = 0.5) +
      geom_line(data = fl, aes(key, rate, colour = fit), linewidth = 0.8) +
      scale_colour_brewer(palette = "Dark2", name = "fit") +
      guides(colour = guide_legend(override.aes = list(linewidth = 1.2, size = 0)))
  }
  p <- p +
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

ins_t <- c("ins_q_novel", "ins_q_run"); ins_l <- c("novel base", "lengthens a run")
del_t <- paste0("del_fq_L", 1:5); del_l <- c("run 1", "run 2", "run 3", "run 4", "run 5+")
ins <- rates(ins_t, ins_l)
del <- rates(del_t, del_l)

lens <- d[d$table %in% c("ins_len", "del_len") & d$events > 0, ]
lens$length <- ifelse(lens$table == "del_len", -lens$key, lens$key)
lens$length <- pmax(pmin(lens$length, opt$max_len), -opt$max_len)
lens <- aggregate(events ~ sample + length + table, data = lens, FUN = sum)
lens$type <- ifelse(lens$table == "del_len", "deletion", "insertion")

pdf(out, width = opt$width, height = opt$height)
print(rate_plot(ins, "quality of the inserted base",
                "Insertions; ticks on the axis mark qualities with no insertions", ins_t, ins_l))
print(rate_plot(del, "lower quality of the read bases either side",
                "Deletions, by centre homopolymer run length", del_t, del_l))
# Stems from below 1 so single events stay visible on the log axis.
len_plot <- if (overlay) {
  # Dodge iterations side by side; colour by sample, shape by type.
  ggplot(lens, aes(length, events, colour = sample, shape = type)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
    scale_colour_brewer(palette = "Dark2")
} else {
  ggplot(lens, aes(length, events, colour = type)) +
    geom_segment(aes(xend = length, y = 0.5, yend = events), linewidth = 1.2) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(deletion = "#b2182b", insertion = "#2166ac"))
}
print(len_plot +
        scale_y_log10() +
        scale_x_continuous(breaks = function(l) {
          b <- scales::breaks_extended()(l)
          unique(round(b[b == round(b)]))
        }) +
        labs(x = sprintf("indel length (deletions < 0; lengths beyond %d pooled)", opt$max_len),
             y = "read-weighted events", title = title, subtitle = "Indel length distribution") +
        theme_bw())
invisible(dev.off())
cat("wrote", out, "\n")
