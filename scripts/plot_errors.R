#!/usr/bin/env Rscript
#
# plot_errors.R — visualise the output of `dada2-rs learn-errors`
#
# Produces a 4x4 panel plot (one panel per nucleotide transition) in the style
# of DADA2's plotErrors() by Benjamin Callahan (see dada2/R/plot-methods.R in
# the DADA2 R package); all credit for the method is his.
#
#   • Points      — observed error rate from the transition counts
#   • Black line  — estimated error rate (err_out, the model fit to trans)
#   • Dashed line — nominal error rate fed into the final DADA run (err_in)
#   • Red line    — theoretical Phred rate: (1/3) x 10^(-Q/10)
#
# Beyond plotErrors(), the points can be weighted by *observation mass*: how
# much of the data actually sits at each quality score. This matters for binned
# quality scores (NovaSeq, Revio, i100), where nearly all of the mass sits on a
# handful of Q values and the visually prominent low-mass points carry almost
# no weight in the fit.
#
# Usage:
#   Rscript plot_errors.R [options] <learn_errors.json> [output.pdf]
#   Rscript plot_errors.R --help
#
# Dependencies: jsonlite, ggplot2, scales, optparse

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(optparse)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
nts <- c("A", "C", "G", "T")

opt_list <- list(
  make_option("--out", type = "character", default = NULL, metavar = "FILE",
              help = "Output path [default: <input_stem>_errors.pdf]"),
  make_option("--width", type = "double", default = 10, metavar = "N",
              help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 9, metavar = "N",
              help = "Figure height in inches [default %default]"),
  make_option("--title", type = "character",
              default = "Error rates by transition type", metavar = "TEXT",
              help = "Plot title [default \"%default\"]"),
  make_option("--nti", type = "character", default = "A,C,G,T",
              metavar = "A,C,G,T",
              help = "Reference (from) nucleotides to show [default %default]"),
  make_option("--ntj", type = "character", default = "A,C,G,T",
              metavar = "A,C,G,T",
              help = "Query (to) nucleotides to show [default %default]"),
  make_option("--mass", type = "character", default = "both", metavar = "MODE",
              help = paste(
                "How observation mass is shown on the points [default %default]:",
                "size = point area scales with the transition count;",
                "colour = point colour is the share of all observations from",
                "the reference base sitting at that quality score;",
                "both; none")),
  make_option("--mass-shading", action = "store_true", default = FALSE,
              dest = "mass_shading",
              help = "Shade each panel background by observation mass per Q"),
  make_option("--min-mass", type = "double", default = 0, metavar = "FRAC",
              dest = "min_mass",
              help = paste("Drop points at quality scores holding less than",
                           "FRAC of the observation mass [default %default]")),
  make_option("--size-trans", type = "character", default = "log10",
              metavar = "T", dest = "size_trans",
              help = paste("Transform for the count->point-size scale:",
                           "log10 | sqrt | identity [default %default].",
                           "log10 keeps rare-quality points visible when one",
                           "binned Q dominates the counts")),
  make_option("--no-obs", action = "store_false", default = TRUE,
              dest = "obs", help = "Omit the observed points"),
  make_option("--no-estimated", action = "store_false", default = TRUE,
              dest = "estimated", help = "Omit the estimated (err_out) line"),
  make_option("--no-nominal", action = "store_false", default = TRUE,
              dest = "nominal", help = "Omit the input (err_in) dashed line"),
  make_option("--no-theoretical", action = "store_false", default = TRUE,
              dest = "theoretical",
              help = "Omit the theoretical Phred (red) line")
)

parser <- OptionParser(
  usage = "usage: %prog [options] <learn_errors.json> [output.pdf]",
  option_list = opt_list,
  description = paste0(
    "\nPlot the error profile from `dada2-rs learn-errors`, in the style of\n",
    "DADA2's plotErrors(), with optional weighting of the observed points by\n",
    "how much of the data actually sits at each quality score.\n")
)
argv <- parse_args(parser, positional_arguments = c(1, 2))
opt  <- argv$options
positional <- argv$args

parse_nts <- function(x, flag) {
  v <- toupper(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
  v <- unique(v[nzchar(v)])
  if (length(v) == 0 || length(setdiff(v, nts)) > 0) {
    stop(sprintf("%s: must be a comma-separated subset of A,C,G,T", flag))
  }
  v
}

nti_sel <- parse_nts(opt$nti, "--nti")
ntj_sel <- parse_nts(opt$ntj, "--ntj")

mass_mode <- tolower(opt$mass)
if (mass_mode == "color") mass_mode <- "colour"
if (!mass_mode %in% c("both", "size", "colour", "none")) {
  stop("--mass must be one of: both, size, colour, none")
}

if (!opt$size_trans %in% c("log10", "sqrt", "identity")) {
  stop("--size-trans must be one of: log10, sqrt, identity")
}

min_mass <- opt$min_mass
if (is.na(min_mass) || min_mass < 0 || min_mass >= 1) {
  stop("--min-mass must be in [0, 1)")
}

width        <- opt$width
height       <- opt$height
title        <- opt$title
mass_shading <- opt$mass_shading
show_obs     <- opt$obs
show_est     <- opt$estimated
show_nom     <- opt$nominal
show_theo    <- opt$theoretical

json_path <- positional[1]
out_path  <- opt$out
if (is.null(out_path)) {
  out_path <- if (length(positional) >= 2) positional[2] else {
    stem <- sub("\\.json$", "", basename(json_path))
    file.path(getwd(), paste0(stem, "_errors.pdf"))
  }
}

# ---------------------------------------------------------------------------
# Load JSON
# ---------------------------------------------------------------------------
dat <- fromJSON(json_path)
if (!is.null(dat$type) && dat$type == "learn_errors" && !is.null(dat$data)) {
  dat <- dat$data
}

nq        <- dat$nq
trans_mat <- dat$trans    # 16 x nq
err_in    <- dat$err_in   # 16 x nq
err_out   <- dat$err_out  # 16 x nq

# fromJSON reads 2-D arrays as matrices with rows = 16 and cols = nq.
trans_mat <- matrix(unlist(trans_mat), nrow = 16, ncol = nq, byrow = FALSE)
err_in    <- matrix(unlist(err_in),    nrow = 16, ncol = nq, byrow = FALSE)
err_out   <- matrix(unlist(err_out),   nrow = 16, ncol = nq, byrow = FALSE)

quality <- seq(0, nq - 1)   # quality score values 0 .. nq-1

# ---------------------------------------------------------------------------
# Build the long-format table
#
# mass = share of all observations from reference base nti that sit at this
# quality score. Sums to 1 within each reference base, so it is comparable
# across panels in a row.
# ---------------------------------------------------------------------------
rows <- list()
for (nti in seq_along(nts)) {
  # Total observations originating from reference nucleotide nti per quality.
  row_total <- colSums(trans_mat[((nti - 1) * 4 + 1):((nti - 1) * 4 + 4), ,
                                 drop = FALSE])
  grand <- sum(row_total)
  mass  <- if (grand > 0) row_total / grand else rep(NA_real_, nq)

  for (ntj in seq_along(nts)) {
    r       <- (nti - 1) * 4 + (ntj - 1) + 1  # 1-indexed row in matrix
    label   <- paste0(nts[nti], "→", nts[ntj])
    is_self <- nti == ntj

    counts   <- trans_mat[r, ]
    obs_rate <- ifelse(row_total > 0, counts / row_total, NA_real_)

    # Theoretical Phred rate: errors = (1/3)*10^(-Q/10), self = 1-10^(-Q/10)
    theoretical <- if (is_self) {
      1 - 10 ^ (-quality / 10)
    } else {
      (1 / 3) * 10 ^ (-quality / 10)
    }

    rows[[length(rows) + 1]] <- data.frame(
      Quality     = quality,
      Transition  = label,
      from        = nts[nti],
      to          = nts[ntj],
      is_self     = is_self,
      count       = counts,
      total       = row_total,
      mass        = mass,
      Observed    = obs_rate,
      Estimated   = err_out[r, ],
      Nominal     = err_in[r, ],
      Theoretical = theoretical,
      stringsAsFactors = FALSE
    )
  }
}
df <- do.call(rbind, rows)

# Self-transitions (A→A etc.) are not errors; blank observed, estimated, and
# nominal so only the theoretical line shows — matching DADA2's plotErrors.
df$Observed[df$is_self]  <- NA_real_
df$Estimated[df$is_self] <- NA_real_
df$Nominal[df$is_self]   <- NA_real_

# Zero counts have no point to draw and break a log-transformed size scale.
df$count[!is.na(df$count) & df$count <= 0] <- NA_real_

# Zero rates cannot be drawn on a log axis.
for (col in c("Observed", "Estimated", "Nominal", "Theoretical")) {
  df[[col]][!is.na(df[[col]]) & df[[col]] <= 0] <- NA_real_
}

# Drop points at quality scores holding a negligible share of the data.
if (min_mass > 0) {
  df$Observed[!is.na(df$mass) & df$mass < min_mass] <- NA_real_
}

# Fix factor order: A→A, A→C, … T→T (row-major, ref varies slowest)
trans_levels <- paste0(rep(nts, each = 4), "→", rep(nts, times = 4))
df$Transition <- factor(df$Transition, levels = trans_levels)

df <- df[df$from %in% nti_sel & df$to %in% ntj_sel, ]
df$Transition <- droplevels(df$Transition)
if (nrow(df) == 0) stop("no transitions selected")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(df, aes(x = Quality))

# Background shading by observation mass (drawn first, under everything).
if (mass_shading) {
  shade <- df[!is.na(df$mass) & df$mass > 0, ]
  # One rect per (Transition, Q); drawn on a sqrt scale so that quality scores
  # carrying a small but non-trivial share stay visible next to a dominant bin.
  shade <- unique(shade[, c("Transition", "Quality", "mass")])
  mx <- max(shade$mass)
  # Finite y bounds: -Inf/Inf are NaN under the log10 transform.
  yv <- unlist(df[, c("Observed", "Estimated", "Nominal", "Theoretical")])
  yv <- yv[is.finite(yv) & yv > 0]
  shade$y_lo <- min(yv)
  shade$y_hi <- max(yv)
  p <- p + geom_rect(
    data = shade,
    aes(xmin = Quality - 0.5, xmax = Quality + 0.5,
        ymin = y_lo, ymax = y_hi, alpha = sqrt(mass / mx)),
    inherit.aes = FALSE,
    fill = "steelblue3"
  ) + scale_alpha_continuous(range = c(0.02, 0.45), limits = c(0, 1),
                             guide = "none")
}

if (show_obs) {
  pt_aes <- switch(
    mass_mode,
    both   = aes(y = Observed, size = count, colour = mass),
    size   = aes(y = Observed, size = count),
    colour = aes(y = Observed, colour = mass),
    none   = aes(y = Observed)
  )
  p <- p + if (mass_mode %in% c("both", "colour")) {
    geom_point(pt_aes, na.rm = TRUE)
  } else {
    geom_point(pt_aes, na.rm = TRUE, colour = "gray40")
  }
}

if (show_theo) {
  p <- p + geom_line(aes(y = Theoretical), colour = "red2",
                     linewidth = 0.6, na.rm = TRUE)
}
if (show_nom) {
  p <- p + geom_line(aes(y = Nominal), colour = "gray30", linewidth = 0.6,
                     linetype = "dashed", na.rm = TRUE)
}
if (show_est) {
  p <- p + geom_line(aes(y = Estimated), colour = "black",
                     linewidth = 0.8, na.rm = TRUE)
}

p <- p + scale_y_log10(labels = scales::label_log())

if (show_obs && mass_mode %in% c("both", "size")) {
  # plotErrors() uses scale_size_area (area proportional to count). With binned
  # quality scores a single Q holds most of the counts, which shrinks every
  # other point to a dot — so transform the size scale by default.
  p <- p + if (opt$size_trans == "identity") {
    scale_size_area(max_size = 3.5, guide = "none")
  } else {
    scale_size(range = c(0.5, 3.5), trans = opt$size_trans, guide = "none")
  }
}
if (show_obs && mass_mode %in% c("both", "colour")) {
  p <- p + scale_colour_viridis_c(
    name   = "Share of\nobservations",
    trans  = "sqrt",
    labels = scales::label_percent(accuracy = 1)
  )
}

caption <- paste0(
  if (show_obs) "points: observed  •  " else "",
  if (show_est) "solid: estimated (err_out)  •  " else "",
  if (show_nom) "dashed: nominal (err_in)  •  " else "",
  if (show_theo) "red: theoretical Phred  •  " else ""
)
caption <- sub("  •  $", "", caption)
if (show_obs && mass_mode %in% c("both", "colour")) {
  caption <- paste0(
    caption,
    "\ncolour: fraction of all observations from the reference base at that Q"
  )
}
if (min_mass > 0) {
  caption <- paste0(caption, sprintf(
    "\nquality scores holding < %s of the mass are not plotted",
    scales::label_percent(accuracy = 0.01)(min_mass)
  ))
}

p <- p +
  facet_wrap(~Transition, nrow = length(nti_sel)) +
  labs(
    x       = "Consensus quality score",
    y       = "Error frequency",
    title   = title,
    caption = caption
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
# cairo_pdf renders the arrow/bullet glyphs the base pdf() device substitutes.
dev <- if (grepl("\\.pdf$", out_path, ignore.case = TRUE) &&
           isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else NULL
ggsave(out_path, plot = p, width = width, height = height, device = dev)
message("Plot written to: ", out_path)
