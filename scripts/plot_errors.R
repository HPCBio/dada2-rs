#!/usr/bin/env Rscript
#
# plot_errors.R — visualise the output of `dada2-rs learn-errors`
#
# Produces a 4×4 panel plot (one panel per nucleotide transition) that mirrors
# the style of DADA2's plotErrors():
#   • Points  — observed error rate from transition counts (sized by count)
#   • Black line  — estimated error rate (err_out, the model fit to trans)
#   • Dashed line — nominal error rate (input to the final DADA run, err_in)
#   • Red line    — theoretical Phred rate: (1/3) × 10^(-Q/10)
#
# Usage:
#   Rscript plot_errors.R <learn_errors.json> [output.pdf]
#
# If the output path is omitted the plot is written to
#   <input_stem>_errors.pdf  in the current directory.
#
# Dependencies: jsonlite, ggplot2

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_errors.R <learn_errors.json> [output.pdf]\n")
  quit(status = 1)
}

json_path  <- args[1]
out_path   <- if (length(args) >= 2) args[2] else {
  stem <- sub("\\.json$", "", basename(json_path))
  file.path(getwd(), paste0(stem, "_errors.pdf"))
}

# ---------------------------------------------------------------------------
# Load JSON
# ---------------------------------------------------------------------------
dat <- fromJSON(json_path)

nq        <- dat$nq
trans_mat <- dat$trans    # 16 × nq (list of 16 vectors from JSON)
err_in    <- dat$err_in   # 16 × nq
err_out   <- dat$err_out  # 16 × nq

# fromJSON reads 2-D arrays as matrices with rows = 16 and cols = nq.
# Ensure we have plain matrices.
trans_mat <- matrix(unlist(trans_mat), nrow = 16, ncol = nq, byrow = FALSE)
err_in    <- matrix(unlist(err_in),    nrow = 16, ncol = nq, byrow = FALSE)
err_out   <- matrix(unlist(err_out),   nrow = 16, ncol = nq, byrow = FALSE)

quality   <- seq(0, nq - 1)   # quality score values 0 .. nq-1

# ---------------------------------------------------------------------------
# Transition labels and indices
# ---------------------------------------------------------------------------
nts <- c("A", "C", "G", "T")

# Row r (0-indexed) = (nti * 4 + ntj), nti = ref, ntj = query.
# Build a data frame with one row per (quality, transition).
rows <- list()
for (nti in seq_along(nts)) {
  for (ntj in seq_along(nts)) {
    r       <- (nti - 1) * 4 + (ntj - 1) + 1  # 1-indexed row in matrix
    label   <- paste0(nts[nti], "\u2192", nts[ntj])
    is_self <- nti == ntj

    counts  <- trans_mat[r, ]
    # Total reads originating from reference nucleotide nti at each quality.
    row_total <- colSums(trans_mat[((nti - 1) * 4 + 1):((nti - 1) * 4 + 4), ,
                                   drop = FALSE])
    obs_rate  <- ifelse(row_total > 0, counts / row_total, NA_real_)

    # Theoretical Phred rate: errors = (1/3)*10^(-Q/10), self = 1-10^(-Q/10)
    theoretical <- if (is_self) {
      1 - 10 ^ (-quality / 10)
    } else {
      (1 / 3) * 10 ^ (-quality / 10)
    }

    rows[[length(rows) + 1]] <- data.frame(
      Quality     = quality,
      Transition  = label,
      is_self     = is_self,
      count       = counts,
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

# Fix factor order: A→A, A→C, … T→T (row-major, ref varies slowest)
trans_levels <- paste0(
  rep(nts, each = 4), "\u2192", rep(nts, times = 4)
)
df$Transition <- factor(df$Transition, levels = trans_levels)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
p <- ggplot(df, aes(x = Quality)) +
  # Observed points (omitted for self-transitions)
  geom_point(
    aes(y = Observed, size = count),
    colour  = "gray40",
    na.rm   = TRUE
  ) +
  # Theoretical Phred rate (red)
  geom_line(
    aes(y = Theoretical),
    colour    = "red2",
    linewidth = 0.6,
    na.rm     = TRUE
  ) +
  # Nominal / err_in (dashed, dark grey)
  geom_line(
    aes(y = Nominal),
    colour    = "gray30",
    linewidth = 0.6,
    linetype  = "dashed",
    na.rm     = TRUE
  ) +
  # Estimated / err_out (solid black)
  geom_line(
    aes(y = Estimated),
    colour    = "black",
    linewidth = 0.8,
    na.rm     = TRUE
  ) +
  scale_y_log10(
    labels = scales::label_log()
  ) +
  scale_size_area(max_size = 3, guide = "none") +
  facet_wrap(~Transition, nrow = 4, ncol = 4) +
  labs(
    x     = "Consensus quality score",
    y     = "Error frequency",
    title = "Error rates by transition type",
    caption = paste0(
      "Points: observed  \u2022  solid: estimated (err_out)  ",
      "\u2022  dashed: nominal (err_in)  \u2022  red: theoretical Phred"
    )
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
ggsave(out_path, plot = p, width = 10, height = 9)
message("Plot written to: ", out_path)
