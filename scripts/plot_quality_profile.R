#!/usr/bin/env Rscript
#
# Reproduce DADA2's plotQualityProfile() figure from one or more `summary`
# JSON outputs produced by dada2-rs. Mirrors the styling of
# dada2::plotQualityProfile (heatmap of cycle x quality, plus per-cycle
# mean / Q25 / median / Q75 / cumulative-reads lines).
#
# Usage:
#   Rscript plot_quality_profile.R [--aggregate] [--out=plot.pdf] \
#                                  [--width=8] [--height=5] \
#                                  summary1.json [summary2.json ...]
#
# Defaults to writing quality_profile.pdf in the current directory.

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

# ---- Argument parsing ---------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("usage: plot_quality_profile.R [--aggregate] [--out=FILE] [--width=N] [--height=N] summary.json [...]")
}

aggregate <- FALSE
out_file <- "quality_profile.pdf"
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

# ---- Load summaries -----------------------------------------------------

read_summary <- function(path) {
  doc <- fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  if (!is.null(doc$type) && doc$type == "summary") {
    doc <- doc$data
  }
  hist <- doc$quality_histogram
  if (is.null(hist)) {
    stop(sprintf(
      "%s: no 'quality_histogram' field — re-run dada2-rs summary to regenerate",
      path
    ))
  }
  # quality_histogram is a list of length n_cycles, each element a numeric
  # vector of length max_quality + 1 (counts at qualities 0..max_quality).
  hist_mat <- do.call(rbind, lapply(hist, as.numeric))
  list(
    total_reads = as.numeric(doc$total_reads),
    reads_per_position = as.numeric(doc$reads_per_position),
    hist = hist_mat,
    file = basename(path)
  )
}

summaries <- lapply(files, read_summary)

# ---- Build long-form data frames ---------------------------------------

# plotdf: one row per (file, cycle, quality) tile with count > 0.
# statdf: one row per (file, cycle) with mean / Q25 / Q50 / Q75 / cum.
plotdf_list <- list()
statdf_list <- list()
anndf_list  <- list()

quantile_from_hist <- function(scores, counts, q) {
  total <- sum(counts)
  if (total == 0) return(NA_real_)
  scores[which(cumsum(counts) / total >= q)[1]]
}

for (s in summaries) {
  n_cycles <- nrow(s$hist)
  n_qual   <- ncol(s$hist) - 1L
  scores   <- 0:n_qual

  # Long-form heatmap rows.
  cycle_vec <- rep(seq_len(n_cycles), each = length(scores))
  score_vec <- rep(scores, times = n_cycles)
  count_vec <- as.numeric(t(s$hist))
  keep <- count_vec > 0
  plotdf_list[[length(plotdf_list) + 1L]] <- data.frame(
    Cycle = cycle_vec[keep],
    Score = score_vec[keep],
    Count = count_vec[keep],
    file  = s$file
  )

  # Per-cycle summary statistics.
  means <- numeric(n_cycles)
  q25   <- numeric(n_cycles)
  q50   <- numeric(n_cycles)
  q75   <- numeric(n_cycles)
  cums  <- numeric(n_cycles)
  for (i in seq_len(n_cycles)) {
    counts_i <- s$hist[i, ]
    cums[i]  <- sum(counts_i)
    if (cums[i] > 0) {
      means[i] <- sum(scores * counts_i) / cums[i]
      q25[i]   <- quantile_from_hist(scores, counts_i, 0.25)
      q50[i]   <- quantile_from_hist(scores, counts_i, 0.50)
      q75[i]   <- quantile_from_hist(scores, counts_i, 0.75)
    }
  }
  rc <- s$total_reads
  statdf_list[[length(statdf_list) + 1L]] <- data.frame(
    Cycle = seq_len(n_cycles),
    Mean  = means,
    Q25   = q25,
    Q50   = q50,
    Q75   = q75,
    Cum   = if (rc > 0) 10 * cums / rc else rep(0, n_cycles),
    file  = s$file
  )

  anndf_list[[length(anndf_list) + 1L]] <- data.frame(
    label = s$file,
    rclabel = sprintf("Reads: %d", as.integer(rc)),
    rc = rc,
    file = s$file
  )
}

plotdf <- do.call(rbind, plotdf_list)
statdf <- do.call(rbind, statdf_list)
anndf  <- do.call(rbind, anndf_list)

# ---- Build plot ---------------------------------------------------------

build_plot <- function(plotdf, statdf, anndf) {
  p <- ggplot(plotdf, aes(x = Cycle, y = Score)) +
    geom_tile(aes(fill = Count), height = 1) +
    scale_fill_gradient(low = "#F5F5F5", high = "black") +
    geom_line(data = statdf, aes(y = Mean), color = "#66C2A5") +
    geom_line(data = statdf, aes(y = Q25), color = "#FC8D62",
              linewidth = 0.25, linetype = "dashed") +
    geom_line(data = statdf, aes(y = Q50), color = "#FC8D62",
              linewidth = 0.25) +
    geom_line(data = statdf, aes(y = Q75), color = "#FC8D62",
              linewidth = 0.25, linetype = "dashed") +
    ylab("Quality Score") + xlab("Cycle") +
    theme_bw() + theme(panel.grid = element_blank()) +
    guides(fill = "none")

  if (length(unique(statdf$Cum)) > 1L) {
    p <- p +
      geom_line(data = statdf, aes(y = Cum), color = "red",
                linewidth = 0.25, linetype = "solid") +
      scale_y_continuous(
        limits = c(0, NA),
        sec.axis = sec_axis(~ . * 10, breaks = c(0, 100),
                            labels = c("0%", "100%"))
      ) +
      theme(axis.text.y.right  = element_text(color = "red"),
            axis.title.y.right = element_text(color = "red"))
  } else {
    p <- p + ylim(c(0, NA))
  }
  p
}

if (aggregate) {
  agg_plotdf <- aggregate(Count ~ Cycle + Score, plotdf, sum)
  agg_plotdf$file <- sprintf("%d files (aggregated)", length(summaries))

  cycles <- sort(unique(agg_plotdf$Cycle))
  agg_stat <- do.call(rbind, lapply(cycles, function(c) {
    sub <- agg_plotdf[agg_plotdf$Cycle == c, ]
    total <- sum(sub$Count)
    data.frame(
      Cycle = c,
      Mean  = if (total > 0) sum(sub$Score * sub$Count) / total else 0,
      Q25   = quantile_from_hist(sub$Score, sub$Count, 0.25),
      Q50   = quantile_from_hist(sub$Score, sub$Count, 0.50),
      Q75   = quantile_from_hist(sub$Score, sub$Count, 0.75),
      Cum   = total,
      file  = agg_plotdf$file[1]
    )
  }))
  total_reads <- sum(anndf$rc)
  agg_stat$Cum <- if (total_reads > 0) 10 * agg_stat$Cum / total_reads else 0

  p <- build_plot(agg_plotdf, agg_stat,
                  data.frame(label = agg_plotdf$file[1],
                             rc    = total_reads,
                             file  = agg_plotdf$file[1])) +
    facet_wrap(~ file)
} else {
  p <- build_plot(plotdf, statdf, anndf) +
    geom_text(data = anndf, aes(x = 0, y = 0, label = rclabel),
              color = "red", hjust = 0, inherit.aes = FALSE) +
    facet_wrap(~ file)
}

ggsave(out_file, plot = p, width = width, height = height)
message(sprintf("wrote %s", out_file))
