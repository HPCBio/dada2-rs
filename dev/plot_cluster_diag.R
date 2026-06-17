#!/usr/bin/env Rscript
# plot_cluster_diag.R
#
# Visualise per-iteration cluster diagnostics written by dada2-rs
# errors-from-sample / learn-errors when --diag-dir is supplied.
#
# Usage:
#   Rscript plot_cluster_diag.R <diag_dir> [output.pdf]
#
# <diag_dir> is the directory containing iter_001.json, iter_002.json, …
# If output path is omitted, writes <diag_dir>/cluster_diag.pdf

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_cluster_diag.R <diag_dir> [output.pdf]\n")
  quit(status = 1)
}

diag_dir <- args[1]
out_path <- if (length(args) >= 2) args[2] else
  file.path(diag_dir, "cluster_diag.pdf")

# ---------------------------------------------------------------------------
# Load all iter_NNN.json files
# ---------------------------------------------------------------------------
json_files <- sort(list.files(diag_dir, pattern = "^iter_\\d+\\.json$",
                               full.names = TRUE))
if (length(json_files) == 0) {
  stop("No iter_NNN.json files found in: ", diag_dir)
}

cat(sprintf("Loading %d iteration file(s) from %s\n", length(json_files), diag_dir))

rows <- lapply(json_files, function(f) {
  d <- fromJSON(f)
  # d$samples is a data frame with one row per sample
  s <- d$samples
  data.frame(
    iter        = d$iter,
    converged   = d$converged,
    max_delta   = d$max_delta,
    sample      = s$sample,
    n_clusters  = s$n_clusters,
    total_reads = s$total_reads,
    n_initial   = s$n_initial,
    n_abundance = s$n_abundance,
    n_prior     = s$n_prior,
    n_singleton = s$n_singleton,
    nalign      = s$nalign,
    nshroud     = s$nshroud,
    stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
df$iter   <- as.integer(df$iter)
df$sample <- factor(paste0("sample_", df$sample))

n_samples <- length(unique(df$sample))
n_iters   <- max(df$iter)

cat(sprintf("  %d sample(s), %d iteration(s)\n", n_samples, n_iters))

# Mark converged iteration(s) for annotations
converged_iters <- unique(df$iter[df$converged])

# ---------------------------------------------------------------------------
# Panel 1 — Total cluster count per iteration
# ---------------------------------------------------------------------------
p1 <- ggplot(df, aes(x = iter, y = n_clusters, colour = sample, group = sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq_len(n_iters)) +
  labs(
    title   = "Cluster count per iteration",
    x       = "Iteration",
    y       = "Number of clusters",
    colour  = "Sample"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = if (n_samples > 1) "right" else "none")

if (length(converged_iters) > 0) {
  p1 <- p1 + geom_vline(xintercept = min(converged_iters),
                         linetype = "dashed", colour = "gray50") +
              annotate("text", x = min(converged_iters), y = Inf,
                       label = "converged", vjust = 1.5, hjust = -0.1,
                       size = 3, colour = "gray40")
}

# ---------------------------------------------------------------------------
# Panel 2 — Birth-type breakdown (stacked bar, summed across samples)
# ---------------------------------------------------------------------------
birth_cols <- c(
  "Initial"   = "#4E79A7",
  "Abundance" = "#F28E2B",
  "Prior"     = "#59A14F",
  "Singleton" = "#E15759"
)

# Aggregate across samples for the stacked bar
bt <- aggregate(
  cbind(n_initial, n_abundance, n_prior, n_singleton) ~ iter,
  data = df, FUN = sum
)

# Reshape to long form
bt_long <- data.frame(
  iter       = rep(bt$iter, 4),
  birth_type = rep(c("Initial", "Abundance", "Prior", "Singleton"),
                   each = nrow(bt)),
  count      = c(bt$n_initial, bt$n_abundance, bt$n_prior, bt$n_singleton)
)
bt_long$birth_type <- factor(bt_long$birth_type,
                              levels = c("Initial", "Abundance", "Prior", "Singleton"))

p2 <- ggplot(bt_long, aes(x = iter, y = count, fill = birth_type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_x_continuous(breaks = seq_len(n_iters)) +
  scale_fill_manual(values = birth_cols) +
  labs(
    title  = "Cluster birth-type breakdown per iteration (all samples)",
    x      = "Iteration",
    y      = "Number of clusters",
    fill   = "Birth type"
  ) +
  theme_bw(base_size = 11)

if (length(converged_iters) > 0) {
  p2 <- p2 + geom_vline(xintercept = min(converged_iters),
                         linetype = "dashed", colour = "gray50")
}

# ---------------------------------------------------------------------------
# Panel 3 — Max |err_in − err_out| convergence trace
# ---------------------------------------------------------------------------
conv_df <- unique(df[, c("iter", "max_delta")])

p3 <- ggplot(conv_df, aes(x = iter, y = max_delta)) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  geom_point(size = 2, colour = "steelblue") +
  scale_x_continuous(breaks = seq_len(n_iters)) +
  scale_y_log10() +
  labs(
    title = "Convergence: max |err_in \u2212 err_out| per iteration",
    x     = "Iteration",
    y     = "max |err_in \u2212 err_out|  (log scale)"
  ) +
  theme_bw(base_size = 11)

if (length(converged_iters) > 0) {
  p3 <- p3 + geom_vline(xintercept = min(converged_iters),
                         linetype = "dashed", colour = "gray50")
}

# ---------------------------------------------------------------------------
# Panel 4 — Alignment counts (total aligned vs shrouded) per iteration
# ---------------------------------------------------------------------------
align_df <- aggregate(cbind(nalign, nshroud) ~ iter, data = df, FUN = sum)
align_long <- data.frame(
  iter  = rep(align_df$iter, 2),
  type  = rep(c("Aligned", "Shrouded (k-mer screened)"), each = nrow(align_df)),
  count = c(align_df$nalign, align_df$nshroud)
)

p4 <- ggplot(align_long, aes(x = iter, y = count, colour = type, group = type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq_len(n_iters)) +
  scale_colour_manual(values = c("Aligned" = "#4E79A7",
                                  "Shrouded (k-mer screened)" = "#E15759")) +
  labs(
    title  = "Alignment work per iteration (all samples)",
    x      = "Iteration",
    y      = "Count",
    colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Save all panels to PDF
# ---------------------------------------------------------------------------
pdf(out_path, width = 9, height = 10)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 1)
invisible(dev.off())

cat("Plot written to:", out_path, "\n")
