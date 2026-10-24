#!/usr/bin/env Rscript
#
# hp_structure_proto.R — does structure explain homopolymer slippage? (issue #291)
#
# Reads hp_structure_features.py output (one row per parent x run x edit:
# variant reads, error-free parent reads, structure covariates) and fits binomial
# GAMs for the slippage rate, scored on held-out data by deviance per variant
# read (lower is better):
#
#   base        edit x base x run length (3..6, 7+)
#   +rfam       + consensus pairing (fraction paired, stem/loop junction) by edit
#   +plf37/72   + local pairing probability (RNAplfold, DNA parameters) by edit
#   +ir         + inverted-repeat stem length by edit
#   +ctx        + flanking 3-mers either side, REML-shrunk random effect by edit
#   +ctx+X      context together with the best structure term
#
# Splits: train on one dataset group and test on the other, and odd vs even
# parents within the first dataset. Rows with fewer than --min-exact error-free
# parent reads are dropped.
#
# Usage: Rscript hp_structure_proto.R features.tsv [min_exact]

suppressPackageStartupMessages(library(mgcv))
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(a[1], stringsAsFactors = FALSE)
min_exact <- if (length(a) >= 2) as.integer(a[2]) else 300
d <- d[!is.na(d$frac) & d$parent_exact >= min_exact, ]
d$Lc <- factor(pmin(d$run_len, 7))
d$edit <- factor(d$edit)
d$cell <- interaction(d$edit, d$base, d$Lc, drop = TRUE)
d$ctx <- factor(paste(d$edit, substring(d$left, nchar(d$left) - 2), d$base, substr(d$right, 1, 3), sep = "_"))
d$ir <- pmin(d$ir_stem, 10)
d$y <- cbind(d$variant_reads, d$parent_exact)
pb <- unique(d$dataset)[1]

dev_per <- function(obs, n, p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  ll <- function(q) ifelse(obs > 0, obs * log(q), 0) + ifelse(n - obs > 0, (n - obs) * log(1 - q), 0)
  sat <- obs / n
  2 * sum(ll(sat) - ll(p)) / max(sum(obs), 1)
}

forms <- list(
  base = y ~ cell,
  `+rfam` = y ~ cell + edit:rfam_paired + edit:rfam_junction,
  `+plf37` = y ~ cell + s(plf37, by = edit, k = 5),
  `+plf72` = y ~ cell + s(plf72, by = edit, k = 5),
  `+ir` = y ~ cell + s(ir, by = edit, k = 5),
  `+ctx` = y ~ cell + s(ctx, bs = "re"),
  `+ctx+rfam` = y ~ cell + s(ctx, bs = "re") + edit:rfam_paired + edit:rfam_junction,
  `+ctx+plf72` = y ~ cell + s(ctx, bs = "re") + s(plf72, by = edit, k = 5),
  `+ctx+ir` = y ~ cell + s(ctx, bs = "re") + s(ir, by = edit, k = 5))

odd <- as.integer(sub("p", "", d$parent)) %% 2 == 1
splits <- list(list(tr = d$dataset == pb, te = d$dataset != pb),
               list(tr = d$dataset != pb, te = d$dataset == pb),
               list(tr = d$dataset == pb & odd, te = d$dataset == pb & !odd))
names(splits) <- c(sprintf("%s -> others", pb), sprintf("others -> %s", pb), sprintf("%s odd -> even parents", pb))

for (sn in names(splits)) {
  tr <- d[splits[[sn]]$tr, ]; te <- d[splits[[sn]]$te, ]
  te <- te[te$cell %in% tr$cell, ]
  cat(sprintf("== %s: train %d rows (%d variant reads), test %d rows (%d variant reads)\n",
              sn, nrow(tr), sum(tr$variant_reads), nrow(te), sum(te$variant_reads)))
  for (name in names(forms)) {
    m <- tryCatch(gam(forms[[name]], family = binomial, data = tr, method = "REML",
                      drop.unused.levels = FALSE), error = function(e) { cat("  ", name, conditionMessage(e), "\n"); NULL })
    if (is.null(m)) next
    p <- as.numeric(predict(m, newdata = te, type = "response"))
    hot <- te$site %in% c("A", "B") & te$edit == "+1"
    cal_hot <- if (any(hot)) sum(te$parent_exact[hot] * p[hot] / (1 - p[hot])) / max(sum(te$variant_reads[hot]), 1) else NA
    cat(sprintf("   %-11s edf=%6.1f  dev/variant read=%7.4f   hotspot rows pred/obs=%5.2f\n",
                name, sum(m$edf), dev_per(te$variant_reads, te$variant_reads + te$parent_exact, p), cal_hot))
  }
}
