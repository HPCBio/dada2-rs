#!/usr/bin/env Rscript
# Ground-truth arm for the LOESS oracle harness: R stats::loess(rlogp ~ q,
# weights = tot) on the exact triples dada2-rs feeds its own loess_predict,
# for both of R's surfaces.
#
# Emits the same (pred, rate) pair as the Rust arms, with the same
# flat-extrapolation and [1e-7, 0.25] clamp, so compare.py can join them
# directly.
#
# NOTE: at three or four populated quality columns R fits an exact quadratic and
# oscillates -- it will emit warnings and can return error rates that RISE with
# quality. That is a real property of the reference, not a bug in this script;
# see dev/loess-oracle/README.md.
#
# Usage: Rscript r_truth.R <triples.tsv> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript r_truth.R <triples.tsv> <out.tsv>")
}
inf <- args[1]
ouf <- args[2]

d <- read.delim(inf, na.strings = c("nan", "NA", "NaN"))

flat_fill <- function(p) {
  ok <- which(is.finite(p))
  if (!length(ok)) return(rep(NA_real_, length(p)))
  lo <- min(ok); hi <- max(ok)
  if (lo > 1) p[1:(lo - 1)] <- p[lo]
  if (hi < length(p)) p[(hi + 1):length(p)] <- p[hi]
  p
}

res <- list()
for (pr in unique(d$pair)) {
  df <- d[d$pair == pr, ]
  qq <- df$q
  # R's loessErrfun NA-s non-finite rlogp; tot==0 rows carry weight 0 anyway.
  df$rlogp[!is.finite(df$rlogp)] <- NA
  for (surf in c("direct", "interpolate")) {
    pred <- rep(NA_real_, nrow(df))
    fit <- try(loess(rlogp ~ q, df, weights = tot, surface = surf), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      p <- try(predict(fit, qq), silent = TRUE)
      if (!inherits(p, "try-error")) pred <- as.numeric(p)
    }
    # Outside the fitted data range: direct extrapolates the polynomial;
    # dada2-rs returns None there and flat-fills. Mirror that.
    rng <- range(df$q[is.finite(df$rlogp) & df$tot > 0])
    pred[qq < rng[1] | qq > rng[2]] <- NA
    pred <- flat_fill(pred)
    res[[length(res) + 1]] <- data.frame(
      pair = pr, q = qq, surface = surf, pred = pred,
      rate = pmin(pmax(10^pred, 1e-7), 0.25)
    )
  }
}
out <- do.call(rbind, res)
write.table(out, ouf, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", ouf, nrow(out), "rows\n")
