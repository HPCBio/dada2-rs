#!/usr/bin/env Rscript
# loess_reference_direct.R
#
# Variant of loess_reference.R that forces `loess(..., surface = "direct")`
# instead of the default `surface = "interpolate"`.  Default R loess fits at
# the vertices of an adaptive kd-tree and interpolates between them;
# `surface = "direct"` evaluates the local polynomial directly at every
# requested x — matching dada2-rs's `loess_predict`.
#
# Purpose: isolate whether the residual ~1e-3 absolute diff between
# dada2-rs's native loess and R DADA2 (issue #14) is the kd-tree
# interpolation surface or something else.  If a full `learn-errors` run
# with `--errfun-cmd "Rscript loess_reference_direct.R"` produces an err
# matrix that matches the native-loess run to machine epsilon, the
# kd-tree is confirmed as the sole remaining source.
#
# Wire format matches loess_reference.R:
#   args[1] = path to input trans TSV  (16 rows: A2A,A2C,…,T2T; nq cols)
#   args[2] = path to output err TSV   (same shape)

args  <- commandArgs(trailingOnly = TRUE)
trans <- as.matrix(read.table(args[1], sep = "\t", header = TRUE,
                              row.names = 1, check.names = FALSE))

loessErrfun <- function(trans) {
  qq  <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs  <- trans[paste0(nti, "2", ntj), ]
        tot   <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df    <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        # surface = "direct" evaluates the local polynomial at every query
        # point — including outside the data range, where the polynomial
        # extrapolates wildly.  surface = "interpolate" naturally returns NA
        # there (via the kd-tree).  Replicate that explicitly so the existing
        # flat-fill logic still triggers, matching dada2-rs's `loess_predict`
        # which short-circuits to None outside [x_min, x_max].
        mod.lo <- loess(rlogp ~ q, df, weights = tot, surface = "direct")
        pred   <- predict(mod.lo, qq)
        valid_q <- df$q[!is.na(df$rlogp) & df$tot > 0]
        pred[qq < min(valid_q) | qq > max(valid_q)] <- NA
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  # Post-fit clamp from dada2 errorModels.R:53-56.
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  err <- rbind(1 - colSums(est[1:3, ]),  est[1:3, ],
               est[4, ],     1 - colSums(est[4:6, ]),  est[5:6, ],
               est[7:8, ],   1 - colSums(est[7:9, ]),  est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2",
                          c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}

err <- loessErrfun(trans)
write.table(err, args[2], sep = "\t", quote = FALSE, col.names = NA)
