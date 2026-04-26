#!/usr/bin/env Rscript
# loess_modified.R
#
# A modified loessErrfun that addresses the non-monotonic behavior often seen
# with NovaSeq / quality-binned data. Two changes vs the stock version:
#
#  1. loess uses span=2 and log10(tot) weights (Salazar/Ruscheweyh, #938).
#  2. Each cell's value is clamped to be no smaller than its column's value
#     at Q40 (Holland-Moritz et al., #791) — enforces monotonicity at the
#     high-quality tail.
#
# Reference: https://github.com/benjjneb/dada2/issues/1307#issuecomment-821190155
#
# Usage from the dada2-rs CLI:
#   dada2-rs learn-errors *.fastq.gz \
#     --errfun external \
#     --errfun-cmd "Rscript examples/external_errfun/loess_modified.R" \
#     -o err_modified.json
#
# Requires: dplyr, magrittr.

suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
})

args  <- commandArgs(trailingOnly = TRUE)
trans <- as.matrix(read.table(args[1], sep = "\t", header = TRUE,
                              row.names = 1, check.names = FALSE))

loessErrfun_mod <- function(trans) {
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

        # Salazar/Ruscheweyh: log10-weighted, span=2.
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)

        pred   <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }

  # Clamp to a sensible probability range.
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  # Holland-Moritz monotonicity: floor every cell at its row's Q40 value.
  estorig <- est
  est <- est %>%
    data.frame(check.names = FALSE) %>%
    mutate(across(everything(), ~ ifelse(. < .data[["40"]], .data[["40"]], .))) %>%
    as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)

  err <- rbind(1 - colSums(est[1:3, ]),  est[1:3, ],
               est[4, ],     1 - colSums(est[4:6, ]),  est[5:6, ],
               est[7:8, ],   1 - colSums(est[7:9, ]),  est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2",
                          c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}

err <- loessErrfun_mod(trans)
write.table(err, args[2], sep = "\t", quote = FALSE, col.names = NA)
