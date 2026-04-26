#!/usr/bin/env Rscript
# loess_reference.R
#
# Reference implementation of R DADA2's stock loessErrfun, intended as a
# back-door comparison test against dada2-rs's native `loess_errfun` Rust port.
# The contract is the dada2-rs --errfun external wire format:
#
#   args[1] = path to input trans TSV  (16 rows: A2A,A2C,…,T2T; nq quality cols)
#   args[2] = path to output err TSV   (same shape)
#
# Usage from the dada2-rs CLI:
#   dada2-rs learn-errors *.fastq.gz \
#     --errfun external \
#     --errfun-cmd "Rscript examples/external_errfun/loess_reference.R" \
#     -o err_R.json
#
# Verbatim port of dada2:::loessErrfun (no DADA2 install required, only
# base R + stats).

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
        mod.lo <- loess(rlogp ~ q, df, weights = tot)
        pred   <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  # Reconstruct the full 16-row matrix with self-transition probabilities.
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
