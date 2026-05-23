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

# from DADA2's errorModels.R
PacBioErrfun <- function(trans) {
  if("93" %in% colnames(trans)) {
    i.93 <- which(colnames(trans) %in% "93")
    if(i.93 != ncol(trans)) stop("Max qual score of 93 not the last column as expected.")
    err <- loessErrfun(trans[,1:(i.93-1)])
    tot93 <- rep(c(sum(trans[1:4,"93"]), sum(trans[5:8,"93"]), sum(trans[9:12,"93"]), sum(trans[13:16,"93"])), each=4)
    err93 <- (trans[,"93"] + 1)/(tot93 + 4)
    err <- cbind(err, "93"=err93)
  } else {
    message("The max qual score of 93 was not detected. Using standard error fitting.")
    err <- loessErrfun(trans)
  }
  return(err)
}

# from DADA2's errorModels.R
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
  # Post-fit clamp (R DADA2 errorModels.R:53-56 — the `# HACKY` step).
  # Off-diagonal rates are pinned to [1e-7, 0.25] before the diagonal is
  # computed.  Omitting this was a long-running bug in this reference
  # script; see issue #14.
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

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

err <- PacBioErrfun(trans)
write.table(err, args[2], sep = "\t", quote = FALSE, col.names = NA)
