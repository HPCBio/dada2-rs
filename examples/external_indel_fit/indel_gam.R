#!/usr/bin/env Rscript
# indel_gam.R
#
# Penalised GAM fit of indel rates for the experimental indel-aware λ
# (issue #291), for `dada2-rs learn-errors --indel-fit-cmd`. The wire format:
#
#   args[1] = input counts TSV  (sample, table, key, opportunities, events;
#                                the --indel-opportunities layout)
#   args[2] = output rates TSV  (table, key, rate)
#
# Usage from the dada2-rs CLI:
#   dada2-rs learn-errors *.fastq.gz -o err.json \
#     --indel-rates-out indel_rates.tsv \
#     --indel-fit-cmd "Rscript examples/external_indel_fit/indel_gam.R"
#
# Models (mgcv, family poisson, smoothness chosen by REML):
#
#   insertions   events ~ s(Q) [+ cap]        per class: ins_q_novel, ins_q_run
#   deletions    events ~ L + s(Q) + s(Q, L, bs = "fs", m = 1)
#                one shared curve in flanking Q plus run-length deviations
#                penalised towards it (L = centre run length 1-4, 5+); used
#                for runs of 1 and 2 only
#
# Runs of 3 or more take dada2-rs's native pooled rates instead (each quality
# widened over its neighbours until it holds 20 events). Their events are
# sparse on Illumina, so REML shrinks their deviations to zero and they inherit
# the bumps of the shared curve, fitted to runs of 1 and 2, at qualities where
# they have no events at all (#291).
#
# `cap` is an indicator for a capped top quality (PacBio HiFi's Q93), fitted
# apart from the smooth; it is used when the top key is 90 or more. These won
# the held-out comparisons on #291 (dev/indel_gam_proto.R): REML lets the fit
# be flexible where events are plentiful (PacBio) and close to flat where they
# are sparse (Illumina). Tables too small to fit fall back to pooled rates.
#
# The offset sits in the formula: mgcv's predict ignores an offset passed as an
# argument.
#
# Dependencies: mgcv (a recommended package, shipped with R).

suppressPackageStartupMessages(library(mgcv))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: indel_gam.R <counts.tsv> <rates.tsv>")

d <- read.delim(args[1], comment.char = "#", stringsAsFactors = FALSE)
d <- aggregate(cbind(opportunities, events) ~ table + key, data = d, FUN = sum)
lens <- d[d$table %in% c("ins_len", "del_len"), ]
d <- d[d$opportunities > 0, ]

out <- data.frame(table = character(), key = integer(), rate = numeric())
emit <- function(table, key, rate) {
  out <<- rbind(out, data.frame(table = table, key = key, rate = pmax(rate, 0)))
}
pooled <- function(x) sum(x$events) / sum(x$opportunities)
GAM_RUNS <- 1:2

# Native pooling (IndelModel::from_counts): the rate at each key from 0 to the
# highest, widening a window over neighbouring keys until it holds `min_events`.
pooled_keys <- function(x, min_events = 20) {
  n <- max(x$key) + 1
  o <- e <- numeric(n)
  o[x$key + 1] <- x$opportunities
  e[x$key + 1] <- x$events
  sapply(seq_len(n), function(k) {
    w <- 0
    repeat {
      lo <- max(1, k - w); hi <- min(n, k + w)
      if (sum(e[lo:hi]) >= min_events || (lo == 1 && hi == n))
        return(if (sum(o[lo:hi]) > 0) sum(e[lo:hi]) / sum(o[lo:hi]) else 0)
      w <- w + 1
    }
  })
}

# Rate at every key from the lowest to the highest seen; keys outside the
# data take the nearest fitted value rather than an extrapolated smooth.
fit_ins <- function(tab) {
  x <- d[d$table == tab, ]
  if (nrow(x) == 0) stop("no `", tab, "` table in ", args[1])
  keys <- 0:max(x$key)
  cap_q <- max(x$key)
  x$cap <- as.numeric(cap_q >= 90 & x$key == cap_q)
  nk <- length(unique(x$key[x$cap == 0]))
  rate <- if (sum(x$events) == 0 || nk < 5) {
    rep(pooled(x), length(keys))
  } else {
    f <- if (sum(x$cap) > 0 && sum(x$events[x$cap == 0]) > 0)
      events ~ s(key, k = min(15, nk - 1)) + cap + offset(log(opportunities)) else
      events ~ s(key, k = min(15, nk - 1)) + offset(log(opportunities))
    m <- gam(f, family = poisson, data = x, method = "REML")
    lo <- min(x$key)
    nd <- data.frame(key = pmax(keys, lo), opportunities = 1,
                     cap = as.numeric(cap_q >= 90 & keys == cap_q))
    as.numeric(predict(m, newdata = nd, type = "response"))
  }
  emit(tab, keys, rate)
}
fit_ins("ins_q_novel")
fit_ins("ins_q_run")

x <- d[grepl("^del_fq_L[1-5]$", d$table), ]
if (!"del_fq_L1" %in% x$table) stop("no `del_fq_L1` table in ", args[1])
x$L <- factor(substring(x$table, nchar(x$table)))
keys <- 0:max(x$key)
nk <- length(unique(x$key))
fit_del <- function(form) {
  # The shared smooth and the factor-smooth share `key` by design; mgcv notes it.
  quiet <- function(w) if (grepl("repeated 1-d smooths", conditionMessage(w))) invokeRestart("muffleWarning")
  tryCatch(withCallingHandlers(gam(form, family = poisson, data = x, method = "REML"), warning = quiet),
           error = function(e) NULL)
}
m <- NULL
if (sum(x$events) > 0 && nk >= 5) {
  if (nlevels(x$L) > 1) {
    m <- fit_del(events ~ L + s(key, k = min(12, nk - 1)) +
                   s(key, L, bs = "fs", m = 1, k = min(8, nk - 1)) + offset(log(opportunities)))
    if (is.null(m)) m <- fit_del(events ~ L + s(key, k = min(12, nk - 1)) + offset(log(opportunities)))
  } else {
    m <- fit_del(events ~ s(key, k = min(12, nk - 1)) + offset(log(opportunities)))
  }
}
prev <- NULL
for (L in 1:5) {
  tab <- paste0("del_fq_L", L)
  xl <- x[x$table == tab, ]
  if (nrow(xl) == 0) {
    # No runs this long: borrow the class below, as the native fit does.
    rate <- prev$rate
    run_rate <- prev$run_rate
  } else {
    run_rate <- pooled(xl)
    rate <- if (!L %in% GAM_RUNS) {
      p <- pooled_keys(xl)
      p[pmin(keys, length(p) - 1) + 1]
    } else if (is.null(m)) rep(run_rate, length(keys)) else {
      nd <- data.frame(key = pmin(pmax(keys, min(xl$key)), max(xl$key)), opportunities = 1,
                       L = factor(as.character(L), levels = levels(x$L)))
      as.numeric(predict(m, newdata = nd, type = "response"))
    }
  }
  emit(tab, keys, rate)
  emit("del_L", L, run_rate)
  prev <- list(rate = rate, run_rate = run_rate)
}

# r of a geometric length model: 1 - 1 / mean length.
length_factor <- function(tab) {
  x <- lens[lens$table == tab, ]
  n <- sum(x$events)
  if (n == 0) return(1e-6)
  max(1 - 1 / (sum(x$key * x$events) / n), 1e-6)
}
emit("ins_r", 0, length_factor("ins_len"))
emit("del_r", 0, length_factor("del_len"))

out$rate <- sprintf("%.6e", out$rate)
write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
