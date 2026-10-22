#!/usr/bin/env Rscript
#
# indel_gam_proto.R — penalised GAM fits of indel rates, scored out of sample (issue #291)
#
# The companion of indel_fit_proto.py: the same tables, the same held-out
# Poisson deviance per event and calibration, for mgcv GAMs whose smoothness is
# chosen by REML. The point is a fit that adapts to how much data a table has:
# flexible where events are plentiful (PacBio), close to flat where they are
# sparse (Illumina), without picking a model per platform by hand.
#
#   insertions   ev ~ s(Q) [+ cap]                     per class (novel, run-lengthening)
#   deletions    additive:     ev ~ L + s(Q)
#                varying:      ev ~ L + s(Q, by = L)
#                hierarchical: ev ~ L + s(Q) + s(Q, L, bs = "fs", m = 1)   shared curve +
#                              penalised run-length deviations
#                surface:      ev ~ te(Q, L)          L numeric: neighbouring run lengths
#                              borrow from each other (also on log L)
#                monotone:     scam, s(Q, L, bs = "tesmi2")  surface constrained to rise with L only
#                pooled:       per-cell pooling within each run length (the native fit)
#   base         ev ~ L + s(base, bs = "re") [+ s(L:base, bs = "re")]   shrinkable base terms
#
# All with offset(log(opportunities)) in the formula (mgcv's predict ignores an
# offset passed as an argument), family poisson. `cap` is an indicator for
# the capped top quality (PacBio Q93), fitted apart from the smooth.
#
# Usage:
#   Rscript indel_gam_proto.R train.tsv test.tsv [train_samples] [test_samples] [cap_q]
#   (samples: comma lists or "-" for all rows)

suppressPackageStartupMessages(library(mgcv))
have_scam <- requireNamespace("scam", quietly = TRUE)
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: indel_gam_proto.R train.tsv test.tsv [train_samples] [test_samples] [cap_q]")
pick <- function(i) if (length(a) >= i && a[i] != "-") a[i] else NULL
cap_q <- if (!is.null(pick(5))) as.integer(a[5]) else NA

load <- function(path, samples) {
  d <- read.delim(path, comment.char = "#", stringsAsFactors = FALSE)
  if (!is.null(samples)) d <- d[d$sample %in% strsplit(samples, ",")[[1]], ]
  aggregate(cbind(opportunities, events) ~ table + key, data = d, FUN = sum)
}
train <- load(a[1], pick(3))
test <- load(a[2], pick(4))

dev_per_event <- function(obs, pred) {
  pred <- pmax(pred, 1e-300)
  t <- ifelse(obs > 0, obs * log(obs / pred), 0) - (obs - pred)
  2 * sum(t) / max(sum(obs), 1)
}
calib <- function(obs, pred, group, labels) {
  out <- c()
  for (g in labels) {
    s <- group == g
    if (sum(obs[s]) >= 5) out <- c(out, sprintf("%s:%.2f", g, sum(pred[s]) / sum(obs[s])))
  }
  paste(out, collapse = " ")
}
qband <- function(q) as.character(cut(q, c(-1, 19, 29, 39, 59, 79, 92, 93),
                                      labels = c("Q0-19", "Q20-29", "Q30-39", "Q40-59", "Q60-79", "Q80-92", "Q93")))
edf <- function(m) sprintf("%.1f", sum(m$edf))

cat("== insertions: GAM s(Q), REML; held-out deviance per event; calibration pred/obs\n")
for (tab in c("ins_q_novel", "ins_q_run")) {
  tr <- train[train$table == tab & train$opportunities > 0, ]
  te <- test[test$table == tab & test$opportunities > 0, ]
  te <- te[te$key >= min(tr$key) & te$key <= max(tr$key), ]
  tr$cap <- as.numeric(!is.na(cap_q) & tr$key == cap_q)
  te$cap <- as.numeric(!is.na(cap_q) & te$key == cap_q)
  f <- if (sum(tr$cap) > 0) events ~ s(key, k = 15) + cap + offset(log(opportunities)) else
    events ~ s(key, k = 15) + offset(log(opportunities))
  m <- gam(f, family = poisson, data = tr, method = "REML")
  pred <- as.numeric(predict(m, newdata = te, type = "response"))
  cat(sprintf("   %-12s (test events %d) gam  edf=%s dev/event=%7.3f   %s\n", tab, sum(te$events), edf(m),
              dev_per_event(te$events, pred), calib(te$events, pred, qband(te$key), unique(qband(te$key)))))
}

joint <- function(d) {
  x <- d[grepl("^del_fq_L", d$table) & d$opportunities > 0, ]
  x$L <- factor(substring(x$table, nchar(x$table)))
  x$Ln <- as.numeric(as.character(x$L))
  x$lL <- log(x$Ln)
  x
}
# Native pooling: widen a window over neighbouring keys until it holds 20 events.
pool_pred <- function(tr, te) {
  out <- numeric(nrow(te))
  for (L in unique(te$L)) {
    x <- tr[tr$L == L, ]
    n <- max(x$key) + 1
    o <- e <- numeric(n); o[x$key + 1] <- x$opportunities; e[x$key + 1] <- x$events
    r <- sapply(seq_len(n), function(k) {
      w <- 0
      repeat {
        lo <- max(1, k - w); hi <- min(n, k + w)
        if (sum(e[lo:hi]) >= 20 || (lo == 1 && hi == n))
          return(if (sum(o[lo:hi]) > 0) sum(e[lo:hi]) / sum(o[lo:hi]) else 0)
        w <- w + 1
      }
    })
    s <- te$L == L
    out[s] <- te$opportunities[s] * r[pmin(te$key[s], n - 1) + 1]
  }
  out
}
cat("== deletions (run length x flanking Q), REML\n")
tr <- joint(train); te <- joint(test)
te <- te[te$key >= min(tr$key) & te$key <= max(tr$key) & te$L %in% levels(tr$L), ]
te$L <- factor(te$L, levels = levels(tr$L))
# hierarchical: a shared smooth plus run-length deviations penalised towards
# zero (factor-smooth, m = 1), so REML can sit anywhere between the two.
forms <- list(additive = events ~ L + s(key, k = 12) + offset(log(opportunities)),
              varying = events ~ L + s(key, by = L, k = 10) + offset(log(opportunities)),
              hierarchical = events ~ L + s(key, k = 12) + s(key, L, bs = "fs", m = 1, k = 8) +
                offset(log(opportunities)))
kL <- min(5, nlevels(tr$L))
forms$surface <- events ~ te(key, Ln, k = c(10, kL)) + offset(log(opportunities))
forms$surface_logL <- events ~ te(key, lL, k = c(10, kL)) + offset(log(opportunities))
report <- function(name, pred, edf_s) {
  cat(sprintf("   %-13s edf=%5s dev/event=%7.3f   %s\n", name, edf_s, dev_per_event(te$events, pred),
              calib(te$events, pred, paste0("L", te$L), paste0("L", levels(te$L)))))
}
for (name in names(forms)) {
  m <- gam(forms[[name]], family = poisson, data = tr, method = "REML")
  report(name, as.numeric(predict(m, newdata = te, type = "response")), edf(m))
}
if (have_scam) {
  m <- tryCatch(scam::scam(events ~ s(key, Ln, bs = "tesmi2", k = c(10, kL)) + offset(log(opportunities)),
                           family = poisson, data = tr), error = function(e) { cat("   monotone: ", conditionMessage(e), "\n"); NULL })
  if (!is.null(m)) report("monotone", as.numeric(predict(m, newdata = te, type = "response")), edf(m))
}
report("pooled", pool_pred(tr, te), "cells")

base_long <- function(d, kind) {
  x <- d[grepl(paste0("^", kind, "_[ACGT]$"), d$table) & d$opportunities > 0, ]
  x$base <- factor(substring(x$table, nchar(x$table)), levels = c("A", "C", "G", "T"))
  x$L <- factor(pmin(x$key, 5), levels = 1:5)
  x <- aggregate(cbind(opportunities, events) ~ base + L, data = x, FUN = sum)
  x$Lbase <- interaction(x$L, x$base)
  x
}
cat("== base identity as REML-shrunk random effects\n")
for (kind in c("del_hp", "ins_hp")) {
  tr <- base_long(train, kind); te <- base_long(test, kind)
  te <- te[te$L %in% tr$L & te$Lbase %in% tr$Lbase, ]
  out <- c()
  off <- " + offset(log(opportunities))"
  for (spec in list(c("L", "events ~ L"), c("L+re(base)", "events ~ L + s(base, bs = 're')"),
                    c("L+re(base)+re(L:base)", "events ~ L + s(base, bs = 're') + s(Lbase, bs = 're')"))) {
    m <- gam(as.formula(paste0(spec[2], off)), family = poisson, data = tr, method = "REML")
    pred <- as.numeric(predict(m, newdata = te, type = "response"))
    out <- c(out, sprintf("%s: %.3f (edf %s)", spec[1], dev_per_event(te$events, pred), edf(m)))
  }
  cat(sprintf("   %s (test events %d): %s\n", kind, sum(te$events), paste(out, collapse = "   ")))
}
