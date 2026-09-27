#!/usr/bin/env Rscript
# learn_errors_r.R
# ---------------------------------------------------------------------------
# Fit ONE error model with R DADA2's learnErrors() on a list of already-filtered
# FASTQ files, and save it as .rds.
#
# This exists for bench_pooled.py's `r-model` error-model arm (#205), which
# compares R's error model against ours on the SAME filtered reads. That is why
# it takes the file list on the command line instead of reading bench_step.R's
# manifest.rds: the whole point is that the reads come from one shared filter
# pass, not from R's own.
#
# It only fits. bench_pooled.py converts the .rds to dada2-rs JSON with the
# shipped scripts/learnerrors_to_dada2rs.R, so there is one conversion
# implementation rather than two.
#
# Usage (key=value args, then the filtered FASTQ paths):
#   Rscript learn_errors_r.R platform=illumina out=/tmp/errF.rds \
#       threads=8 nbases=1e8 sample1_F_filt.fastq.gz sample2_F_filt.fastq.gz
#   Rscript learn_errors_r.R platform=pacbio out=/tmp/err.rds \
#       threads=8 nbases=1e8 band=32 homo_gap=-1 s1_filt.fastq.gz
#
# DADA2 is by Benjamin Callahan et al.; learnErrors() and PacBioErrfun are R
# DADA2's, called here as the reference implementation.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
filts <- character(0)
for (a in args) {
  if (grepl("^[A-Za-z_]+=", a)) {
    parts <- strsplit(a, "=", fixed = TRUE)[[1]]
    kv[[parts[1]]] <- paste(parts[-1], collapse = "=")
  } else {
    filts <- c(filts, a)
  }
}
getv <- function(k, default = NULL) if (!is.null(kv[[k]])) kv[[k]] else default
getn <- function(k, default) if (!is.null(kv[[k]])) as.numeric(kv[[k]]) else default

platform <- getv("platform")
out      <- getv("out")
if (is.null(platform) || is.null(out)) stop("required: platform= out=")
if (length(filts) == 0L) stop("no filtered FASTQ files given")
missing <- filts[!file.exists(filts)]
if (length(missing) > 0L)
  stop("missing filtered file(s): ", paste(missing, collapse = ", "))

threads     <- getn("threads", 1)
multithread <- if (threads > 1) threads else FALSE
nbases      <- getn("nbases", 1e8)

# Files are passed in the order bench_pooled.py hands them over, and nbases
# takes whole samples in that order -- so the arm's training set is the same
# set of samples ours trains on only because the order matches. See the
# --nbases caveat in docs/commands/learn-errors.md.
if (identical(platform, "pacbio")) {
  band <- getn("band", 32)
  homo <- getn("homo_gap", NA)
  err <- learnErrors(filts, nbases = nbases,
                     errorEstimationFunction = PacBioErrfun,
                     BAND_SIZE = band,
                     HOMOPOLYMER_GAP_PENALTY = if (is.na(homo)) NULL else homo,
                     multithread = multithread)
} else {
  err <- learnErrors(filts, nbases = nbases, multithread = multithread)
}

saveRDS(err, out)
cat(sprintf("learn_errors_r: %s, %d file(s), nbases=%g -> %s\n",
            platform, length(filts), nbases, out))
