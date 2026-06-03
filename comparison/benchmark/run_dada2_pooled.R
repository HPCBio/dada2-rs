#!/usr/bin/env Rscript
# run_dada2_pooled.R
# ---------------------------------------------------------------------------
# Reference R DADA2 full pipeline with POOLED denoising (pool=TRUE), used as
# the R side of the dada2-rs pooled benchmark (see bench_pooled.py).
#
# pool=TRUE is the worst case for runtime/memory: all per-sample uniques are
# pooled and denoised together, so this is the most demanding comparison.
#
# Two platforms:
#   illumina : paired-end. filterAndTrim -> learnErrors(F,R) ->
#              dada(pool=TRUE) -> mergePairs -> makeSequenceTable ->
#              removeBimeraDenovo
#   pacbio   : single-end long reads. filterAndTrim(minLen/maxLen) ->
#              learnErrors(PacBioErrfun, BAND_SIZE=32) ->
#              dada(pool=TRUE, BAND_SIZE=32, HOMOPOLYMER_GAP_PENALTY=-1) ->
#              makeSequenceTable -> removeBimeraDenovo
#
# Per-step elapsed wall time is printed in a parser-friendly form:
#   BENCH_STEP<TAB><name><TAB><elapsed_seconds>
# so the Python driver can break the run down by step. Overall process wall
# time and peak RSS are captured by the driver wrapping this whole script.
#
# Usage (key=value args, all but platform/input/outdir optional):
#   Rscript run_dada2_pooled.R platform=illumina input=/path/raw outdir=/path/out \
#       threads=1 nbases=1e8 fwd_pattern='_R1' rev_pattern='_R2' \
#       trunc_len=240,160 max_ee=2,2 trunc_q=2 max_n=0
#
#   Rscript run_dada2_pooled.R platform=pacbio input=/path/raw outdir=/path/out \
#       threads=1 nbases=1e8 min_len=1000 max_len=1600 max_ee=2 trunc_q=0 \
#       band=32 homo_gap=-1
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(dada2))

## ---- argument parsing (key=value) ----------------------------------------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (!grepl("=", a, fixed = TRUE)) stop(sprintf("bad arg (need key=value): %s", a))
  parts <- strsplit(a, "=", fixed = TRUE)[[1]]
  kv[[parts[1]]] <- paste(parts[-1], collapse = "=")
}
getv  <- function(k, default = NULL) if (!is.null(kv[[k]])) kv[[k]] else default
getn  <- function(k, default) if (!is.null(kv[[k]])) as.numeric(kv[[k]]) else default
getnv <- function(k, default) if (!is.null(kv[[k]])) as.numeric(strsplit(kv[[k]], ",")[[1]]) else default

platform <- getv("platform")
input    <- getv("input")
outdir   <- getv("outdir")
if (is.null(platform) || is.null(input) || is.null(outdir))
  stop("required: platform=, input=, outdir=")

threads <- getn("threads", 1)
nbases  <- getn("nbases", 1e8)
multithread <- if (threads > 1) threads else FALSE
pool_flag <- !identical(getv("pool", "true"), "false")   # TRUE pooled, FALSE per-sample

filt_dir <- file.path(outdir, "filtered_R")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)

## ---- timing helper -------------------------------------------------------
# Runs expr, prints BENCH_STEP line, returns expr's value.
timed <- function(name, expr) {
  t <- system.time(val <- force(expr))
  cat(sprintf("BENCH_STEP\t%s\t%.2f\n", name, t[["elapsed"]]))
  flush(stdout())
  val
}

## =========================================================================
if (platform == "illumina") {
  fwd_pat <- getv("fwd_pattern", "_R1")
  rev_pat <- getv("rev_pattern", "_R2")
  trunc_len <- getnv("trunc_len", c(240, 160))
  max_ee    <- getnv("max_ee", c(2, 2))
  trunc_q   <- getn("trunc_q", 2)
  max_n     <- getn("max_n", 0)

  fnFs <- sort(list.files(input, pattern = fwd_pat, full.names = TRUE))
  fnFs <- fnFs[grepl("\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$", fnFs)]
  fnRs <- sub(fwd_pat, rev_pat, fnFs, fixed = TRUE)
  if (length(fnFs) == 0) stop("no forward reads found")
  if (!all(file.exists(fnRs))) stop("missing reverse reads for: ",
                                    paste(fnFs[!file.exists(fnRs)], collapse = ", "))
  sample.names <- sapply(strsplit(basename(fnFs), fwd_pat, fixed = TRUE), `[`, 1)
  cat(sprintf("illumina: %d paired samples\n", length(fnFs)))

  filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

  timed("filter", filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                                truncLen = trunc_len, maxN = max_n, maxEE = max_ee,
                                truncQ = trunc_q, rm.phix = FALSE,
                                compress = TRUE, multithread = multithread))

  errF <- timed("learn_fwd", learnErrors(filtFs, nbases = nbases, multithread = multithread))
  errR <- timed("learn_rev", learnErrors(filtRs, nbases = nbases, multithread = multithread))

  ddF <- timed("dada_fwd", dada(filtFs, err = errF, pool = pool_flag, multithread = multithread))
  ddR <- timed("dada_rev", dada(filtRs, err = errR, pool = pool_flag, multithread = multithread))

  mergers <- timed("merge", mergePairs(ddF, filtFs, ddR, filtRs))
  seqtab  <- timed("make_table", makeSequenceTable(mergers))

} else if (platform == "pacbio") {
  min_len <- getn("min_len", 1000)
  max_len <- getn("max_len", 1600)
  max_ee  <- getn("max_ee", 2)
  trunc_q <- getn("trunc_q", 0)
  max_n   <- getn("max_n", 0)
  band    <- getn("band", 32)
  homo    <- getn("homo_gap", -1)
  pfwd    <- getv("primer_fwd"); prev <- getv("primer_rev")
  mm      <- getn("max_mismatch", 2)
  if (is.null(pfwd) || is.null(prev)) stop("pacbio needs primer_fwd= and primer_rev=")

  fns <- sort(list.files(input, pattern = "\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$",
                         full.names = TRUE))
  if (length(fns) == 0) stop("no reads found")
  sample.names <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fns))
  cat(sprintf("pacbio: %d samples\n", length(fns)))

  nop_dir <- file.path(outdir, "noprimers_R"); dir.create(nop_dir, showWarnings = FALSE)
  nops  <- file.path(nop_dir, paste0(sample.names, "_noprimer.fastq.gz"))
  filts <- file.path(filt_dir, paste0(sample.names, "_filt.fastq.gz"))
  rc <- utils::getFromNamespace("rc", "dada2")

  # primer.rev supplied 5'->3' (catalog); removePrimers wants it RC'd.
  timed("remove_primers", removePrimers(fns, nops, primer.fwd = pfwd,
                                        primer.rev = rc(prev), orient = TRUE,
                                        max.mismatch = mm, verbose = TRUE))

  timed("filter", filterAndTrim(nops, filts, minLen = min_len, maxLen = max_len,
                                maxN = max_n, maxEE = max_ee, truncQ = trunc_q,
                                rm.phix = FALSE, compress = TRUE,
                                multithread = multithread))

  err <- timed("learn", learnErrors(filts, nbases = nbases,
                                    errorEstimationFunction = PacBioErrfun,
                                    BAND_SIZE = band, multithread = multithread))

  dd <- timed("dada", dada(filts, err = err, pool = pool_flag, BAND_SIZE = band,
                           HOMOPOLYMER_GAP_PENALTY = homo, multithread = multithread))

  seqtab <- timed("make_table", makeSequenceTable(dd))

} else {
  stop(sprintf("unknown platform: %s (use illumina or pacbio)", platform))
}

seqtab.nochim <- timed("remove_bimera",
                       removeBimeraDenovo(seqtab, method = "consensus",
                                          multithread = multithread, verbose = TRUE))

## ---- report --------------------------------------------------------------
cat(sprintf("R pipeline done: %d samples x %d ASVs (%d after chimera removal)\n",
            nrow(seqtab.nochim), ncol(seqtab), ncol(seqtab.nochim)))
saveRDS(seqtab.nochim, file.path(outdir, "seqtab_nochim_R.rds"))
cat(sprintf("BENCH_RESULT\tn_asv\t%d\n", ncol(seqtab.nochim)))
