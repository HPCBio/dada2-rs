#!/usr/bin/env Rscript
# bench_step.R
# ---------------------------------------------------------------------------
# Runs ONE step of the R DADA2 pooled pipeline as a standalone process, so the
# Python driver (bench_pooled.py) can capture per-step wall time AND per-step
# peak RSS via os.wait4()/ru_maxrss — symmetric with the dada2-rs side, where
# every step is already its own process.
#
# State is passed between steps as .rds files in <statedir>:
#   manifest.rds  — sample names + file paths (written by `filter`)
#   errF/errR.rds — error models;  ddF/ddR.rds — dada results
#   mergers.rds, seqtab.rds, seqtab_nochim.rds
#
# NOTE on RSS: each step is a fresh R process, so its peak RSS includes the R
# interpreter + dada2 library baseline (~150-200 MB) on top of the step's own
# working set. That is the honest cost of "running this step in R"; small steps
# (filter, make_table) will therefore floor at that baseline.
#
# Steps:
#   illumina : filter learn_fwd learn_rev dada_fwd dada_rev merge make_table remove_bimera
#   pacbio   : filter learn dada make_table remove_bimera
#
# Usage (key=value args):
#   Rscript bench_step.R step=filter platform=illumina input=/raw statedir=/out/Rstate \
#       threads=1 nbases=1e8 fwd_pattern=_R1 rev_pattern=_R2 \
#       trunc_len=240,160 max_ee=2,2 trunc_q=2 max_n=0
#   Rscript bench_step.R step=dada platform=pacbio statedir=/out/Rstate \
#       threads=1 band=32 homo_gap=-1
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(dada2))

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

step     <- getv("step")
platform <- getv("platform")
statedir <- getv("statedir")
if (is.null(step) || is.null(platform) || is.null(statedir))
  stop("required: step=, platform=, statedir=")
dir.create(statedir, showWarnings = FALSE, recursive = TRUE)

threads     <- getn("threads", 1)
nbases      <- getn("nbases", 1e8)
multithread <- if (threads > 1) threads else FALSE
# pool: "true" -> TRUE (pooled), "false" -> FALSE (per-sample), "pseudo" -> "pseudo"
pool_in     <- getv("pool", "true")
pool_flag   <- if (pool_in == "pseudo") "pseudo" else !identical(pool_in, "false")
pseudo_prev <- getn("pseudo_prevalence", 2)
pseudo_abund <- getn("pseudo_min_abundance", NA)
pseudo_abund <- if (is.na(pseudo_abund)) Inf else pseudo_abund   # R PSEUDO_ABUNDANCE default Inf

sp <- function(name) file.path(statedir, name)   # state path helper

timed <- function(name, expr) {
  t <- system.time(val <- force(expr))
  # Leading \n: dada(pool="pseudo") prints progress without a trailing newline,
  # which would otherwise prepend onto this line and defeat the ^-anchored parser.
  cat(sprintf("\nBENCH_STEP\t%s\t%.2f\n", name, t[["elapsed"]]))
  flush(stdout())
  val
}

## =========================================================================
## ILLUMINA
## =========================================================================
if (platform == "illumina") {
  if (step == "filter") {
    input   <- getv("input"); if (is.null(input)) stop("filter needs input=")
    fwd_pat <- getv("fwd_pattern", "_R1"); rev_pat <- getv("rev_pattern", "_R2")
    trunc_len <- getnv("trunc_len", c(240, 160)); max_ee <- getnv("max_ee", c(2, 2))
    trunc_q   <- getn("trunc_q", 2); max_n <- getn("max_n", 0)
    filt_dir  <- sp("filtered_R"); dir.create(filt_dir, showWarnings = FALSE)

    fnFs <- sort(list.files(input, pattern = fwd_pat, full.names = TRUE))
    fnFs <- fnFs[grepl("\\.f(ast)?q(\\.gz)?$", fnFs)]
    fnRs <- sub(fwd_pat, rev_pat, fnFs, fixed = TRUE)
    if (length(fnFs) == 0) stop("no forward reads found")
    if (!all(file.exists(fnRs))) stop("missing reverse mate(s)")
    sample.names <- sapply(strsplit(basename(fnFs), fwd_pat, fixed = TRUE), `[`, 1)
    filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

    timed("filter", filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                                  truncLen = trunc_len, maxN = max_n, maxEE = max_ee,
                                  truncQ = trunc_q, rm.phix = FALSE,
                                  compress = TRUE, multithread = multithread))
    saveRDS(list(sample.names = sample.names, filtFs = filtFs, filtRs = filtRs),
            sp("manifest.rds"))

  } else if (step == "learn_fwd") {
    m <- readRDS(sp("manifest.rds"))
    errF <- timed("learn_fwd", learnErrors(m$filtFs, nbases = nbases, multithread = multithread))
    saveRDS(errF, sp("errF.rds"))

  } else if (step == "learn_rev") {
    m <- readRDS(sp("manifest.rds"))
    errR <- timed("learn_rev", learnErrors(m$filtRs, nbases = nbases, multithread = multithread))
    saveRDS(errR, sp("errR.rds"))

  } else if (step == "dada_fwd") {
    m <- readRDS(sp("manifest.rds")); errF <- readRDS(sp("errF.rds"))
    ddF <- timed("dada_fwd", dada(m$filtFs, err = errF, pool = pool_flag,
                                  PSEUDO_PREVALENCE = pseudo_prev, PSEUDO_ABUNDANCE = pseudo_abund,
                                  multithread = multithread))
    saveRDS(ddF, sp("ddF.rds"))

  } else if (step == "dada_rev") {
    m <- readRDS(sp("manifest.rds")); errR <- readRDS(sp("errR.rds"))
    ddR <- timed("dada_rev", dada(m$filtRs, err = errR, pool = pool_flag,
                                  PSEUDO_PREVALENCE = pseudo_prev, PSEUDO_ABUNDANCE = pseudo_abund,
                                  multithread = multithread))
    saveRDS(ddR, sp("ddR.rds"))

  } else if (step == "merge") {
    m <- readRDS(sp("manifest.rds")); ddF <- readRDS(sp("ddF.rds")); ddR <- readRDS(sp("ddR.rds"))
    mergers <- timed("merge", mergePairs(ddF, m$filtFs, ddR, m$filtRs))
    saveRDS(mergers, sp("mergers.rds"))

  } else if (step == "make_table") {
    mergers <- readRDS(sp("mergers.rds"))
    seqtab <- timed("make_table", makeSequenceTable(mergers))
    saveRDS(seqtab, sp("seqtab.rds"))

  } else if (step == "remove_bimera") {
    seqtab <- readRDS(sp("seqtab.rds"))
    nochim <- timed("remove_bimera", removeBimeraDenovo(seqtab, method = "consensus",
                                                        multithread = multithread, verbose = TRUE))
    saveRDS(nochim, sp("seqtab_nochim.rds"))
    cat(sprintf("\nBENCH_RESULT\tn_asv\t%d\n", ncol(nochim)))

  } else stop(sprintf("unknown illumina step: %s", step))

## =========================================================================
## PACBIO
## =========================================================================
} else if (platform == "pacbio") {
  band <- getn("band", 32); homo <- getn("homo_gap", -1)
  if (step == "remove_primers") {
    input   <- getv("input"); if (is.null(input)) stop("remove_primers needs input=")
    pfwd <- getv("primer_fwd"); prev <- getv("primer_rev")
    if (is.null(pfwd) || is.null(prev)) stop("remove_primers needs primer_fwd= and primer_rev=")
    mm <- getn("max_mismatch", 2)
    nop_dir <- sp("noprimers_R"); dir.create(nop_dir, showWarnings = FALSE)
    rc <- utils::getFromNamespace("rc", "dada2")   # dada2's reverse-complement helper

    fns <- sort(list.files(input, pattern = "\\.f(ast)?q(\\.gz)?$", full.names = TRUE))
    if (length(fns) == 0) stop("no reads found")
    sample.names <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fns))
    nops <- file.path(nop_dir, paste0(sample.names, "_noprimer.fastq.gz"))
    # primer.rev is supplied 5'->3' (catalog); removePrimers expects it in the
    # orientation it appears in reads, i.e. reverse-complemented.
    timed("remove_primers", removePrimers(fns, nops, primer.fwd = pfwd,
                                          primer.rev = rc(prev), orient = TRUE,
                                          max.mismatch = mm, verbose = TRUE))
    saveRDS(list(sample.names = sample.names, nops = nops), sp("manifest.rds"))

  } else if (step == "filter") {
    m <- readRDS(sp("manifest.rds"))
    min_len <- getn("min_len", 1000); max_len <- getn("max_len", 1600)
    max_ee  <- getn("max_ee", 2); trunc_q <- getn("trunc_q", 0); max_n <- getn("max_n", 0)
    filt_dir <- sp("filtered_R"); dir.create(filt_dir, showWarnings = FALSE)

    filts <- file.path(filt_dir, paste0(m$sample.names, "_filt.fastq.gz"))
    timed("filter", filterAndTrim(m$nops, filts, minLen = min_len, maxLen = max_len,
                                  maxN = max_n, maxEE = max_ee, truncQ = trunc_q,
                                  rm.phix = FALSE, compress = TRUE, multithread = multithread))
    m$filts <- filts
    saveRDS(m, sp("manifest.rds"))

  } else if (step == "learn") {
    m <- readRDS(sp("manifest.rds"))
    err <- timed("learn", learnErrors(m$filts, nbases = nbases,
                                      errorEstimationFunction = PacBioErrfun,
                                      BAND_SIZE = band, multithread = multithread))
    saveRDS(err, sp("err.rds"))

  } else if (step == "dada") {
    m <- readRDS(sp("manifest.rds")); err <- readRDS(sp("err.rds"))
    dd <- timed("dada", dada(m$filts, err = err, pool = pool_flag, BAND_SIZE = band,
                             HOMOPOLYMER_GAP_PENALTY = homo,
                             PSEUDO_PREVALENCE = pseudo_prev, PSEUDO_ABUNDANCE = pseudo_abund,
                             multithread = multithread))
    saveRDS(dd, sp("dd.rds"))

  } else if (step == "make_table") {
    dd <- readRDS(sp("dd.rds"))
    seqtab <- timed("make_table", makeSequenceTable(dd))
    saveRDS(seqtab, sp("seqtab.rds"))

  } else if (step == "remove_bimera") {
    seqtab <- readRDS(sp("seqtab.rds"))
    nochim <- timed("remove_bimera", removeBimeraDenovo(seqtab, method = "consensus",
                                                        multithread = multithread, verbose = TRUE))
    saveRDS(nochim, sp("seqtab_nochim.rds"))
    cat(sprintf("\nBENCH_RESULT\tn_asv\t%d\n", ncol(nochim)))

  } else stop(sprintf("unknown pacbio step: %s", step))

} else stop(sprintf("unknown platform: %s", platform))
