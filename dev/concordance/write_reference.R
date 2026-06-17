#!/usr/bin/env Rscript
# write_reference.R — generate the static R DADA2 reference for the concordance
# guardrail. Run this ONCE on the small fixture data; commit the resulting CSV.
# CI never runs R — it only diffs dada2-rs output against the committed CSV.
#
# The parameters here MUST match the dada2-rs runner (run_illumina.sh /
# run_pacbio.sh) exactly, or the comparison is apples-to-oranges.
#
# Output schema (the contract compare_to_reference.py reads) — long CSV:
#   sequence,sample,count
#   <ASV nt seq>,<sample name>,<integer count>
# one row per (ASV, sample) with count > 0, from the post-chimera seqtab.nochim.
#
# Usage:
#   # Illumina (paired): <data-dir> holds <sample>F.fastq.gz / <sample>R.fastq.gz
#   Rscript write_reference.R illumina <data-dir> reference/illumina_seqtab_nochim.csv
#
#   # PacBio (single-end, primered raw reads)
#   Rscript write_reference.R pacbio <data-dir> reference/pacbio_seqtab_nochim.csv \
#       AGRGTTYGATYMTGGCTCAG RGYTACCTTGTTACGACTT

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: write_reference.R <illumina|pacbio> <data-dir> <out.csv> [primer_fwd primer_rev]")
platform <- args[1]; data_dir <- args[2]; out_csv <- args[3]

write_long <- function(seqtab, path) {
  # seqtab: matrix rows = samples, cols = sequences (colnames = ASV seqs)
  seqs <- colnames(seqtab)
  rows <- list()
  for (si in seq_len(nrow(seqtab))) {
    sample <- rownames(seqtab)[si]
    for (j in seq_along(seqs)) {
      cnt <- seqtab[si, j]
      if (cnt > 0) rows[[length(rows) + 1]] <- data.frame(
        sequence = seqs[j], sample = sample, count = as.integer(cnt),
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, path, row.names = FALSE, quote = FALSE)
  cat(sprintf("wrote %s: %d ASVs, %d sample(s), %d rows\n",
              path, length(seqs), nrow(seqtab), nrow(df)))
}

if (platform == "illumina") {
  # --- Parameters: keep in sync with run_illumina.sh ---
  TRUNC_LEN <- c(240, 160); MAX_EE <- c(2, 2); TRUNC_Q <- 2; MAX_N <- 0

  fnFs <- sort(list.files(data_dir, pattern = "F\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sub("F\\.fastq\\.gz$", "R\\.fastq\\.gz", fnFs)
  if (length(fnFs) == 0) stop("no *F.fastq.gz in ", data_dir)
  sample.names <- sub("F\\.fastq\\.gz$", "", basename(fnFs))

  filt_dir <- file.path(tempdir(), "filtered")
  filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
  filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = TRUNC_LEN,
                maxN = MAX_N, maxEE = MAX_EE, truncQ = TRUNC_Q,
                rm.phix = FALSE, compress = TRUE, multithread = TRUE)

  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)
  ddF <- dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
  ddR <- dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)
  mergers <- mergePairs(ddF, filtFs, ddR, filtRs)
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = TRUE, verbose = TRUE)
  if (length(sample.names) == 1) rownames(seqtab.nochim) <- sample.names
  write_long(seqtab.nochim, out_csv)

} else if (platform == "pacbio") {
  # --- Parameters: keep in sync with run_pacbio.sh ---
  if (length(args) < 5) stop("pacbio needs primer_fwd primer_rev")
  primer_fwd <- args[4]; primer_rev <- args[5]
  MIN_LEN <- 1000; MAX_LEN <- 1600; MAX_EE <- 2; TRUNC_Q <- 0; MAX_N <- 0
  rc <- getFromNamespace("rc", "dada2")

  fns <- sort(list.files(data_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
  if (length(fns) == 0) stop("no *.fastq.gz in ", data_dir)
  sample.names <- sub("\\.fastq\\.gz$", "", basename(fns))

  nop_dir <- file.path(tempdir(), "noprimers")
  nops <- file.path(nop_dir, paste0(sample.names, "_noprimer.fastq.gz"))
  removePrimers(fns, nops, primer.fwd = primer_fwd, primer.rev = rc(primer_rev),
                orient = TRUE, verbose = TRUE)

  filt_dir <- file.path(tempdir(), "filtered")
  filts <- file.path(filt_dir, paste0(sample.names, "_filt.fastq.gz"))
  filterAndTrim(nops, filts, minLen = MIN_LEN, maxLen = MAX_LEN, maxN = MAX_N,
                maxEE = MAX_EE, truncQ = TRUNC_Q, rm.phix = FALSE,
                compress = TRUE, multithread = TRUE)

  err <- learnErrors(filts, errorEstimationFunction = PacBioErrfun,
                     BAND_SIZE = 32, multithread = TRUE)
  dd <- dada(filts, err = err, pool = FALSE, BAND_SIZE = 32, multithread = TRUE)
  seqtab <- makeSequenceTable(dd)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = TRUE, verbose = TRUE)
  if (length(sample.names) == 1) rownames(seqtab.nochim) <- sample.names
  write_long(seqtab.nochim, out_csv)

} else {
  stop("unknown platform: ", platform, " (expected illumina or pacbio)")
}
