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
#
#   # PacBio against run_pacbio.sh's own filtered output (the like-for-like form:
#   # ONE filtering pass feeds both tools, and the training set is pinned)
#   Rscript write_reference.R pacbio <out-dir>/filtered pacbio_ref.csv \
#       --prefiltered --train-samples=train.txt --nbases=1e12 --threads=24

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop(paste(
  "usage: write_reference.R <illumina|pacbio> <data-dir> <out.csv>",
  "[primer_fwd primer_rev] [--pool=false|pseudo|true] [--errfun=loess|binned-qual]",
  "[--binned-quals=2,11,25,37] [--prefiltered] [--threads=N] [--nbases=N]",
  "[--train-samples=FILE]",
  "(--pool/--prefiltered/--train-samples/--nbases apply to BOTH platforms)"))
platform <- args[1]; data_dir <- args[2]; out_csv <- args[3]

# --pool=pseudo generates a reference for `dada(pool="pseudo")` instead of the
# per-sample default, so run_illumina.sh POOL=pseudo can be compared against a
# matching R run. Anywhere in the args; default per-sample.
flag <- function(name, default = NA_character_) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}

pool_raw <- flag("pool", "false")
POOL <- switch(pool_raw, pseudo = "pseudo", true = TRUE, TRUE_ = TRUE, FALSE)
cat(sprintf("pool mode: %s\n",
            if (identical(POOL, "pseudo")) "pseudo"
            else if (isTRUE(POOL)) "TRUE (full pooling)" else "FALSE (per-sample)"))

# Error function, mirroring run_illumina.sh's ERRFUN / ERRFUN_ARGS so the two
# sides can be pointed at the same model. `binned-qual` needs its anchors.
ERRFUN <- flag("errfun", "loess")
BINNED <- flag("binned-quals")
ERRFUN_FN <- if (identical(ERRFUN, "binned-qual")) {
  if (is.na(BINNED)) stop("--errfun=binned-qual requires --binned-quals=a,b,c")
  bins <- as.numeric(strsplit(BINNED, ",")[[1]])
  cat(sprintf("errfun: binned-qual, anchors %s\n", paste(bins, collapse = ",")))
  makeBinnedQualErrfun(bins)
} else {
  # NB: the pacbio branch ignores this and uses PacBioErrfun. Say so here
  # rather than print "loess" into a PacBio log, which reads as though the
  # platform default had been overridden.
  if (identical(platform, "pacbio")) {
    cat("errfun: PacBioErrfun (fixed for this platform; --errfun is ignored)\n")
  } else {
    cat("errfun: loess (R default)\n")
  }
  loessErrfun
}

# --prefiltered: the inputs are ALREADY trimmed and filtered, so use them as-is,
# mirroring run_illumina.sh's PREFILTERED. This is what lets ONE filtering pass
# feed both sides, rather than each tool filtering separately and the comparison
# silently carrying that difference too.
PREFILTERED <- any(args == "--prefiltered")
if (PREFILTERED) cat("inputs treated as pre-filtered; skipping filterAndTrim\n")

# Thread count. `multithread = MT` makes DADA2 call parallel::detectCores(),
# which reports the PHYSICAL machine rather than a cgroup or cpuset -- so under
# SLURM it happily spawns one thread per host core against a much smaller
# allocation, oversubscribing the node and making any timing meaningless.
# Prefer an explicit count: --threads=N, else $SLURM_CPUS_PER_TASK, else the old
# TRUE so off-cluster behaviour is unchanged.
THREADS <- flag("threads")
MT <- if (!is.na(THREADS)) {
  as.integer(THREADS)
} else if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else {
  TRUE
}
cat(sprintf("threads: %s\n", if (isTRUE(MT)) "TRUE (detectCores; NOT cgroup-aware)" else MT))

# learnErrors' subsampling budget. R defaults to 1e8; run_illumina.sh defaults to
# 2e7 and run_pacbio.sh to 2e8 -- so left alone the two sides train on DIFFERENT
# amounts of data, and any arm comparison is measuring the budget as well as
# whatever it meant to test.
# Set both, and set them above the run total if the surface or errfun is the
# thing under test.
# --train-samples=FILE: learn the error model from ONLY these samples (one name
# per line), then denoise everything. Pinning the training set is what makes an
# arm comparison mean anything once a run is large enough to subsample: with a
# budget alone, each tool draws its own reads and the arms differ in their
# training data as well as in whatever was under test. It is also closer to real
# use -- nobody trains on a whole run.
TRAIN_SAMPLES <- flag("train-samples")

NBASES <- flag("nbases")
NBASES <- if (!is.na(NBASES)) as.numeric(NBASES) else 1e8
cat(sprintf("nbases: %s%s\n", format(NBASES, scientific = TRUE),
            if (is.na(flag("nbases")))
              " (R default; run_illumina.sh defaults to 2e7, run_pacbio.sh to 2e8 -- match them)"
            else ""))

args <- args[!grepl("^--", args)]

write_long <- function(seqtab, path) {
  # seqtab: matrix rows = samples, cols = sequences (colnames = ASV seqs)
  seqs <- colnames(seqtab)
  # Vectorised: the old form allocated one data.frame per non-zero cell and
  # rbound them, which is fine for a fixture and hopeless for a real table --
  # a 2,823 x 362 pre-chimera matrix has ~10^5-10^6 non-zero cells, and the
  # pre-chimera CSV silently failed to appear on the first 362-sample run
  # because of it.
  nz <- which(seqtab > 0, arr.ind = TRUE)
  df <- data.frame(
    sequence = seqs[nz[, "col"]],
    sample   = rownames(seqtab)[nz[, "row"]],
    count    = as.integer(seqtab[nz]),
    stringsAsFactors = FALSE)
  df <- df[order(df$sample, df$sequence), , drop = FALSE]
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, path, row.names = FALSE, quote = FALSE)
  cat(sprintf("wrote %s: %d ASVs, %d sample(s), %d rows\n",
              path, length(seqs), nrow(seqtab), nrow(df)))
}

# Save what the post-chimera CSV cannot answer later. Two things have already
# been wanted and were not there:
#   * the error models, which `scripts/learnerrors_to_dada2rs.R` converts so our
#     dada can run on R's model -- the arm that separates the error model from
#     everything downstream of it;
#   * the PRE-chimera table, because "R does not have this ASV" and "R had it and
#     called it chimeric" are different findings and the post-chimera CSV cannot
#     tell them apart (three extras on the 362-sample MiSeq run turned out to be
#     clean two-parent bimeras, and confirming that needed R's pre-chimera set).
# --- Training-set pinning, shared by both platforms ----------------------------
# Kept as functions rather than inlined per branch: the Illumina path grew these
# first and the PacBio path silently lacked them, so a PacBio arm comparison was
# not pinned and the nbases guard below could never fire on it.

# Indices into sample.names for the pinned training set, or all of them.
training_index <- function(sample.names) {
  if (is.na(TRAIN_SAMPLES)) return(seq_along(sample.names))
  keep <- trimws(readLines(TRAIN_SAMPLES))
  keep <- keep[nzchar(keep)]
  idx <- match(keep, sample.names)
  if (anyNA(idx)) stop("--train-samples names not found in ", data_dir, ": ",
                       paste(keep[is.na(idx)], collapse = ", "))
  cat(sprintf("training on %d of %d samples (pinned via %s); denoising all %d\n",
              length(idx), length(sample.names), TRAIN_SAMPLES, length(sample.names)))
  idx
}

# nbases must NOT truncate a pinned training set.  The whole point of
# --train-samples is that both tools see exactly these samples; if nbases bites
# first, R silently trains on a PREFIX of the manifest and the arms differ in
# training data as well as in whatever is under test.  This cost a full
# 362-sample comparison once: the manifest held 3.0e8 bases, nbases sat at its
# 1e8 default, and R trained on the first third while dada2-rs (run with
# --nbases 1e12) used all of it -- a 2.99x gap in the transition matrix that
# looked like a self-consistency-loop difference.
assert_nbases_covers <- function(files, label) {
  if (is.na(TRAIN_SAMPLES)) return(invisible(NULL))
  # countFastq reads the headers only -- do NOT pull the records into memory
  # just to total them.
  train_bases <- sum(as.numeric(ShortRead::countFastq(files)$nucleotides))
  cat(sprintf("training set (%s): %d samples, %s bases; nbases = %s\n",
              label, length(files), format(train_bases, big.mark = ","),
              format(NBASES, scientific = TRUE)))
  if (NBASES < train_bases) {
    stop(sprintf(paste0("--nbases (%s) is smaller than the pinned %s training set (%s bases), ",
                        "so learnErrors would train on only a prefix of it. Pass ",
                        "--nbases=%s or larger, or shrink the manifest."),
                 format(NBASES, scientific = TRUE), label,
                 format(train_bases, big.mark = ","),
                 format(ceiling(train_bases * 1.1), scientific = TRUE)))
  }
  invisible(NULL)
}

save_artifacts <- function(out_csv, seqtab_pre, errs) {
  stem <- sub("\\.csv$", "", out_csv)
  for (nm in names(errs)) {
    path <- paste0(stem, ".", nm, ".rds")
    saveRDS(errs[[nm]], path)
    cat(sprintf("wrote %s (learnErrors object; feed to scripts/learnerrors_to_dada2rs.R)\n", path))
  }
  saveRDS(seqtab_pre, paste0(stem, ".prechimera.rds"))
  write_long(seqtab_pre, paste0(stem, ".prechimera.csv"))
}

if (platform == "illumina") {
  # --- Parameters: keep in sync with run_illumina.sh ---
  TRUNC_LEN <- c(240, 160); MAX_EE <- c(2, 2); TRUNC_Q <- 2; MAX_N <- 0

  fnFs <- sort(list.files(data_dir, pattern = "F\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sub("F\\.fastq\\.gz$", "R\\.fastq\\.gz", fnFs)
  if (length(fnFs) == 0) stop("no *F.fastq.gz in ", data_dir)
  sample.names <- sub("F\\.fastq\\.gz$", "", basename(fnFs))

  if (PREFILTERED) {
    filtFs <- fnFs; filtRs <- fnRs
  } else {
    filt_dir <- file.path(tempdir(), "filtered")
    filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
    filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = TRUNC_LEN,
                  maxN = MAX_N, maxEE = MAX_EE, truncQ = TRUNC_Q,
                  rm.phix = FALSE, compress = TRUE, multithread = MT)
  }

  idx <- training_index(sample.names)
  trainFs <- filtFs[idx]; trainRs <- filtRs[idx]
  # Both directions: learnErrors is called separately on each, each with its own
  # nbases, and the reverse set is the smaller of the two -- so checking only
  # the forward would let nbases truncate the reverse unnoticed.
  assert_nbases_covers(trainFs, "forward")
  assert_nbases_covers(trainRs, "reverse")

  errF <- learnErrors(trainFs, errorEstimationFunction = ERRFUN_FN, nbases = NBASES, multithread = MT)
  errR <- learnErrors(trainRs, errorEstimationFunction = ERRFUN_FN, nbases = NBASES, multithread = MT)
  ddF <- dada(filtFs, err = errF, pool = POOL, multithread = MT)
  ddR <- dada(filtRs, err = errR, pool = POOL, multithread = MT)
  mergers <- mergePairs(ddF, filtFs, ddR, filtRs)
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = MT, verbose = TRUE)
  if (length(sample.names) == 1) rownames(seqtab.nochim) <- sample.names
  write_long(seqtab.nochim, out_csv)
  save_artifacts(out_csv, seqtab, list(errF = errF, errR = errR))
  # Per-DIRECTION tables, before mergePairs. Without these a merged-table
  # difference cannot be attributed: "R never called this ASV" may mean R's
  # forward or reverse dada never called the component, or that both were called
  # and the pair failed to merge. Chasing one extra ASV on the 362-sample run
  # ran out of evidence at exactly this point.
  stem <- sub("\\.csv$", "", out_csv)
  for (nm in c("F", "R")) {
    st_dir <- makeSequenceTable(if (nm == "F") ddF else ddR)
    saveRDS(st_dir, paste0(stem, ".dada", nm, ".rds"))
    write_long(st_dir, paste0(stem, ".dada", nm, ".csv"))
  }

} else if (platform == "pacbio") {
  # --- Parameters: keep in sync with run_pacbio.sh ---
  # --prefiltered means the reads are ALREADY primer-stripped and length/EE
  # filtered, so both removePrimers and filterAndTrim are skipped and the
  # primers are not needed. Without it this branch filtered into tempdir(),
  # which is destroyed on exit -- so dada2-rs could never consume the same
  # reads and every PacBio comparison silently carried two filtering passes.
  MIN_LEN <- 1000; MAX_LEN <- 1600; MAX_EE <- 2; TRUNC_Q <- 0; MAX_N <- 0

  fns <- sort(list.files(data_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
  if (length(fns) == 0) stop("no *.fastq.gz in ", data_dir)
  sample.names <- sub("\\.fastq\\.gz$", "", basename(fns))
  # run_pacbio.sh writes its filtered reads as <sample>_filt.fastq.gz into a
  # PERSISTENT $OUT/filtered, which is what --prefiltered is meant to consume.
  # Strip the suffix so these agree with pick_training_subset.py's manifest
  # (run it with --suffix _filt.fastq.gz). Only under --prefiltered, so raw
  # reads that happen to end in _filt are left alone.
  if (PREFILTERED) sample.names <- sub("_filt$", "", sample.names)

  if (PREFILTERED) {
    filts <- fns
  } else {
    if (length(args) < 5) stop("pacbio needs primer_fwd primer_rev (or --prefiltered)")
    primer_fwd <- args[4]; primer_rev <- args[5]
    rc <- getFromNamespace("rc", "dada2")

    nop_dir <- file.path(tempdir(), "noprimers")
    nops <- file.path(nop_dir, paste0(sample.names, "_noprimer.fastq.gz"))
    removePrimers(fns, nops, primer.fwd = primer_fwd, primer.rev = rc(primer_rev),
                  orient = TRUE, verbose = TRUE)

    filt_dir <- file.path(tempdir(), "filtered")
    filts <- file.path(filt_dir, paste0(sample.names, "_filt.fastq.gz"))
    filterAndTrim(nops, filts, minLen = MIN_LEN, maxLen = MAX_LEN, maxN = MAX_N,
                  maxEE = MAX_EE, truncQ = TRUNC_Q, rm.phix = FALSE,
                  compress = TRUE, multithread = MT)
  }

  idx <- training_index(sample.names)
  trains <- filts[idx]
  assert_nbases_covers(trains, "reads")

  err <- learnErrors(trains, errorEstimationFunction = PacBioErrfun, nbases = NBASES,
                     BAND_SIZE = 32, multithread = MT)
  dd <- dada(filts, err = err, pool = POOL, BAND_SIZE = 32, multithread = MT)
  seqtab <- makeSequenceTable(dd)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = MT, verbose = TRUE)
  if (length(sample.names) == 1) rownames(seqtab.nochim) <- sample.names
  write_long(seqtab.nochim, out_csv)
  save_artifacts(out_csv, seqtab, list(err = err))

} else {
  stop("unknown platform: ", platform, " (expected illumina or pacbio)")
}
