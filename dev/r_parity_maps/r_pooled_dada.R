#!/usr/bin/env Rscript
# r_pooled_dada.R -- pooled dada() on one read direction with a SAVED error
# model, keeping the dada objects for r_pooled_maps.R (#284). Skips learnErrors,
# the other direction, merging and chimera removal.
#
# Usage:
#   Rscript r_pooled_dada.R <err.rds> <ref.dada.rds> <out.dd.rds> <threads> <filt.fastq.gz>...
# Pass the filtered FASTQs in the reference's order (sort(list.files()) order).
# <ref.dada.rds> is write_reference.R's per-direction table (.dadaF.rds or
# .dadaR.rds); the run must reproduce it exactly, or its maps say nothing about
# the reference.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: r_pooled_dada.R <err.rds> <ref.dadaF.rds> <out.rds> <threads> <fastq>...")
err <- readRDS(args[1]); ref <- readRDS(args[2]); out <- args[3]
threads <- as.integer(args[4]); fqs <- args[-(1:4)]
stopifnot(all(file.exists(fqs)))

cat(sprintf("dada2 %s, %s, %d files, first %s\n", packageVersion("dada2"),
            R.version$platform, length(fqs), basename(fqs[1])))
dd <- dada(fqs, err = err, pool = TRUE, multithread = threads)
saveRDS(dd, out)
cat(sprintf("wrote %s\n", out))

# Same table as the reference, cell for cell?
st <- makeSequenceTable(dd)
if (!setequal(rownames(st), rownames(ref))) stop("sample names differ from the reference: check file names")
if (!setequal(colnames(st), colnames(ref))) {
  stop(sprintf("ASV sets differ: %d only here, %d only in the reference",
               length(setdiff(colnames(st), colnames(ref))), length(setdiff(colnames(ref), colnames(st)))))
}
d <- st[rownames(ref), colnames(ref)] - ref
cat(sprintf("vs reference: %d cells differ, L1 %d\n", sum(d != 0), sum(abs(d))))
if (any(d != 0)) stop("does not reproduce the reference; its maps would not explain it")
cat("reproduces the reference exactly\n")
