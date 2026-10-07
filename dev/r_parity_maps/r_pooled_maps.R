#!/usr/bin/env Rscript
# r_pooled_maps.R — per-sample unique -> cluster tables from a SAVED R dada run
# (write_reference.R --save-dada), in r_dada_maps.R's format, so compare_maps.py
# can read a pooled run too (#277).
#
# Each object's $map indexes that sample's derepFastq() uniques, which the run
# does not keep, so each FASTQ is dereplicated again here. derepFastq is
# deterministic; the abundance check below catches a wrong file.
#
# Usage:
#   Rscript r_pooled_maps.R <dd.rds> <fastq-dir> <out-dir>
# Writes <out-dir>/<s>.r.uniques.tsv (sequence, abundance, center); no trans,
# since a pooled run's $trans is not per sample.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: r_pooled_maps.R <dd.rds> <fastq-dir> <out-dir>")
dd <- readRDS(args[1]); fq_dir <- args[2]; out_dir <- args[3]
if (inherits(dd, "dada")) stop("one dada object, not a list: use r_dada_maps.R for one sample")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (nm in names(dd)) {
  s <- sub("_filt$", "", sub("\\.fastq\\.gz$", "", nm))
  drp <- derepFastq(file.path(fq_dir, nm))
  d <- dd[[nm]]
  if (length(d$map) != length(drp$uniques)) {
    stop(sprintf("%s: $map has %d entries, derep has %d uniques", nm, length(d$map), length(drp$uniques)))
  }
  # Each ASV's reads rebuilt from the map must equal $denoised; a re-derep in a
  # different order would scramble the map and fail here.
  rebuilt <- tapply(as.integer(drp$uniques), factor(d$map, levels = seq_along(d$denoised)), sum)
  rebuilt[is.na(rebuilt)] <- 0L
  if (!all(rebuilt == d$denoised)) stop(sprintf("%s: map + derep do not reproduce $denoised", nm))
  write.table(data.frame(sequence = names(drp$uniques), abundance = as.integer(drp$uniques),
                         center = d$sequence[d$map]),
              file.path(out_dir, paste0(s, ".r.uniques.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s: %d uniques, %d ASVs\n", s, length(drp$uniques), length(d$denoised)))
}
