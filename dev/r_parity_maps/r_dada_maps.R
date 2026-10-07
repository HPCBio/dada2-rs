#!/usr/bin/env Rscript
# r_dada_maps.R — per-sample R dada() under a FIXED error model, keeping what
# a sequence table throws away: which unique went to which cluster ($map) and
# the transition counts ($trans) (#277).
#
# The table comparison can say "one read moved between these two ASVs" but not
# which read, and a learned model can say "trans differs by +-1" but not whose
# alignment. With the model held fixed on both sides, any map or trans
# difference comes from screening, alignment or the divisive loop, not the fit.
#
# Usage:
#   Rscript r_dada_maps.R <err.rds> <out-dir> <fastq>... [--band=32] [--threads=N]
#
# Per sample <s> (FASTQ name minus .fastq.gz and a trailing _filt), in <out-dir>:
#   <s>.r.uniques.tsv  sequence, abundance, center (cluster centre sequence, or NA)
#   <s>.r.trans.tsv    16 rows (A2A..T2T) x quality columns, R's $trans as-is
#   <s>.r.dada.rds     the dada-class object
# compare_maps.py reads the first two against dada2-rs's dada JSON.

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NA_character_) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
BAND <- as.integer(flag("band", "32"))
THREADS <- flag("threads")
MT <- if (!is.na(THREADS)) as.integer(THREADS) else 1L
pos <- args[!grepl("^--", args)]
if (length(pos) < 3) stop("usage: r_dada_maps.R <err.rds> <out-dir> <fastq>... [--band=32] [--threads=N]")
err <- readRDS(pos[1])
out_dir <- pos[2]
fqs <- pos[-(1:2)]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("dada2 %s, BAND_SIZE=%d, %d sample(s), threads %d\n",
            packageVersion("dada2"), BAND, length(fqs), MT))

for (fq in fqs) {
  s <- sub("_filt$", "", sub("\\.fastq\\.gz$", "", basename(fq)))
  drp <- derepFastq(fq)
  # Defaults otherwise, as write_reference.R's PacBio branch calls dada().
  dd <- dada(drp, err = err, BAND_SIZE = BAND, multithread = MT, verbose = FALSE)

  seqs <- names(drp$uniques)
  center <- dd$sequence[dd$map]  # NA where a unique was not assigned
  write.table(data.frame(sequence = seqs, abundance = as.integer(drp$uniques), center = center),
              file.path(out_dir, paste0(s, ".r.uniques.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  tr <- dd$trans
  write.table(data.frame(trans = rownames(tr), tr, check.names = FALSE),
              file.path(out_dir, paste0(s, ".r.trans.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(dd, file.path(out_dir, paste0(s, ".r.dada.rds")))
  cat(sprintf("%s: %d uniques, %d reads, %d ASVs\n",
              s, length(seqs), sum(drp$uniques), length(dd$sequence)))
}
