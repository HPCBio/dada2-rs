#!/usr/bin/env Rscript
# seqtab_to_r.R — run R DADA2's removeBimeraDenovo on a dada2-rs seqtab JSON.
#
# The companion to rtab_to_seqtab.py, which goes the other way. Together they
# make a 2x2 that neither pipeline can produce alone:
#
#              | our chimera | R chimera
#   R input    |     ?       |  R's own
#   our input  |  our own    |     ?
#
# A same-pipeline comparison cannot separate an input difference from an
# implementation difference at this step. Cross-feeding can: if both
# implementations agree on the same input, the gap is upstream of chimera
# removal. On the 95-sample PacBio run that is exactly what happened -- 2046 vs
# 2046 on R's input and 2078 vs 2078 on ours, exact set identity both ways --
# which moved 13 ASVs out of "chimera difference" and into "denoising".
#
# Usage:
#   Rscript seqtab_to_r.R <dada2-rs seqtab.json> <out_asvs.txt> [method] [threads]
#
# Parameters must match whatever the dada2-rs side used, or the comparison is
# measuring the parameters. Defaults mirror run_pacbio.sh / run_illumina.sh.

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: seqtab_to_r.R <seqtab.json> <out_asvs.txt> [method] [threads]")
}
method <- if (length(args) >= 3) args[3] else "consensus"
mt <- if (length(args) >= 4) as.integer(args[4]) else TRUE

j <- fromJSON(args[1])
# `make-sequence-table` output is a plain object; a Tagged wrapper nests it.
if (!is.null(j$data)) j <- j$data

m <- as.matrix(j$counts)
rownames(m) <- j$samples
colnames(m) <- j$sequences
storage.mode(m) <- "integer"
cat(sprintf("input: %d ASVs x %d samples, %s reads\n",
            ncol(m), nrow(m), format(sum(m), big.mark = ",")))

nochim <- removeBimeraDenovo(m, method = method, multithread = mt, verbose = TRUE)
cat(sprintf("after R removeBimeraDenovo(method=%s): %d ASVs, %s reads\n",
            method, ncol(nochim), format(sum(nochim), big.mark = ",")))
writeLines(colnames(nochim), args[2])
cat(sprintf("wrote %s\n", args[2]))
