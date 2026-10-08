#!/usr/bin/env bash
# trace_read.sh — follow one read through R DADA2 and dada2-rs, event by event,
# and show where the two first disagree (#277).
#
# Usage:
#   R_TRACE_LIB=<lib> trace_read.sh <binary> <err.rds> <fastq> <sequence> <out-dir>
# Env:
#   R_TRACE_LIB  R library holding the traced DADA2 (install_r_trace.sh); required
#   ITER         learnErrors pass whose input model to use (default: the last)
#   OMEGA_C (0, as learnErrors), BAND (32), KMER (5), THREADS (1)
#
# Both sides run per-sample dada on <fastq> under $err_in[[ITER]] of <err.rds>
# (dada2-rs reads the same matrix via export_err_in.R) and print TRACE lines for
# the read: every cluster birth, the read's comparison with each new cluster
# (lambda, stored or not, lock, E_minmax), and its cluster after each shuffle
# (S) and p-update (P). Doubles are raw IEEE-754 bits, so equal means equal.
# The read is found by sequence on each side, so derep order need not match.
set -euo pipefail
BIN=${1:?binary}; ERR=${2:?err.rds}; FQ=${3:?fastq}; SEQ=${4:?sequence}; OUT=${5:?out-dir}
: "${R_TRACE_LIB:?set R_TRACE_LIB to the library holding the traced DADA2 (install_r_trace.sh)}"
HERE=$(cd "$(dirname "$0")" && pwd)
OMEGA_C=${OMEGA_C:-0}; BAND=${BAND:-32}; KMER=${KMER:-5}; THREADS=${THREADS:-1}
mkdir -p "$OUT"
printf '%s\n' "$SEQ" > "$OUT/sequence.txt"

Rscript "$HERE/export_err_in.R" "$ERR" "$OUT/models" > /dev/null
n=$(ls "$OUT"/models/err_in_*.json | wc -l | tr -d ' ')
ITER=${ITER:-$n}

# --- R ---
Rscript - "$ERR" "$FQ" "$OUT" "$ITER" "$OMEGA_C" "$BAND" "$THREADS" "$R_TRACE_LIB" <<'EOF' 2> "$OUT/r.log"
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[8], .libPaths()))
suppressPackageStartupMessages(library(dada2))
if (!startsWith(find.package("dada2"), normalizePath(a[8]))) stop("traced DADA2 not loaded from ", a[8])
e <- readRDS(a[1]); ins <- if (is.list(e$err_in)) e$err_in else list(e$err_in)
drp <- derepFastq(a[2])
seq <- readLines(file.path(a[3], "sequence.txt"))
idx <- match(seq, names(drp$uniques))
if (is.na(idx)) stop("sequence not in this FASTQ's uniques")
Sys.setenv(DADA2RS_TRACE_RAW = idx - 1L)
message(sprintf("R: dada2 %s from %s; read at derep index %d (%d reads)",
                packageVersion("dada2"), find.package("dada2"), idx - 1L, drp$uniques[idx]))
dd <- dada(drp, err = ins[[as.integer(a[4])]], OMEGA_C = as.numeric(a[5]),
           BAND_SIZE = as.integer(a[6]), multithread = as.integer(a[7]), verbose = FALSE)
message(sprintf("R: read ends in cluster with centre index %d",
                match(dd$sequence[dd$map[idx]], names(drp$uniques)) - 1L))
EOF

# --- dada2-rs ---
"$BIN" derep "$FQ" -o "$OUT/derep.json" 2> /dev/null
idx=$(python3 -c '
import json, sys
u = json.load(open(sys.argv[1]))["uniques"]
s = open(sys.argv[2]).read().strip()
print(next(i for i, x in enumerate(u) if x["sequence"] == s))' "$OUT/derep.json" "$OUT/sequence.txt")
DADA2RS_TRACE_RAW=$idx "$BIN" dada "$OUT/derep.json" --error-model "$OUT/models/err_in_$ITER.json" \
    --band "$BAND" --kmer-size "$KMER" --omega-c "$OMEGA_C" --threads "$THREADS" \
    -o "$OUT/rs.json" 2> "$OUT/rs.log"
echo "dada2-rs: $("$BIN" --version); read at derep index $idx" >> "$OUT/rs.log"

grep '^TRACE' "$OUT/r.log" > "$OUT/r.trace" || true
grep '^TRACE' "$OUT/rs.log" > "$OUT/rs.trace" || true
grep -v '^TRACE' "$OUT/r.log" | grep '^R:' || true
tail -n 1 "$OUT/rs.log"
echo "pass $ITER: R $(wc -l < "$OUT/r.trace") lines, dada2-rs $(wc -l < "$OUT/rs.trace") lines"
if cmp -s "$OUT/r.trace" "$OUT/rs.trace"; then
  echo "==> traces IDENTICAL"
else
  echo "==> first difference (R left, dada2-rs right):"
  diff "$OUT/r.trace" "$OUT/rs.trace" | head -n 12
fi
