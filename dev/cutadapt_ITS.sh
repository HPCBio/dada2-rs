#!/usr/bin/env bash
# Primer removal for gITS7-ITS4 (ITS2 region) fungal amplicon data.
#
# Written for a NovaSeq 6000 soil dataset (BioProject PRJNA1504839) whose NCBI
# metadata claims "ITS1FI2-ITS2 (ITS2)". VERIFIED AGAINST THE READS, and the
# primer names are wrong -- the region label is not:
#
#   R1: gITS7 GTGARTCATCGARTCTTTG  92.5%, offset 0, no pad
#       (ITS86F is the non-degenerate form of the same primer and matches the
#        49.3% subset where the degenerate base is A; ITS1F/ITS1FI2/ITS1 = 0%)
#   R2: ITS4  TCCTCCGCTTATTGATATGC 90.9%, behind a CONSTANT 4-nt pad "GATA"
#
# gITS7+ITS4 amplifies the ITS2 region, so "(ITS2)" is right while
# "ITS1FI2-ITS2" is not. Note "ITS2" names both a region and a primer;
# ITS1FI2+ITS2-the-primer would amplify ITS1, contradicting the region label.
#
# Unlike the 16S arm of this BioProject, the offsets here are CONSTANT (R1 at 0,
# R2 at 4), so a fixed --trim-left would in principle work. cutadapt is still
# preferred: it also removes 3' read-through, which matters for ITS (below).
#
# Why cutadapt and not `dada2-rs filter-and-trim`: filter-and-trim removes
# leading bases by fixed offset (--trim-left) and does not match primer
# sequence. With the primer at 5 different offsets no single value is correct
# (see issue #113). Primer removal must be sequence-aware and happen here;
# filter-and-trim is for quality filtering afterwards.
#
# Usage:
#   ./cutadapt_ITS.sh <in_dir> <out_dir> [threads]
#
# Expects paired files named <sample>_1.fastq.gz / <sample>_2.fastq.gz
# (SRA fasterq-dump convention). Adjust SUF1/SUF2 below if yours differ.

set -euo pipefail

# Cluster module load; skipped silently off-cluster so the script stays testable.
if type module &>/dev/null; then
    module load cutadapt/3.7-IGB-gcc-8.2.0-Python-3.7.2
fi

IN_DIR=${1:?usage: cutadapt_ITS.sh <in_dir> <out_dir> [threads]}
OUT_DIR=${2:?usage: cutadapt_ITS.sh <in_dir> <out_dir> [threads]}
THREADS=${3:-8}

SUF1="_1.fastq.gz"
SUF2="_2.fastq.gz"

# gITS7 / ITS4 primers (IUPAC degenerate; cutadapt understands these natively)
FWD="GTGARTCATCGARTCTTTG"      # gITS7 (5.8S)
REV="TCCTCCGCTTATTGATATGC"     # ITS4  (LSU)
# Reverse complements, for 3' read-through removal
FWD_RC="CAAAGAYTCGATGAYTCAC"
REV_RC="GCATATCAATAAGCGGAGGA"

ERROR_RATE=0.15   # tolerant: degenerate primers + spacers
MIN_LEN=100       # drop reads that lose too much

mkdir -p "$OUT_DIR"
LOG_DIR="$OUT_DIR/logs"; mkdir -p "$LOG_DIR"

command -v cutadapt >/dev/null || { echo "ERROR: cutadapt not on PATH" >&2; exit 1; }
echo "cutadapt: $(cutadapt --version)"

# ---------------------------------------------------------------------------
# Step 1 — DIAGNOSE before trimming.
#
# Do not assume both primers are present, or that the metadata names the right
# ones -- on this BioProject it did not. Measure first. (Both primers ARE
# retained here, at ~92% / ~91%, so DISCARD_UNTRIMMED=1 is safe for this data.)
# ---------------------------------------------------------------------------
# One pass per file, counting both primers and the ITS1f offset distribution.
#
# NOTE: `head` closing the pipe makes zcat exit on SIGPIPE, which `set -o
# pipefail` turns into a script-killing failure. pipefail is disabled around
# the read and restored afterwards -- do not "simplify" this away.
probe() {  # probe <fastq.gz> <label>
    local f=$1 label=$2
    set +o pipefail
    zcat -f -- "$f" 2>/dev/null | head -n 200000 | awk -v label="$label" '
        NR % 4 == 2 {
            n++
            if (match($0, /GTGA[AG]TCATCGA[AG]TCTTTG/) && RSTART <= 9) {
                gits7++; off[RSTART - 1]++
            }
            if (match($0, /TCCTCCGCTTATTGATATGC/) && RSTART <= 9) {
                its4++; roff[RSTART - 1]++
            }
        }
        END {
            if (n == 0) { printf "  %s: NO READS READ (bad path or not gzipped?)\n", label; exit }
            printf "  %s: %d reads sampled\n", label, n
            printf "     gITS7 %7d  (%5.1f%%)\n", gits7, 100 * gits7 / n
            printf "     ITS4  %7d  (%5.1f%%)\n", its4, 100 * its4 / n
            if (gits7 > 0) {
                printf "     gITS7 start offsets:"
                for (o = 0; o <= 9; o++)
                    if (o in off) printf " %d:%.0f%%", o, 100 * off[o] / gits7
                printf "\n"
            }
            if (its4 > 0) {
                printf "     ITS4  start offsets:"
                for (o = 0; o <= 9; o++)
                    if (o in roff) printf " %d:%.0f%%", o, 100 * roff[o] / its4
                printf "\n"
            }
        }'
    local rc=${PIPESTATUS[2]}
    set -o pipefail
    return "$rc"
}

FIRST_R1=$(find "$IN_DIR" -name "*$SUF1" | sort | head -1)
[[ -n "$FIRST_R1" ]] || { echo "ERROR: no *$SUF1 files under $IN_DIR" >&2; exit 1; }
FIRST_R2=${FIRST_R1%$SUF1}$SUF2
[[ -f "$FIRST_R2" ]] || { echo "ERROR: no mate $FIRST_R2" >&2; exit 1; }

echo "=== Primer diagnostic on $(basename "${FIRST_R1%$SUF1}") (first 50k reads)"
probe "$FIRST_R1" "R1"
probe "$FIRST_R2" "R2"
echo
echo "  Interpretation:"
echo "    gITS7 high on R1, ITS4 high on R2  -> primers retained on both reads"
echo "                                         (EXPECTED here); use DISCARD_UNTRIMMED=1"
echo "    gITS7 high on R1, ITS4 ~0 on R2    -> asymmetric prep; leave the default"
echo "    both near zero                     -> already trimmed; skip this script"
echo "    gITS7 high on R2 instead of R1     -> R1/R2 are swapped; fix SUF1/SUF2"
echo
echo "  A SINGLE start offset means a constant pad, which a fixed --trim-left"
echo "  could remove. SEVERAL offsets mean heterogeneity spacers, and then no"
echo "  fixed --trim-left value is correct (see #113). Expected here: R1 at 0,"
echo "  R2 at 4 (a constant GATA pad)."
echo

# Set to 1 ONLY if the diagnostic shows both primers present at high rates.
DISCARD_UNTRIMMED=${DISCARD_UNTRIMMED:-0}
if [[ "$DISCARD_UNTRIMMED" == "1" ]]; then
    DISCARD_FLAG="--discard-untrimmed"
    echo ">>> --discard-untrimmed ENABLED"
else
    DISCARD_FLAG=""
    echo ">>> --discard-untrimmed disabled (reads without a primer are kept, untrimmed)"
fi
echo

# ---------------------------------------------------------------------------
# Step 2 — trim.
#
# -g / -G are NON-anchored 5' adapters: cutadapt removes the primer *and*
# anything preceding it, which is exactly what clears the 0-4 nt spacer.
# -a / -A remove read-through into the opposite primer at the 3' end.
# ---------------------------------------------------------------------------
shopt -s nullglob
n=0
for R1 in "$IN_DIR"/*"$SUF1"; do
    R2=${R1%$SUF1}$SUF2
    S=$(basename "$R1" "$SUF1")
    [[ -f "$R2" ]] || { echo "WARN: no mate for $S, skipping" >&2; continue; }

    cutadapt \
        -g "$FWD" -G "$REV" \
        -a "$REV_RC" -A "$FWD_RC" \
        -e "$ERROR_RATE" -n 2 --max-n 0 \
        --minimum-length "$MIN_LEN" \
        $DISCARD_FLAG \
        -j "$THREADS" \
        -o "$OUT_DIR/${S}${SUF1}" \
        -p "$OUT_DIR/${S}${SUF2}" \
        "$R1" "$R2" \
        > "$LOG_DIR/${S}.cutadapt.log" 2>&1

    n=$((n+1))
    printf '  [%3d] %-24s %s\n' "$n" "$S" \
        "$(grep -m1 'Pairs written' "$LOG_DIR/${S}.cutadapt.log" | sed 's/^ *//')"
done
shopt -u nullglob

echo
echo "=== Done: $n pair(s) -> $OUT_DIR"
echo
echo "Summary of retention:"
grep -h 'Pairs written' "$LOG_DIR"/*.cutadapt.log \
  | sed -E 's/.*\(([0-9.]+)%\).*/\1/' \
  | awk '{s+=$1; n++; if($1<min||n==1)min=$1; if($1>max)max=$1}
         END{printf "  mean %.1f%%  min %.1f%%  max %.1f%%  over %d samples\n", s/n, min, max, n}'

cat <<'NOTE'

Next steps
----------
1. Check the retention summary above. A sharp drop means the primer/orientation
   assumption is wrong -- re-read the diagnostic rather than lowering -e.

2. DO NOT use fixed-length truncation on ITS. The ITS2 region varies in length
   by hundreds of nt ACROSS TAXA, so any --trunc-len is a taxonomic filter: it
   silently deletes every fungus whose amplicon is longer than the cutoff, and
   the loss is systematic rather than random. This is the one place ITS differs
   hardest from 16S, where a fixed insert makes truncation safe. Quality-filter
   with maxEE only.

   Same reasoning applies to the merge step: long-ITS taxa may fail to overlap
   at all on 2x250 and drop out. Compare the taxonomic composition of merged
   vs unmerged reads before trusting the table -- an ITS merge rate is not
   taxonomically neutral the way a 16S one is.

3. Read-through is the mirror image of the same problem. SHORT ITS2 amplicons
   (< read length) run past the far primer, so -a/-A above are load-bearing
   here, not belt-and-braces. Measured before trimming on this data it is only
   ~0.1%, but that is a property of this community, not of the assay.

4. Then quality-filter (NOT primer-trim) with filter-and-trim, and run
   learn-errors + dada-pooled for both errfun arms to rebuild the binned-quality
   A/B on the ITS amplicon.
NOTE
