#!/usr/bin/env bash
# Primer removal for V3-V4 (341F / 805R) amplicon data with heterogeneity spacers.
#
# Written for a NovaSeq 6000 soil dataset (BioProject PRJNA1504839) whose NCBI
# metadata claims "515F-806R (V4)" but which is demonstrably V3-V4:
#   - 341F present in 99.96% of merged reads at the 5' end, at 5 distinct
#     offsets (0-4 nt heterogeneity spacers)
#   - the 515F site sits ~150 nt INTERNALLY
#   - no reverse primer detected in the merged product
#
# Why cutadapt and not `dada2-rs filter-and-trim`: filter-and-trim removes
# leading bases by fixed offset (--trim-left) and does not match primer
# sequence. With the primer at 5 different offsets no single value is correct
# (see issue #113). Primer removal must be sequence-aware and happen here;
# filter-and-trim is for quality filtering afterwards.
#
# Usage:
#   ./cutadapt_v3v4.sh <in_dir> <out_dir> [threads]
#
# Expects paired files named <sample>_1.fastq.gz / <sample>_2.fastq.gz
# (SRA fasterq-dump convention). Adjust SUF1/SUF2 below if yours differ.

set -euo pipefail

module load cutadapt/3.7-IGB-gcc-8.2.0-Python-3.7.2

IN_DIR=${1:?usage: cutadapt_v3v4.sh <in_dir> <out_dir> [threads]}
OUT_DIR=${2:?usage: cutadapt_v3v4.sh <in_dir> <out_dir> [threads]}
THREADS=${3:-8}

SUF1="_1.fastq.gz"
SUF2="_2.fastq.gz"

# V3-V4 primers (IUPAC degenerate; cutadapt understands these natively)
FWD="CCTACGGGNGGCWGCAG"        # 341F
REV="GACTACHVGGGTATCTAATCC"    # 805R
# Reverse complements, for 3' read-through removal
FWD_RC="CTGCWGCCNCCCGTAGG"
REV_RC="GGATTAGATACCCBDGTAGTC"

ERROR_RATE=0.15   # tolerant: degenerate primers + spacers
MIN_LEN=100       # drop reads that lose too much

mkdir -p "$OUT_DIR"
LOG_DIR="$OUT_DIR/logs"; mkdir -p "$LOG_DIR"

command -v cutadapt >/dev/null || { echo "ERROR: cutadapt not on PATH" >&2; exit 1; }
echo "cutadapt: $(cutadapt --version)"

# ---------------------------------------------------------------------------
# Step 1 — DIAGNOSE before trimming.
#
# Do not assume both primers are present. In this dataset the forward primer
# was retained and the reverse was not; running with --discard-untrimmed on
# both reads would throw away nearly everything. Measure first.
# ---------------------------------------------------------------------------
probe() {  # probe <fastq.gz> <pattern-name> <iupac-regex>
    local f=$1 name=$2 rx=$3
    local n hit
    n=$(zcat < "$f" | head -n 200000 | awk 'NR%4==2' | wc -l)
    hit=$(zcat < "$f" | head -n 200000 | awk 'NR%4==2' | grep -cE "^.{0,8}$rx" || true)
    printf '    %-8s %6d / %6d  (%.1f%%)\n' "$name" "$hit" "$n" \
        "$(awk -v a="$hit" -v b="$n" 'BEGIN{print (b?100*a/b:0)}')"
}

FWD_RX="CCTACGGG.GGC[AT]GCAG"
REV_RX="GACTAC[ACT][ACG]GGGTATCTAATCC"

FIRST_R1=$(find "$IN_DIR" -name "*$SUF1" | sort | head -1)
FIRST_R2=${FIRST_R1%$SUF1}$SUF2
echo "=== Primer diagnostic on $(basename "$FIRST_R1") (first 50k reads)"
echo "  R1:"; probe "$FIRST_R1" "341F" "$FWD_RX"; probe "$FIRST_R1" "805R" "$REV_RX"
echo "  R2:"; probe "$FIRST_R2" "341F" "$FWD_RX"; probe "$FIRST_R2" "805R" "$REV_RX"
echo
echo "  Interpretation:"
echo "    both high      -> primers retained on both reads; set DISCARD_UNTRIMMED=1"
echo "    R1 only high   -> asymmetric (this dataset); leave DISCARD_UNTRIMMED=0"
echo "    both near zero -> already trimmed; skip this script entirely"
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
        -e "$ERROR_RATE" -n 2 \
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

2. Watch the overlap budget. This library is ~443 nt with 240 nt reads, i.e.
   ~37 nt of overlap BEFORE primer removal. Stripping 341F (17 nt + up to 4 nt
   spacer) from R1 leaves roughly 33 nt. That is still workable but thin, so
   keep truncation conservative in filter-and-trim -- aggressive --trunc-len
   will push pairs under minOverlap and silently cost you the merge.

3. Then quality-filter (NOT primer-trim) with filter-and-trim, and re-run
   learn-errors + dada-pooled for both errfun arms so the A/B is rebuilt on
   properly trimmed reads.
NOTE
