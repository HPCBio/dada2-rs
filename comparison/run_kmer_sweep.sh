#!/usr/bin/env bash
# run_kmer_sweep.sh
# ---------------------------------------------------------------------------
# Sweep the k-mer pre-alignment screen size (--kmer-size) and measure its
# effect on (a) the inferred ASVs, (b) how aggressively the screen rejects
# pairs ("shrouding"), and (c) runtime / iterations. Built for issue #15.
#
# The pre-alignment screen rejects candidate pairs whose k-mer distance
# exceeds --kdist-cutoff, so only similar sequences reach the expensive
# Needleman-Wunsch alignment. On long reads (e.g. PacBio HiFi ~1.4 kb) the
# screen is effectively a no-op at the k=5 default (it shrouds ~nothing, so
# every pair is fully aligned); larger k makes it engage. The open question
# (issue #15) is whether a larger k changes the final ASVs vs. k=5 / k=6 —
# the hypothesis being that a larger k better isolates highly-similar
# sequences and may be an acceptable memory/speed tradeoff.
#
# This version accepts MULTIPLE samples (a directory or a list of files) and
# denoises them POOLED (R DADA2 pool=TRUE): all per-sample uniques are merged
# into one table and DADA2 runs once on the combined set. Pooling is the right
# mode for a large-scale ASV-impact test — it maximizes the unique count fed to
# the screen and surfaces rare/cross-sample variants, which is exactly where a
# larger k could change calls.
#
# Usage:
#   bash run_kmer_sweep.sh <input> [outdir] [errfun] [band] [k-list]
#
#   <input> may be:
#     - a directory      : all *.fastq / *.fastq.gz (or *.json / *.json.gz)
#                          inside are used as pooled samples
#     - a glob/file list  : quote it, e.g. "data/*.fastq.gz"
#     - a single file     : FASTQ or a derep/sample JSON
#
# Examples:
#   # All HiFi samples in a directory, pooled, default k list (5 6 8):
#   bash run_kmer_sweep.sh /path/to/fastq_dir ./ksweep_out pacbio 32
#
#   # Explicit file list, custom k list:
#   bash run_kmer_sweep.sh "data/*.filt.fastq.gz" ./ksweep_out pacbio 32 "5 6 7 8"
#
#   # Single Illumina sample:
#   bash run_kmer_sweep.sh F3D0_R1.fastq.gz ./ksweep_out loess 16 "5 6 7 8"
#
# Defaults: outdir=./kmer_sweep_out  errfun=pacbio  band=32  k-list="5 6 8"
#
# Environment overrides (held CONSTANT across all k):
#   MAX_CONSIST=10  KDIST_CUTOFF=0.42  THREADS=1  FILE_GLOB="*.fastq.gz *.fastq *.json *.json.gz"
# Set THREADS>1 to speed up large runs (wall times then are not directly
# comparable across machines, but the ASV/shroud results are unaffected).
#
# IMPORTANT for PacBio raw reads: reads must already be primer-trimmed AND
# uniformly oriented before this script. A sequence and its reverse complement
# share almost no k-mers, so mixed-orientation input defeats the screen at ANY
# k. If your reads are not oriented, run `remove-primers --orient` first.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DADA2RS="${SCRIPT_DIR}/../target/release-native/dada2-rs"

INPUT="${1:?Usage: run_kmer_sweep.sh <dir|glob|file> [outdir] [errfun] [band] [k-list]}"
OUTDIR="${2:-${SCRIPT_DIR}/kmer_sweep_out}"
ERRFUN="${3:-pacbio}"
BAND="${4:-32}"
KLIST="${5:-5 6 8}"

MAX_CONSIST="${MAX_CONSIST:-10}"
KDIST_CUTOFF="${KDIST_CUTOFF:-0.42}"
THREADS="${THREADS:-1}"
FILE_GLOB="${FILE_GLOB:-*.fastq.gz *.fastq *.json *.json.gz}"
# Bases used to LEARN the error model, matching R DADA2's learnErrors(nbases=1e8).
# Error learning runs on this subsample; denoising (dada/dada-pooled) always runs
# on the FULL data. Set NBASES=0 to disable subsampling (learn on everything).
NBASES="${NBASES:-100000000}"

if [[ ! -x "$DADA2RS" ]]; then
    echo "ERROR: binary not found/executable at $DADA2RS" >&2
    echo "       Run 'cargo build --release' first." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Step 0: resolve <input> into an array of sample files.
# ---------------------------------------------------------------------------
declare -a INPUTS=()
# Split FILE_GLOB into patterns with globbing OFF, so the pattern strings
# themselves are not expanded against the cwd (a nullglob footgun).
read -r -a GLOB_PATS <<< "$FILE_GLOB"
if [[ -d "$INPUT" ]]; then
    # Directory: collect by FILE_GLOB. nullglob makes non-matching patterns
    # expand to nothing instead of to the literal pattern.
    shopt -s nullglob
    for pat in "${GLOB_PATS[@]}"; do
        for f in "$INPUT"/$pat; do
            INPUTS+=("$f")
        done
    done
    shopt -u nullglob
elif [[ -f "$INPUT" ]]; then
    INPUTS+=("$INPUT")
else
    # Treat as a glob (already expanded by the caller, or expand here).
    shopt -s nullglob
    for f in $INPUT; do INPUTS+=("$f"); done
    shopt -u nullglob
fi
# Sort for deterministic, reproducible ordering across runs.
if [[ ${#INPUTS[@]} -gt 1 ]]; then
    IFS=$'\n' INPUTS=($(printf '%s\n' "${INPUTS[@]}" | sort)); unset IFS
fi

if [[ ${#INPUTS[@]} -eq 0 ]]; then
    echo "ERROR: no input files found for '$INPUT'" >&2
    echo "       (directory glob was: $FILE_GLOB)" >&2
    exit 1
fi

NSAMPLES=${#INPUTS[@]}
echo "==> ${NSAMPLES} sample(s):"
for f in "${INPUTS[@]}"; do echo "      $f"; done

# Paired-end footgun guard: pooling forward (R1) and reverse (R2) reads into one
# table is wrong — they cover different strands/regions and don't align, so the
# k-mer screen would shroud nearly everything and the ASVs would be meaningless.
# DADA2 denoises each orientation SEPARATELY. Warn (don't abort — some naming
# schemes legitimately contain "R1"/"R2" substrings) if both appear.
n_r1=0; n_r2=0
for f in "${INPUTS[@]}"; do
    b="$(basename "$f")"
    [[ "$b" == *_R1[._]* || "$b" == *_R1.* ]] && n_r1=$((n_r1+1))
    [[ "$b" == *_R2[._]* || "$b" == *_R2.* ]] && n_r2=$((n_r2+1))
done
if [[ $n_r1 -gt 0 && $n_r2 -gt 0 ]]; then
    echo "WARNING: input mixes forward (R1: ${n_r1}) and reverse (R2: ${n_r2}) reads." >&2
    echo "         DADA2 denoises each orientation separately; pooling them gives" >&2
    echo "         meaningless ASVs. Run the sweep once per orientation, e.g.:" >&2
    echo "           FILE_GLOB='*_R1_001.fastq.gz' bash $0 <dir> <out>_R1 ..." >&2
    echo "           FILE_GLOB='*_R2_001.fastq.gz' bash $0 <dir> <out>_R2 ..." >&2
    echo "         Continuing anyway in 3s (Ctrl-C to abort)..." >&2
    sleep 3
fi

POOLED=0
[[ $NSAMPLES -gt 1 ]] && POOLED=1

# ---------------------------------------------------------------------------
# Step 0b: build the JSON inputs. `errors-from-sample`, `dada`, and
# `dada-pooled` all take derep/sample JSON, not FASTQ. We build TWO sets:
#
#   DEREPS  — FULL dereplication of each sample (every unique). Used for
#             denoising (dada / dada-pooled): the ASV output we measure must
#             reflect the whole dataset.
#   LEARN   — SUBSAMPLE for error learning, matching R DADA2's
#             learnErrors(nbases=1e8): `sample --nbases NBASES`. Without this
#             the error fit (errors-from-sample) runs on every unique, doing
#             far more alignments than R would and inflating runtime.
#
# Both are derived once and reused across all k. JSON inputs are taken as-is
# for both sets (already dereplicated; subsampling a pre-derep'd JSON isn't
# supported, so pass FASTQ if you want learn-time subsampling).
# ---------------------------------------------------------------------------
DEREP_DIR="${OUTDIR}/derep"
LEARN_DIR="${OUTDIR}/learn_input"
mkdir -p "$DEREP_DIR" "$LEARN_DIR"
declare -a DEREPS=()   # full uniques, per sample (for denoising)
declare -a FASTQS=()   # FASTQ inputs only (for subsampled learning)
declare -a JSON_IN=()  # JSON inputs passed straight through

for f in "${INPUTS[@]}"; do
    case "$f" in
        *.json|*.json.gz)
            DEREPS+=("$f")
            JSON_IN+=("$f")
            ;;
        *)
            base="$(basename "$f")"
            base="${base%.gz}"; base="${base%.fastq}"
            dj="${DEREP_DIR}/${base}.derep.json"
            if [[ ! -s "$dj" ]]; then
                echo "==> derep (full) $f"
                "$DADA2RS" derep "$f" -o "$dj"
            fi
            DEREPS+=("$dj")
            FASTQS+=("$f")
            ;;
    esac
done

# Build the LEARN set: subsample FASTQ inputs to NBASES (R-style), pass any
# JSON inputs through unchanged. With NBASES=0, learn on the full derep set.
declare -a LEARN=()
if [[ "$NBASES" != "0" && ${#FASTQS[@]} -gt 0 ]]; then
    echo "==> subsampling ${#FASTQS[@]} FASTQ(s) to NBASES=${NBASES} for error learning"
    "$DADA2RS" sample "${FASTQS[@]}" --output-dir "$LEARN_DIR" \
        --nbases "$NBASES" --threads "$THREADS"
    # `sample` writes one JSON per input; collect them.
    shopt -s nullglob
    for j in "$LEARN_DIR"/*.json "$LEARN_DIR"/*.json.gz; do LEARN+=("$j"); done
    shopt -u nullglob
    # Plus any JSON inputs that bypassed subsampling (guard empty array under set -u).
    [[ ${#JSON_IN[@]} -gt 0 ]] && LEARN+=("${JSON_IN[@]}")
else
    [[ "$NBASES" == "0" ]] && echo "==> NBASES=0: learning errors on FULL data (no subsampling)"
    LEARN=("${DEREPS[@]}")
fi

if [[ ${#LEARN[@]} -eq 0 ]]; then
    echo "ERROR: no JSON inputs assembled for error learning." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 1: per-k learn-errors (pooled over all samples) + dada / dada-pooled.
# All samples and parameters are held constant; only --kmer-size changes.
# ---------------------------------------------------------------------------
# Wall-time + peak-RSS capture. Three tiers, picked once at startup:
#
#   Wall time  — ALWAYS captured via the bash clock (date +%s.%N / SECONDS); no
#                external binary needed (many clusters have only the bash `time`
#                keyword and no /usr/bin/time).
#   Peak RSS   — in priority order:
#                1. External timer: GNU `time -v` (its own pkg, NOT coreutils) or
#                   BSD `time -l` (macOS). Most accurate when present.
#                2. Linux /proc/<pid>/status VmHWM polling — no external binary,
#                   works on any Linux. VmHWM is the kernel's monotonic peak-RSS
#                   high-water mark, so polling and keeping the max catches it.
#                3. None (e.g. macOS without an external timer): RSS stays blank.
# The Step-2 parser understands all of these plus the WALL_S= line below.
EXT_TIME_BIN=""
EXT_TIME_FLAG=""
for cand in /usr/bin/time /opt/homebrew/bin/gtime /usr/local/bin/gtime "$(command -v gtime 2>/dev/null)"; do
    [[ -n "$cand" && -x "$cand" ]] || continue
    if "$cand" -v true >/dev/null 2>&1; then EXT_TIME_BIN="$cand"; EXT_TIME_FLAG="-v"; break; fi
    if "$cand" -l true >/dev/null 2>&1; then EXT_TIME_BIN="$cand"; EXT_TIME_FLAG="-l"; break; fi
done
# /proc VmHWM availability (Linux). Used only when no external timer was found.
HAVE_PROC_RSS=0
if [[ -z "$EXT_TIME_BIN" && -r /proc/self/status ]] && grep -q '^VmHWM:' /proc/self/status 2>/dev/null; then
    HAVE_PROC_RSS=1
fi
if [[ -n "$EXT_TIME_BIN" ]]; then
    : # external timer gives RSS
elif [[ $HAVE_PROC_RSS -eq 1 ]]; then
    echo "NOTE: no external timer; using /proc VmHWM polling for peak RSS." >&2
else
    echo "NOTE: no external timer and no /proc VmHWM; wall_s captured, maxrss_kb blank." >&2
fi

# High-resolution-ish wall clock: GNU `date +%s.%N` where available, else the
# integer bash SECONDS counter. (BSD date lacks %N and returns a literal "N",
# which we detect and fall back from.)
_now() {
    local t
    t=$(date +%s.%N 2>/dev/null)
    if [[ "$t" == *N* || -z "$t" ]]; then echo "$SECONDS"; else echo "$t"; fi
}

# Run "$@", capturing the command's stderr AND a wall-time line into $1.
# An external timer's report (incl. peak RSS) is appended when available;
# otherwise on Linux a PROC_MAXRSS_KB= line is appended via VmHWM polling.
# Usage: run_timed <logfile> <cmd...>
run_timed() {
    local tf="$1"; shift
    local start end rc
    start=$(_now)
    if [[ -n "$EXT_TIME_BIN" ]]; then
        rc=0; "$EXT_TIME_BIN" "$EXT_TIME_FLAG" "$@" >/dev/null 2>"$tf" || rc=$?
    elif [[ $HAVE_PROC_RSS -eq 1 ]]; then
        "$@" >/dev/null 2>"$tf" &
        local pid=$! peak=0 cur
        while kill -0 "$pid" 2>/dev/null; do
            cur=$(awk '/^VmHWM:/{print $2}' "/proc/$pid/status" 2>/dev/null || true)
            [[ -n "$cur" && "$cur" -gt "$peak" ]] && peak=$cur
            sleep 0.2
        done
        rc=0; wait "$pid" || rc=$?
        echo "PROC_MAXRSS_KB=$peak" >> "$tf"
    else
        rc=0; "$@" >/dev/null 2>"$tf" || rc=$?
    fi
    end=$(_now)
    # Append a parser-friendly wall-time line. awk handles float subtraction
    # portably (bash can't do float arithmetic).
    awk -v s="$start" -v e="$end" 'BEGIN{ printf "WALL_S=%.2f\n", e - s }' >> "$tf"
    return $rc
}

SUMMARY_CSV="${OUTDIR}/summary.csv"
echo "k,learn_iters,dada_aligns,dada_shrouded,shroud_pct,n_asv_total,wall_s_dada,maxrss_kb_dada" > "$SUMMARY_CSV"

for k in $KLIST; do
    echo ""
    echo "================  k = $k  ================"
    ERR_JSON="${OUTDIR}/errors_k${k}.json"
    LEARN_LOG="${OUTDIR}/learn_k${k}.log"
    DADA_OUTDIR="${OUTDIR}/dada_k${k}"
    # Single combined file: dada's --verbose stderr (the "ALIGN:" line) AND the
    # timer report both land here, so Step 2 reads aligns/shrouded and wall/RSS
    # from the one file. (run_timed sends both to this path.)
    DADA_LOG="${OUTDIR}/dada_k${k}.log"
    mkdir -p "$DADA_OUTDIR"

    echo "==> learn-errors (k=$k, ${#LEARN[@]} subsampled input(s), NBASES=${NBASES})"
    "$DADA2RS" errors-from-sample "${LEARN[@]}" \
        --errfun "$ERRFUN" --band "$BAND" \
        --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
        --max-consist "$MAX_CONSIST" --threads "$THREADS" \
        --verbose -o "$ERR_JSON" 2> "$LEARN_LOG"

    if [[ $POOLED -eq 1 ]]; then
        echo "==> dada-pooled (k=$k)"
        run_timed "$DADA_LOG" "$DADA2RS" dada-pooled "${DEREPS[@]}" \
            --error-model "$ERR_JSON" --output-dir "$DADA_OUTDIR" \
            --band "$BAND" --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
            --threads "$THREADS" --verbose
    else
        echo "==> dada (k=$k, single sample)"
        run_timed "$DADA_LOG" "$DADA2RS" dada "${DEREPS[0]}" \
            --error-model "$ERR_JSON" \
            --band "$BAND" --kmer-size "$k" --kdist-cutoff "$KDIST_CUTOFF" \
            --threads "$THREADS" --verbose -o "${DADA_OUTDIR}/sample.dada.json"
    fi
    echo "    errors -> $ERR_JSON"
    echo "    dada   -> $DADA_OUTDIR/"
done

# ---------------------------------------------------------------------------
# Step 2: analyze — fill the CSV and compare the POOLED ASV set across k.
# The ASV set is the UNION of ASV sequences across all per-sample dada outputs
# (pooled mode writes one JSON per sample, each containing only that sample's
# ASVs; the union is the full pooled ASV catalogue).
# ---------------------------------------------------------------------------
python3 - "$OUTDIR" "$SUMMARY_CSV" "$KLIST" <<'PY'
import json, os, re, sys, glob

outdir, csv_path, klist = sys.argv[1], sys.argv[2], sys.argv[3].split()

def asvs_from_file(path):
    """Set of ASV sequences from one dada output JSON (schema-tolerant)."""
    try:
        d = json.load(open(path))
    except Exception:
        return set()
    if isinstance(d, dict) and isinstance(d.get("asvs"), list):
        return {el["sequence"].upper() for el in d["asvs"]
                if isinstance(el, dict) and isinstance(el.get("sequence"), str)}
    # Fallback: largest list of dicts carrying a sequence-like field.
    best = set()
    def walk(node):
        nonlocal best
        if isinstance(node, list):
            seqs = set()
            for el in node:
                if isinstance(el, dict):
                    for key in ("sequence", "seq", "denoised", "asv"):
                        v = el.get(key)
                        if isinstance(v, str) and set(v.upper()) <= set("ACGTN-"):
                            seqs.add(v.upper()); break
            if len(seqs) > len(best): best = seqs
            for el in node: walk(el)
        elif isinstance(node, dict):
            for v in node.values(): walk(v)
    walk(d)
    return best

def pooled_asvs(k):
    """Union of ASVs across every per-sample JSON written for this k."""
    ddir = f"{outdir}/dada_k{k}"
    union = set()
    for path in glob.glob(f"{ddir}/*.json"):
        union |= asvs_from_file(path)
    return union

def parse_iters(log):
    if not os.path.exists(log): return ""
    converged = None; rounds = set()
    for ln in open(log, errors="replace"):
        m = re.search(r"converged after\s+(\d+)\s+iteration", ln)
        if m: converged = int(m.group(1)); continue
        m = re.search(r"\biter=(\d+)", ln)
        if m: rounds.add(int(m.group(1))); continue
        m = re.search(r"Selfconsist round\s+(\d+)", ln) or re.search(r"Iteration\s+(\d+)", ln)
        if m: rounds.add(int(m.group(1)))
    if converged is not None: return converged
    return max(rounds) if rounds else 0

def parse_align(log):
    """Sum (aligns, shrouded) over all 'ALIGN:' lines (one per sample in a
    pooled run is possible; sum gives the total screen workload)."""
    if not os.path.exists(log): return (0, 0)
    a = s = 0
    for ln in open(log, errors="replace"):
        m = re.search(r"ALIGN:\s*(\d+)\s+aligns,\s*(\d+)\s+shrouded", ln)
        if m: a += int(m.group(1)); s += int(m.group(2))
    return (a, s)

def parse_time(tf):
    if not os.path.exists(tf): return ("", "")
    wall = rss = ""
    txt = open(tf, errors="replace").read()
    # Wall time: prefer the always-present "WALL_S=" line that run_timed appends
    # (bash builtin clock, no external binary needed). Fall back to an external
    # timer's report if for some reason the line is missing.
    m = re.search(r"^WALL_S=([\d.]+)", txt, re.M)
    if m:
        wall = m.group(1)
    else:
        m = re.search(r"^\s*([\d.]+)\s+real", txt, re.M)   # BSD: "3.14 real"
        if m:
            wall = m.group(1)
        else:
            # GNU: "Elapsed (wall clock) time (h:mm:ss or m:ss): [h:]mm:ss[.ss]"
            m = re.search(r"wall clock.*?:\s*([\d:.]+)\s*$", txt, re.M)
            if m:
                secs = 0.0
                for p in m.group(1).split(":"):   # h:m:s or m:s or s
                    secs = secs * 60 + float(p)
                wall = f"{secs:g}"
    # Peak RSS (kB), in priority order: BSD time (bytes), GNU time (kB), or the
    # Linux /proc VmHWM polling line appended by run_timed (kB).
    m = re.search(r"(\d+)\s+maximum resident set size", txt)   # BSD: bytes
    if m:
        rss = str(round(int(m.group(1))/1024))
    else:
        m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)  # GNU: kB
        if m:
            rss = m.group(1)
        else:
            m = re.search(r"^PROC_MAXRSS_KB=(\d+)", txt, re.M)               # /proc VmHWM: kB
            if m and int(m.group(1)) > 0:
                rss = m.group(1)
    return (wall, rss)

asv_sets = {}; rows = []
for k in klist:
    iters = parse_iters(f"{outdir}/learn_k{k}.log")
    aligns, shrouded = parse_align(f"{outdir}/dada_k{k}.log")
    wall, rss = parse_time(f"{outdir}/dada_k{k}.log")
    asvs = pooled_asvs(k)
    asv_sets[k] = asvs
    pct = f"{100*shrouded/(aligns+shrouded):.1f}" if (aligns+shrouded) else ""
    rows.append((k, iters, aligns, shrouded, pct, len(asvs), wall, rss))

with open(csv_path, "w") as fh:
    fh.write("k,learn_iters,dada_aligns,dada_shrouded,shroud_pct,n_asv_total,wall_s_dada,maxrss_kb_dada\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")

print("\n================  SUMMARY  ================")
hdr = f"{'k':>3} {'iters':>5} {'aligns':>10} {'shroud':>10} {'shroud%':>7} {'#ASV':>6} {'wall_s':>8} {'RSS_MB':>7}"
print(hdr); print("-"*len(hdr))
for k, iters, aligns, shrouded, pct, n_asv, wall, rss in rows:
    rss_mb = f"{int(rss)/1024:.0f}" if rss else ""
    print(f"{k:>3} {iters:>5} {aligns:>10} {shrouded:>10} {pct:>7} {n_asv:>6} {wall:>8} {rss_mb:>7}")

ks = list(klist); base = ks[0]
b = asv_sets.get(base) or set()
print(f"\nPooled ASV-set comparison (baseline k={base}, {len(b)} ASVs):")
for k in ks[1:]:
    s = asv_sets.get(k) or set()
    only_base = len(b - s); only_k = len(s - b); shared = len(b & s)
    status = "IDENTICAL" if (only_base == 0 and only_k == 0) else "DIFFERS"
    print(f"  k={k}: shared={shared}  only_in_k{base}={only_base}  only_in_k{k}={only_k}   [{status}]")

print(f"\nWrote {csv_path}")
print("Per-k outputs: errors_k*.json, dada_k*/ (one JSON per sample), learn_k*.log, dada_k*.log")
PY

echo ""
echo "Done. CSV: ${SUMMARY_CSV}"
