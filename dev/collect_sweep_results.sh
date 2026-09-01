#!/usr/bin/env bash
# collect_sweep_results.sh — package a screen-sweep run into a small bundle.
#
# A sweep output directory is mostly bulk that does not need to travel: each arm
# carries filtered FASTQs, per-sample `dada` JSONs (full read->ASV maps), merged
# reads and an intermediate seqtab. Of that, the analysis needs only:
#
#   * seqtab.nochim.json per arm   -- ASV set + the sample x ASV count matrix
#   * two integers per dada JSON   -- nalign / nshroud, for alignment work
#   * timings.tsv                  -- wall clock, if the sweep timed anything
#   * phase_split.txt              -- screen vs align share, which says whether a
#                                     timing result was screen- or align-driven
#   * errF.json / err.json per arm -- error models, IF they differ between arms
#
# The per-sample dada JSONs are the bulk and are reduced here to one TSV row
# each, cluster-side, so gigabytes become kilobytes.
#
# Usage:
#   collect_sweep_results.sh <sweep-out-dir> [bundle.tar.gz]
#
# Then copy the bundle over and unpack it anywhere.
set -euo pipefail

SRC="${1:?usage: collect_sweep_results.sh <sweep-out-dir> [bundle.tar.gz]}"
BUNDLE="${2:-$(basename "$SRC")-bundle.tar.gz}"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
NAME="$(basename "$SRC")"
mkdir -p "$STAGE/$NAME"

shopt -s nullglob
echo "==> collecting from $SRC"

# Where the arms live. `run_screen_sweep.sh` puts them at the top level; the
# PacBio script nests them under `arms/`. Collecting only the top level is how
# the first pooled PacBio bundle arrived with an EMPTY align_stats.tsv and no
# seqtab at all -- the timing and phase-split halves were fine, so nothing looked
# broken until the accuracy table was missing. Accept both layouts.
ARMS=("$SRC"/*/ "$SRC"/arms/*/)
is_arm() {  # a directory of derep/model bulk is not an arm
  case "$(basename "$1")" in derep|models|arms|.timing|.verbose) return 1 ;; esac
  return 0
}

# One row per (arm, sample): the alignment-work numbers, extracted here so the
# multi-GB dada JSONs stay put.
{
  printf "arm\tstage\tsample\tnalign\tnshroud\taligned\tn_asvs\ttotal_reads\n"
  for arm in "${ARMS[@]}"; do
    is_arm "$arm" || continue
    a=$(basename "$arm")
    # "." = the arm directory itself: a per-sample `dada` sweep writes its
    # sample JSONs straight into the arm, with no `dada/` stage subdirectory.
    for stage in . dada dada_fwd dada_rev; do
      for j in "$arm$stage"/*.json; do
        python3 - "$a" "$stage" "$j" <<'PY'
import json, os, sys
arm, stage, path = sys.argv[1], sys.argv[2], sys.argv[3]
j = json.load(open(path))
s = j.get("stats")
# Scanning the arm directory itself also sees seqtab / error-model JSONs. Those
# have no `stats`, and emitting them as zero rows would quietly pad the table.
if not isinstance(s, dict) or "nalign" not in s:
    sys.exit(0)
na, ns = s.get("nalign", 0), s.get("nshroud", 0)
print(f"{arm}\t{stage}\t{j.get('sample', os.path.basename(path))}\t{na}\t{ns}\t{na-ns}"
      f"\t{len(j.get('asvs', []))}\t{j.get('total_reads', 0)}")
PY
      done
    done
  done
} > "$STAGE/$NAME/align_stats.tsv"
echo "    align_stats.tsv: $(( $(wc -l < "$STAGE/$NAME/align_stats.tsv") - 1 )) rows"

# The chimera-filtered tables -- the actual comparison inputs.
n=0
for arm in "${ARMS[@]}"; do
  is_arm "$arm" || continue
  a=$(basename "$arm")
  for f in "$arm/seqtab.nochim.json"; do
    [ -f "$f" ] || continue
    mkdir -p "$STAGE/$NAME/$a"
    cp "$f" "$STAGE/$NAME/$a/"
    n=$((n+1))
  done
done
echo "    seqtab.nochim.json: $n arms"

# Error models. Small, and needed to check whether the arms actually shared one
# -- if they did not, alignment counts between arms are not comparable.
for arm in "${ARMS[@]}"; do
  is_arm "$arm" || continue
  a=$(basename "$arm")
  for f in "$arm"/err*.json; do
    [ -f "$f" ] || continue
    mkdir -p "$STAGE/$NAME/$a"
    cp "$f" "$STAGE/$NAME/$a/"
  done
done
# Models and the calibration SUMMARY -- but not the calibration CSVs themselves.
# `cal_*.csv` is one row per sampled pair: 2.3 GB on the PacBio run, which is why
# the first bundle was 545 MB of which ~1 MB was needed. The note at the end of
# this script always said the CSVs were excluded; the `models/*.csv` glob was
# quietly including them anyway.
for f in "$SRC"/models/*.json "$SRC"/models/*.txt; do
  [ -f "$f" ] || continue
  mkdir -p "$STAGE/$NAME/models"
  cp "$f" "$STAGE/$NAME/models/"
done
for f in "$SRC"/models/*.csv; do
  [ -f "$f" ] || continue
  echo "    skipping $(basename "$f") ($(du -h "$f" | cut -f1); one row per sampled pair)"
done

# Timings and logs (logs are small and carry the parameters each arm ran with).
for f in "$SRC"/timings.tsv "$SRC"/phase_split.txt "$SRC"/*.log; do
  [ -f "$f" ] && cp "$f" "$STAGE/$NAME/"
done

tar czf "$BUNDLE" -C "$STAGE" "$NAME"
echo "==> wrote $BUNDLE ($(du -h "$BUNDLE" | cut -f1))"
echo
echo "Contents:"
tar tzf "$BUNDLE" | head -25
echo
echo "NOTE: kdist-calibrate CSVs are NOT included -- they are one row per sampled"
echo "pair and can be very large. Either gzip and send them separately, or run"
echo "  python3 dev/analyze_kdist_curves.py kmer=k.csv minimizer=m.csv > cal.txt"
echo "cluster-side and send cal.txt instead."
