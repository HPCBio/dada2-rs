# train_pin.sh -- error-model training-set pinning for the concordance runners.
# Sourced by run_illumina.sh and run_pacbio.sh; not executable on its own.
#
# A pinned manifest (TRAIN_SAMPLES: one sample name per line, the contract of
# write_reference.R --train-samples and pick_training_subset.py) makes both
# tools learn from exactly the same samples. NBASES must then cover the whole
# set: learn-errors takes whole samples in order until the budget is passed, so
# a smaller NBASES silently trains on a prefix of the manifest (#225; the trap
# write_reference.R's assert_nbases_covers() guards on the R side).
#
# Plain bash 3.2 (macOS /bin/bash): no namerefs or associative arrays.

# train_pin_nbases <value>: print NBASES as an integer. learn-errors parses
# --nbases as u64 and rejects "1e12"; R and pick_training_subset.py accept it.
train_pin_nbases() {
  awk -v v="$1" 'BEGIN {
    if (v !~ /^[0-9]+(\.[0-9]+)?([eE][+]?[0-9]+)?$/) exit 1
    printf "%.0f\n", v
  }' || { echo "NBASES is not a number: $1" >&2; return 1; }
}

# train_pin_indices <caller> <manifest> <name>...: print the 0-based index of
# each manifest name among <name>..., one per line, in manifest order. Blank
# lines and CR are ignored. Fails on a name that is not present.
train_pin_indices() {
  local caller=$1 manifest=$2 want name i found
  shift 2
  [ -r "$manifest" ] || { echo "$caller: cannot read TRAIN_SAMPLES manifest $manifest" >&2; return 1; }
  while IFS= read -r want || [ -n "$want" ]; do
    want=$(printf '%s' "$want" | tr -d '[:space:]')
    [ -z "$want" ] && continue
    found=""
    i=0
    for name in "$@"; do
      if [ "$name" = "$want" ]; then found=$i; break; fi
      i=$((i + 1))
    done
    if [ -z "$found" ]; then
      echo "$caller: TRAIN_SAMPLES name not present in the inputs: $want" >&2
      return 1
    fi
    echo "$found"
  done < "$manifest"
}

# train_pin_guard <caller> <nbases> <label> <fastq.gz>...: exit if the files
# hold more bases than <nbases>. Run once per read direction on paired data;
# the two directions are learned separately and differ in length.
train_pin_guard() {
  local caller=$1 nbases=$2 label=$3 bases
  shift 3
  bases=$(gzip -cd "$@" | awk 'NR % 4 == 2 { n += length($0) } END { printf "%.0f\n", n }')
  echo "==> training set ($label): $# samples, $bases bases; nbases = $nbases"
  if awk -v b="$bases" -v n="$nbases" 'BEGIN { exit !(b > n) }'; then
    echo "$caller: NBASES ($nbases) is smaller than the pinned $label training set" >&2
    echo "  ($bases bases), so learn-errors would train on only a prefix." >&2
    echo "  Re-run with NBASES=$(awk -v b="$bases" 'BEGIN { printf "%.0f", b * 1.1 }') or larger." >&2
    exit 1
  fi
}
