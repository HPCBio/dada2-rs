# `remove-bimera-denovo`

Remove bimeric sequences from a sequence table. Mirrors R's
`removeBimeraDenovo()`.

```bash
dada2-rs remove-bimera-denovo seqtab.json --threads 24 -o seqtab.nochim.json
```

Reads a [`make-sequence-table`](make-sequence-table.md) JSON and removes
sequences identified as bimeras — chimeras of two more-abundant parents.

## Input

**`<INPUT>`** — sequence table JSON.

## Chimera

**`--method`** (default `consensus`) — `consensus`, `pooled`, or `per-sample`.
`consensus` flags per sample and then votes, which is R's default and the safest
choice on multi-sample runs.

**`--min-fold-parent-over-abundance`** (default 1.5) — minimum fold-difference
in abundance for a sequence to be considered a parent.

**`--min-parent-abundance`** (default 2) — minimum abundance for a sequence to
be a parent.

**`--allow-one-off`** — also flag sequences one mismatch or indel away from an
exact bimera. Off by default, and it raises false positives on data with real
few-SNP variants.

**`--min-one-off-parent-distance`** (default 4) — minimum mismatches to a parent
required before one-off detection applies.

**`--min-sample-fraction`** (default 0.9) — `consensus` only: the fraction of
samples a sequence must be flagged in.

**`--ignore-n-negatives`** (default 1) — `consensus` only: how many unflagged
samples to ignore in the fraction vote.

## Alignment

**`--max-shift`** (default 16) — maximum shift in the ends-free alignment to
potential parents.

**`--match`** (5), **`--mismatch`** (−4), **`--gap-p`** (−8) — parent alignment
scoring. R's `removeBimeraDenovo` honors the global `setDadaOpt` gap penalty,
and dada2-rs matches that.

**`--align-backend`** — `nw` (default) or the experimental `wfa2`, which
requires a build with `--features wfa`.

## Performance

**`--threads`** (default 1) — threads for parallel bimera detection.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--compact`** — minified JSON.

## Diagnostics

**`--verbose`** — progress to stderr.

## Experimental

**`--wfa-max-edits`** (default 50, `0` = unbounded) — WFA edit-budget cap, used
with `--align-backend wfa2`. See [`dada`](dada.md#experimental).

## See also

- [`chimera-diagnostics`](chimera-diagnostics.md) — higher-order chimeras
