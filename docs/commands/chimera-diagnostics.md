# `chimera-diagnostics`

Screen a sequence table for higher-order chimeras (trimeras).

```bash
dada2-rs chimera-diagnostics seqtab.json --threads 24 -o trimera.tsv
```

Reads a [`make-sequence-table`](make-sequence-table.md) or
[`remove-bimera-denovo`](remove-bimera-denovo.md) JSON and emits a per-sequence
TSV of bimera **coverage** metrics.

Unlike the boolean `remove-bimera-denovo` decision, this retains how much of
each sequence a single two-parent junction explains. A sequence that survives
bimera removal yet is nearly covered — high `cover_frac` — leaves a small
internal gap that a third parent could fill, making it a trimera suspect.

This is most useful on long amplicons (full-length 16S, nodA) and low-biomass
samples, where complex chimeras are more likely. Coverage uses the pooled,
across-sample abundance model.

!!! warning "Over-calls few-SNP variants"
    On real data this screen flags close variants of one abundant parent, which
    leave a coverage gap without being chimeric at all. The
    `--trimera-min-parent-dist` gate below exists to reject them; treat every
    hit as a suspect to investigate, never as a call.

## Input

**`<INPUT>`** — sequence table JSON.

## Chimera

**`--min-fold-parent-over-abundance`** (default 1.5) and
**`--min-parent-abundance`** (default 2) — parent eligibility, as in
[`remove-bimera-denovo`](remove-bimera-denovo.md#chimera).

**`--trimera-min-parent-dist`** (default 15) — minimum distance to the nearest
single parent before a sequence is flagged. This is the gate that rejects
few-SNP variants of one abundant parent.

**`--trimera-min-gap`** (default 20) — minimum residual gap length in bp for a
credible third segment. Rejects one-off bimeras, whose gap is ~1 base.

**`--trimera-max-gap-error`** (default 0.10) — maximum third-parent mismatch
fraction across the gap, i.e. how cleanly a third parent must fit.

**`--trimera-min-flank`** (default 30) — minimum length in bp of each end flank.
A genuine three-segment mosaic needs two substantial flanks; this rejects
tiny-flank divergent singletons.

## Alignment

**`--max-shift`** (default 16), **`--match`** (5), **`--mismatch`** (−4),
**`--gap-p`** (−8), **`--align-backend`** — as in
[`remove-bimera-denovo`](remove-bimera-denovo.md#alignment).

## Performance

**`--threads`** (default 1) — threads for parallel diagnostics.

## Output

**`--output` / `-o`** — write the TSV here instead of stdout.

## Experimental

**`--wfa-max-edits`** (default 50, `0` = unbounded) — see
[`dada`](dada.md#experimental).
