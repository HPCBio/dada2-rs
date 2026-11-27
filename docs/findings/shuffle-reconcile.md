# Inside `b_shuffle`: the reconcile's wasted passes

**Verdict:** two exact cuts to the shuffle's reconcile take NovaSeq
`dada-pooled` from 1,335 s to 1,082 s (**−19% wall**), with output identical to
`main` on all three test pools. Shuffle falls from 379 s to 144 s, from 30% of
the main loop to 14%. MiSeq gains 10% and PacBio 4%. Both cuts are on by
default.

- **Skip fallen clusters in pass A2** ([#358][358], `7c80b49`): −83 s on
  NovaSeq. A2's walk over a cluster whose reads fell cannot change any raw's
  best cluster, and that walk was 97.8% of A2's pairs.
- **Settle fallen incumbents by a carried runner-up bound** (`7e77acf`):
  another −152 s. Rescans fall from 2.51 billion to 815,000, and rescan time
  from 160 s to 0.2 s.

Changes since [the build-scan closure](shuffle-build-scan.md) halved
`b_compare` on NovaSeq, which made shuffle 29% of `run_dada` there. The
accounting showed that most of the reconcile's work proved things already
known.

## Background

The incremental `b_shuffle_converge` keeps `compmax` (each raw's best cluster)
and `emax` (its score, `λ × reads`) across buds ([#139][139]). After each move
pass, the reconcile restores them:

- **A1** walks the comp vecs of clusters whose reads fell and marks every raw
  whose incumbent is that cluster.
- **A2** walks the comp vecs of every changed cluster (fallen and risen) and
  reassigns a raw when that cluster now beats its incumbent under the
  `(score, −ci)` order.
- **Pass B** rescans each marked raw over its candidates via the raw-major
  index.

## The accounting

Measured first, with nothing built (instrumentation branch
`perf/shuffle-account-358`, job 33763, NovaSeq, 48 threads):

| part | time | work | what it changed |
|---|---|---|---|
| collect (A1 + A2) | 189 s | 29.9 G A2 pairs | A2 on fallen clusters: 29.3 G pairs, **0 effects** |
| rescan (pass B) | 160 s | 2.51 G raws, 2.0 candidates each | **99.72%** returned the same incumbent; 80.7% had one candidate |

[#358][358]'s 71 ns/comp divided the whole reconcile by the rescan's
comparisons. The two halves are about the same size, at ~6.3 ns per A2 pair
and ~65 ns per rescanned raw, and need different cuts.

## Cut 1: A2 skips fallen clusters

A fallen cluster's score for any raw only drops. If that cluster was the raw's
incumbent, A1 has already marked the raw for a rescan. Otherwise the incumbent
already beat the cluster under `(score, −ci)`, and a lower score cannot reverse
that. A2 therefore has nothing to do on fallen clusters. The accounting
confirmed it empirically, with 0 effects over 29.3 G pairs. The cut is
identical by construction.

## Cut 2: a carried runner-up bound

`emax` carries, beside each raw's best score, an upper bound on the score of
every other candidate. When A1 finds a fallen incumbent, it reprices it. If the
new score still strictly exceeds the bound, the raw keeps its incumbent and
skips the rescan. A tie goes to the rescan, which applies the `−ci`
tie-break.

The bound must stay an upper bound through every update:

- the build and pass B compute the exact runner-up;
- A2 folds in the old score of a displaced incumbent;
- A2 folds in the score of a risen cluster that did not win;
- an incumbent whose own cluster rose keeps its bound.

`DADA2RS_RECONCILE_VERIFY=1` (always on in debug builds) checks after each
reconcile that `compmax` and `emax` equal a full rebuild and that each carried
bound is at least the exact runner-up. Three mutation checks each made the test
suite panic:

- settling unconditionally;
- dropping the losing-candidate update;
- dropping the displaced-incumbent update.

In the accounting run the carried bound settled 99.983% of rescans, with 0
unsound settles.

## A/B

Job 33765. `dada-pooled` at 48 threads with NUMA interleave and 2 reps. The
arms are `main` (8eac43f), L1 (cut 1) and L2 (both cuts), each its own binary
with its version tag checked. Reps ran in reverse order. Every L1 and L2 output
matched `main` apart from the version string. The verify run passed on MiSeq
and PacBio.

NovaSeq, 60 soil 16S samples, 2.31 M uniques, 27,464 clusters (seconds,
rep 1 / rep 2):

| arm | wall | loop | compare (control) | shuffle | collect | rescan |
|---|---|---|---|---|---|---|
| main | 1338 / 1332 | 1265 / 1258 | 791 / 785 | 379 / 379 | 189 / 189 | 161 / 160 |
| L1 | 1228 / 1261 | 1153 / 1187 | 770 / 792 | 290 / 302 | 100 / 102 | 161 / 171 |
| L2 | 1075 / 1090 | 1002 / 1017 | 765 / 782 | 145 / 144 | 116 / 115 | 0.22 / 0.22 |

Compare, which neither cut touches, moves by ±13 s across arms and reps. Under
L2, collect rises by ~15 s over L1 because A1 now reprices and settles inside
it.

Controls (mean of 2 reps):

| pool | shuffle main → L2 | wall main → L2 |
|---|---|---|
| MiSeq 362 | 5.8 → 1.9 s | 55.9 → 50.5 s (−10%) |
| PacBio 95 | 9.2 → 3.1 s | 200.9 → 192.9 s (−4%) |

Peak RSS is unchanged on every pool.

## What this dictates

- **Shuffle is no longer the NovaSeq target.** Shuffle is 144 s against 765 s
  of compare and 87 s of `p_update`. Most of what remains is collect, 115 s.
  A1 still walks every fallen cluster's comp vec. A possible follow-up walks
  the fallen cluster's members instead, since after a move pass a non-centre
  raw's incumbent is the cluster that holds it. That is unpriced.
- **Account before cutting, and in the units of the work.** The issue's single
  ns/comp figure mixed two costs that differ by ~10× per unit. Split by pass,
  each cut was priced before it was built and landed within a few percent of
  its price: 87 s predicted for cut 1, 83 s measured.
- **An exact bound is affordable when it rides an existing sequential walk.**
  [#124][124]'s bounded build scan pruned 56% of comparisons and still lost,
  because the survivors ran at the scattered rate. This bound costs one compare
  inside A1's cluster-major walk, and each settle saves ~65 ns of scattered
  rescan work.

## Reproducing

```sh
dada2-rs dada-pooled ... --metrics-json m.json   # runs[0].reconcile: collect, rescan, affected, settled
DADA2RS_A2_WALK_FALLEN=1    dada2-rs dada-pooled ...  # cut 1 off
DADA2RS_RECONCILE_NO_BOUND=1 dada2-rs dada-pooled ... # cut 2 off
DADA2RS_RECONCILE_VERIFY=1  dada2-rs dada-pooled ...  # exactness check (slow)
```

Turning either cut off leaves the output unchanged.

[124]: https://github.com/HPCBio/dada2-rs/issues/124
[139]: https://github.com/HPCBio/dada2-rs/issues/139
[358]: https://github.com/HPCBio/dada2-rs/issues/358
