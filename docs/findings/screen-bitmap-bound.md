# An exact bound in front of the k-mer screen

**Verdict: on by default for every dense screen (k ≤ 7).** A per-raw presence
bitmap rejects most far pairs without reading their k-mer count vectors, and
never rejects a pair the full screen passes. Outputs are byte-identical.
`dada-pooled` wall drops 9% on MiSeq, 11% on PacBio and 30% on NovaSeq;
NovaSeq `learn-errors` drops 46%. Issue
[#295](https://github.com/HPCBio/dada2-rs/issues/295).

## Starting point

After batched alignment ([#337](batched-alignment.md)), the k-mer screen was the
largest part of the compare step everywhere. On the PacBio 95 pool at 48
threads its pass took ~80 s against ~57 s for alignment; on the NovaSeq pool,
~1,380 s against ~100 s. The screen reads a `4^k`-byte count vector per pair
(1 KB at k=5, 16 KB at k=7), and ~89% of PacBio pairs and ~94% of NovaSeq pairs
fail it.

## The bound

Per raw, store a `4^k`-bit presence bitmap and
`excess = nkmers − distinct`. For two raws with `D` shared distinct k-mers:

```
Σ min(a_i, b_i)  ≤  D + min(excess_a, excess_b)
```

because over the shared k-mers `min(a, b) ≤ 1 + (a − 1)`, and `Σ (a − 1)` is
at most `excess`. So `1 − (D + min excess) / scale` is a lower bound on the
k-mer distance; a pair whose bound already exceeds `--kdist-cutoff` fails the
screen unread. `D` is one AND-popcount over the bitmaps, 8× fewer bytes than
the counts (with an AVX2+POPCNT build, since baseline x86-64 has no POPCNT).
`excess` is counted from the sequence, not from the 8-bit counts, which
saturate at 255 and would make the bound unsound on the overflow path.

The tests guard both: `bound_is_exact_and_bites` fails if `excess` is dropped
or taken from the saturating counts, and `kmer_bound_matches_off` runs
`dada-pooled` at k=5 and k=7 with the bound on and off and fails on an unsound
bound.

## The gate, and why it misled at k=5

An offline check (`examples/kmer_bound_check.rs`, 100 centres against every
unique, 48 threads) measured tightness — the share of screen rejections the
bound also makes — and timed the screen with and without the bound:

| Pool | k | Tightness | Screen time, bound / full |
|---|---|---|---|
| PacBio 95 | 7 | 96.2% | 0.43 |
| MiSeq 362 | 5 | 83.2% | 0.69 |
| NovaSeq | 5 | 88.5% | 0.90 |
| PacBio 95 | 6 | 43.2% | 1.19 |

The bound was exact everywhere (zero unsound pairs; dropping `excess` produced
9,481 on MiSeq SOP). The issue's planned gate, ≥99% tightness, was missed
everywhere, and the timed arms suggested building k=7 only. The end-to-end A/B
overturned that for k=5: NovaSeq's screen pass halved where the gate predicted
10%. The 100-centre loop is not a faithful model of the production screen
there, so the decision was taken from end-to-end runs, k=6 included.

## Results

One binary, bound on (default) vs off (`DADA2RS_KMER_BOUND=0`), 48 threads,
NUMA interleave, `hpcbioamd`. Two reps with the arm order swapped (k=6: one
rep). Outputs identical in every pair of runs.

| | MiSeq k=5 | PacBio k=7 | NovaSeq k=5 | PacBio k=6 |
|---|---|---|---|---|
| `dada-pooled` wall | 64.1 → 58.2 s (−9%) | 240.3 → 214.7 s (−11%) | 2,453 → 1,707 s (−30%) | 229.8 → 227.8 s (−1%) |
| screen pass | 20.8 → 14.7 s | 81.3 → 41.9 s | 1,379 → 610 s | 40.4 → 35.3 s |
| `learn-errors` wall | 9.2 → 9.0 s | 12.7 → 10.7 s (−16%) | 109 → 58.5 s (−46%) | 10.1 → 10.0 s |
| peak RSS | +2% | +1.0 GB (+5%) | +0.4 GB (+3%) | +0.1 GB (+1%) |

Jobs 33554 (gate), 33556 and 33603 (A/B); data in
`results/2026-11-25-issue295/` on the cluster.

## Consequences

- **On by default for k ≤ 7.** `DADA2RS_KMER_BOUND=0` turns it off for A/Bs; it
  is reported as the `kmer bound` tuning gate. The sparse screen (k ≥ 8) does
  not use it.
- **It is skipped under `--screen-audit`**, which records each backend's
  verdict from real distances.
- **Gate on end-to-end runs, not on the offline loop.** The loop got k=7 right
  and k=5 wrong by a factor of five.
- **The screen is no longer the largest cost on PacBio** (42 s screen against
  58 s align), so overlapping the two passes (#342 item 3) is worth less than
  it was. On NovaSeq the screen (610 s) still dominates.
