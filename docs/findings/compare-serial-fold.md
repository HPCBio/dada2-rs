# `b_compare`'s serial remainder: attributing 329 s, then deleting it

*Issue [#143](https://github.com/HPCBio/dada2-rs/issues/143). Two NovaSeq soil
pools (16S and ITS2, R1 + R2), pooled `dada`, 64 threads, exclusive node,
`numactl --interleave=all`. Byte-identical output — no algorithm or result
changed.*

## Why this measurement

[#127](compare-screen-vs-align.md) established that `b_compare`'s parallel map
is execution-bound below ~48 threads and bandwidth-bound above, so the map had
no headroom left. What it did not explain was the rest of `b_compare`. On the
soil 16S pool the phase cost 946 s against a 426 s map and a 182 s store,
leaving **329 s — 48% of `run_dada` — unaccounted for**.

The issue was explicit that this was an attribution task, not an optimisation
one:

> Only then decide whether anything is worth building. Do not design an
> optimisation against a share-of-wall figure before re-measuring what it is
> made of.

That order is what made the fix obvious, and it is the transferable lesson here.
A 329 s hole invites theories; measuring it first showed there was nothing to
theorise about.

## The split

`--verbose` now reports `b_compare` as map / reduction / store / free / setup,
with an explicit `unattributed` remainder so the split cannot silently stop
adding up.

| soil 16S R1 | | |
|---|---|---|
| map | 446.3 s | parallel |
| reduction | 337.0 s | serial, 28.3 ns/raw-visit |
| store | 181.5 s | serial, 15.3 ns/raw-visit |
| free | 4.9 s | serial |
| setup | 0.0 s | serial |
| **unattributed** | **0.01 s** | |

The remainder closes to 0.01 s on both pools. The missing time was a single
thing: **six serial passes over the map's result vector**, summing per-item
costs and denominators, plus the vector's free.

## The fix

Those six passes walk the same vector the store loop already walks. Folding them
in makes it one pass instead of seven, removing ~288 B of memory traffic per
raw-visit.

| | pre-fold | post-fold | |
|---|---|---|---|
| 16S R1 `run_dada` | 1101.3 s | 781.5 / 751.6 s | **−30.4%** |
| 16S R2 `run_dada` | 1262.1 s | 891.2 / 879.7 s | **−29.8%** |
| ITS2 R1 `run_dada` | 312.1 / 316.3 s | 235.2 / 234.7 s | **−25.2%** |
| ITS2 R2 `run_dada` | 407.0 / 426.5 s | 315.2 / 308.6 s | **−25.2%** |
| 16S R1 serial in compare | 518.5 s | 217.6 s | **−58%** |
| 16S R2 serial in compare | 622.2 s | 271.1 s | **−56%** |

The soil 16S baselines are single-rep; reading the most conservative post-fold
value against them instead gives −24.0% / −27.6%. The honest range is **−24 to
−30%**, and every arm lands inside it.

120 per-sample outputs byte-identical against the pre-fold baselines. The
control channel (`shuffle`, untouched by this change) moves −1.9% to +2.9%
across the two reads — noise in both directions rather than a rig shift.

## What was falsified: cache residency

The reduction's rate first looked like a cache story. It ran at 28.3 ns/raw-visit
on soil 16S, against 15.4 ns on MiSeq, and the result vectors sit either side of
the EPYC 7713's 32 MB per-CCD L3 — 59 MB for soil, 13 MB for MiSeq. A clean
account, and wrong.

ITS2's vector is far smaller than soil 16S's, and its reduction ran at
**~29 ns/raw-visit** — indistinguishable from soil's, at a fraction of the
working set. The six passes were expensive because they were six passes, not
because of where the data lived. MiSeq's 15.4 ns still wants an explanation;
residency is no longer it.

!!! note "The store did not get free"
    The store absorbs part of the reduction's work — 15.3 → ~19 ns/raw-visit —
    so the saving is not the full 337 s. The rest was pure re-walking, which is
    the part worth remembering when a "cheap" serial pass is added to a loop
    that already exists.

## What this leaves

The store is now the entire serial block inside `b_compare`: **27–30% of
`run_dada` on soil 16S, 23–24% on ITS2**, at a strikingly consistent
**18.4–19.6 ns/raw-visit** across both pools and all four reads. That
consistency across a 5× range in pool size is the useful lead — it suggests a
per-visit cost that does not depend on the working set, which is a different
kind of target from the one just removed.

Occupancy frames the same point. Reconstructed from the phase timers and
cross-checked against SLURM accounting to within 9%, soil 16S post-fold runs at
**~30 of 64 cores**. The deficit is serial `b_compare` plus `b_shuffle`, and
after this change the store is the larger half of it.

## A measurement note that cost two runs

An A/B in this series compared a configuration against itself: the gate had been
renamed from `DADA2RS_SHUFFLE_CARRY` to `DADA2RS_SHUFFLE_NO_CARRY` when the
carry became the default, and the stale variable was accepted in silence. Both
arms ran identically, completed cleanly, and came back byte-identical and within
1% — **which is exactly what a correct null result looks like.**

Every subsequent A/B was verified against a positive control in the log itself
(`shuffle scan split` reading `1 builds` versus thousands) *before* any timing
was read. [#145](https://github.com/HPCBio/dada2-rs/issues/145) tracks reporting
resolved tuning gates so this fails loudly instead.
