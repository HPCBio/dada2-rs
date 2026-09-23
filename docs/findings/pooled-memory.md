# Pooled memory: peeling the peak, one dominant term at a time

**Pooled `dada` on 93 PacBio samples went from 52.3 GB to 13.6 GB at k=7, and
made k=8 possible at all** — 42.8 GB down to 13.6 GB, where before it simply did
not fit. Every step was byte-identical at the ASV level.

The useful thing about this arc is not any single change. It is that **each fix
moved the peak somewhere else**, so the work was a repeated measurement of "what
dominates *now*" rather than one optimisation. Four rounds, in the order the
terms surfaced:

| # | What dominated | Fix | Result |
|---|---|---|---|
| [#32](https://github.com/HPCBio/dada2-rs/issues/32) | resident u16 k-mer frequency vector | drop it; recompute on the rare fallback | k7 **52.3 → 35.7 GB** |
| [#39](https://github.com/HPCBio/dada2-rs/issues/39)/[#40](https://github.com/HPCBio/dada2-rs/pull/40) | dead intermediates held live through `dada` | lifetime management | peak drops to the derep floor |
| [#41](https://github.com/HPCBio/dada2-rs/issues/41)/[#42](https://github.com/HPCBio/dada2-rs/pull/42) | all per-sample dereps resident at once | stream derep → merge | k5 **17.0 → 9.8 GB** (−42.3%) |
| [#43](https://github.com/HPCBio/dada2-rs/issues/43)/[#45](https://github.com/HPCBio/dada2-rs/pull/45) | dense `4^k` k-mer screen array | sparse, gated to k≥8 | k8 **42.8 → 13.6 GB** (−68.3%) |

## The instrument came first

Round one shipped a `--verbose` line reporting the resident `Raw` footprint per
build — `seq+qual` and k-mer bytes, with a per-raw figure. Everything after it is
downstream of that print:

```text
[dada-pooled] peak RSS after derep+merge: 5967 MB
[dada]        resident Raw footprint: 547243 raws, seq+qual 1518.7 MB,
              k-mer vectors 10063.1 MB (19282 B/raw) [k=7]
[dada-pooled] peak RSS after dada:        17626 MB
```

That is what turns "the peak is 17.6 GB" into "10.0 GB of it is k-mer vectors,
57% of the total, and here is the gap to non-`Raw` structures." The same pattern
as the effective-cores column in [the threading work](threading-serial-steps.md):
the cheap instrument found the expensive problem.

It also corrected an earlier mistake. Two previous attempts had gone after
`Derep.map`, which phase-boundary profiling showed to be **~15 MB — negligible**.
The real cost was quals, seqs and the merge accumulator.

## Round 1 — the biggest array was a fallback

`raw_assign_kmers` populated three vectors per unique. At k=7 they cost ~51 KB
each, against ~6 KB at k=5 — a 16× blow-up on the dominant resident term, since
pooled mode holds a `Raw` for every unique across all samples simultaneously.

The `kmer` vector (u16, `4^k`, 32 KB at k=7) was the largest **and was only a
fallback**: screening uses the u8 `kmer8`, consulting u16 only when the u8
distance is borderline. Dropping it and recomputing on that rare path:

| mode | k | before → after | Δ |
|---|---|---|---:|
| pooled | 5 | 28.72 → 27.59 GB | −3.9% |
| pooled | 7 | 52.33 → 35.69 GB | **−31.8%** |
| pseudo, cached | 7 | 14.29 → 10.27 GB | −28.1% |

The win scales exactly with `4^k` (−1.1 GB at k5 against −16.6 GB at k7 ≈ 16×),
plus a k7-only **wall-clock bonus** of −5.6% from the reduced memory traffic.
Concordance 2082 = 2082, churn 0.

### Two follow-ups, both closed as non-viable

Worth recording so they are not re-proposed:

- **Free the k-mer vectors between rounds.** There is no idle window. Pooled
  `run_dada` is *all* `b_compare` — every bud round re-screens every raw, so the
  vectors are needed continuously and the peak is mid-loop. In pseudo/cached
  mode there is nothing to free: the cache holds `RawInput` without k-mers, and
  `Raw`s are built transiently per call.
- **Drop `kord`.** Load-bearing under the default `gapless = true`, read in the
  screening hot loop. Unlike the u16 vector, which is consulted rarely, `kord` is
  per-comparison — dropping it means disabling the gapless fast path, which is an
  alignment-performance trade rather than a free win.

## Rounds 2 and 3 — the merge phase was mostly dead weight

Phase-boundary profiling put the merge phase as the single biggest jump: **+9.2
GB**, against +3.9 GB for dada itself. Two causes, fixed in order.

**Dead intermediates held live.** `dereps` and the merge intermediates stayed in
scope through the `dada_uniques` call, which needs only `raw_inputs`. The
per-sample output needs each sample's unique read *counts*, not its quals and
seqs. Extracting a small `Vec<Vec<u32>>` and dropping the rest before the dada
call is pure lifetime management — and took the peak down to the derep floor.

**All dereps resident at once.** That floor was ~15 GB at 93 samples: every
per-sample derep loaded up front and held through the merge loop. Streaming it —
load one sample, fold it into the accumulator, build its small index vectors,
drop it — collapses the floor to roughly the dada working set:

**k5 dada peak: 17.04 → 9.82 GB (−42.3%)**, which is about −64% against
pre-#39 `main`.

Fold order stays strict sample index, identical to the old loop, so the merged
table and every output byte are unchanged *by construction* rather than by
measurement.

### The regression that was a misreading

A first look at k=7 showed **+4.2%** and got written down as a regression. A
paired thread sweep showed what was actually happening:

| threads | pre | post |
|--------:|-------:|-------:|
| 1 | 21.22 GB | **16.60 GB** |
| 2 | 17.89 GB | 16.75 GB |
| 4 | 16.94 GB | 16.94 GB |
| 8 | 17.40 GB | 17.21 GB |
| 16 | 17.39 GB | 17.56 GB |
| 24 | 17.01 GB | 17.88 GB |

The story is in the absolute columns, not the deltas. **`post` is a clean curve**
— 16.60 → 17.88 GB, smoothly monotonic in thread count, about 1.3 GB of
concurrent scratch across the whole range. **`pre` is erratic** — non-monotonic,
ranging 16.94 → 21.22 GB, with a 21.2 GB spike at one thread, because its peak
was the allocation-timing-sensitive all-dereps floor.

The original "+4.2%" compared `pre` near its *floor* against `post` at its *curve
top*. The tell that it was not structural: a structural offset would appear at
N=1 too, and N=1 is the single biggest **win** (−21.8%). What remains is an
allocator-under-concurrency transient worth about +5% at full saturation.

**The lesson: when a baseline is allocation-timing sensitive, compare curves, not
points.** A single paired measurement cannot distinguish a real regression from
two different points on two differently-shaped curves.

## Round 4 — the dense `4^k` array, and where sparse actually wins

With the merge collapsed, the peak was the k-mer screen itself: dense `kmer8` at
`4^7` = 16,384 B/raw, **~90% zeros** — a 1.5 kb read contains only ~1,489
distinct 7-mers. At 547k uniques that is ~8.7 GB, 57% of the 17.6 GB peak.

A sparse `(index, count)` representation makes the footprint track *sequence
length* instead of `4^k`. The distance becomes a merge-join over two sorted
vectors rather than a dense sum-of-min sweep — bit-identical, overflow sentinel
included, which unit tests assert directly.

The catch is that this is the screening hot loop, so it trades memory against
throughput on the most-travelled path in pooled dada. Sparse costs a fixed
~6 KB/raw; dense costs `4^k`. The A/B pins the crossover exactly where `4^k`
passes ~6 KB:

| k | dense `4^k` | wall Δ | RSS Δ | RSS before → after |
|---|---|---|---|---|
| 5 | 1 KB | −0.1% | −0.1% | 9.85 → 9.84 GB (dense both sides, sanity arm) |
| 6 | 4 KB | **+20.1%** | **+20.8%** | 11.65 → 14.07 GB — sparse is *bigger and slower* |
| 7 | 16 KB | +25.7% | −22.6% | 17.80 → 13.77 GB — memory win, wall loss |
| 8 | 64 KB | **−17.3%** | **−68.3%** | 42.80 → 13.58 GB — wins both, and fits at all |

So the gate is `SPARSE_KMER_MIN = 8`: sparse only where it is a pure win, and
k ≤ 7 keeps the byte-identical dense path. Churn 0 at every k on the 93-sample
pooled set.

Two things in that table are worth more than the headline. **At k=6 sparse is
larger than dense** — the fixed per-entry overhead loses to a 4 KB array, which
is exactly the kind of result an argument-from-first-principles misses. And at
k=7 the dense sweep still wins on wall despite losing on memory, because the
array fits in cache; only at k=8 does dense spill and the merge-join win both.

**The footprint is now decoupled from `4^k`** — k7 at 13.77 GB against k8 at
13.58 GB. The residual is the sequence term.

## What this dictates

- **Instrument the phase boundaries before optimising a peak.** Two earlier
  attempts went after a ~15 MB term. The per-phase RSS print is what made the
  rest of this arc addressable, and it cost almost nothing.
- **Gate on the measured crossover, not on the idea.** Sparse representation is
  obviously better in the abstract and is *worse* at k=6. The constant encodes a
  measurement, which is why it has a number in it.
- **Do not re-propose freeing k-mer vectors between rounds, or dropping `kord`.**
  Both were investigated and closed; the reasons are above.
- **The next memory term is sequence encoding**
  ([#28](https://github.com/HPCBio/dada2-rs/issues/28)), now that `4^k` no longer
  drives the footprint.
- **k=8 is feasible, which makes the `KMER_SIZE_MAX` cap an index-width artifact
  rather than a memory one** — worth raising experimentally
  ([#44](https://github.com/HPCBio/dada2-rs/issues/44)), but note that
  [k-mer size has ~no effect on the final table](kmer-size-screening.md), so any
  extra ASVs are a hypothesis to disprove rather than a win to claim.
- **The per-sample memory dial is `--sample-jobs`**, which bounds how many
  samples are in flight; see [threading](threading-serial-steps.md).
