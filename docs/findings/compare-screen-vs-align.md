# Inside `b_compare`: the screen, the DP kernel, and where the time actually goes

*Issue [#127](https://github.com/HPCBio/dada2-rs/issues/127). Measurement only —
no algorithm or output changed. 362-sample MiSeq (R1 + R2) and 95-sample PacBio
HiFi, pooled `dada`, 24 threads, both platforms pinned to the same node.*

## Why this measurement

Closing [#124](shuffle-build-scan.md) left `b_compare` as the only phase big
enough to matter: **77% / 63% / 88% of `run_dada`** on MiSeq R1, MiSeq R2 and
PacBio respectively. Its parallel map runs at **95–99% efficiency on 24
threads**, so there is no threading headroom — any win has to be algorithmic.

We had comparison *counts* but no *time* split, and the counts alone pointed two
incompatible ways. The k-mer screen is paid on **100%** of comparisons; the
alignment only on the 12–25% that pass it. So a screen that is individually
cheap can still dominate by volume, and the ratio decides the target:

- **screen-dominated** → the lever is the screen itself: an inverted index or
  minimizers, replacing the ESPRIT-era linear scan (never the partition — see
  [k-mer screen provenance](kmer-size-screening.md)).
- **align-dominated** → the lever is the aligner, which is where most prior
  effort already went (#51 WFA, band size, #15 k-mer size).

`--verbose` now reports the split directly, gated on the verbose flag because an
`Instant` pair costs ~40–50 ns and the screen itself is sub-µs.

## The result

Percentages are of `b_compare`'s **busy** time (summed across worker threads);
`ns/comp` divides by the number of units each row is actually paid on.

| | **MiSeq R1** | **MiSeq R2** | **PacBio HiFi** |
|---|---|---|---|
| uniques (`nraw`) | 272,574 | 296,893 | 547,273 |
| clusters | 2,564 | 2,007 | 2,844 |
| comparisons | 482.6 M | 400.9 M | 771.7 M |
| passed the screen | 25.0% | 25.5% | **11.9%** |
| `b_compare` busy | 2,427 s | 1,427 s | 13,426 s |
| map / parallel eff. | 105 s / 97% | 62 s / 95% | 564 s / 99% |
| **k-mer screen** | 405.8 s (**16.7%**) | 327.7 s (**23.0%**) | 2,038.5 s (**15.2%**) |
| — per comparison | 841 ns | 818 ns | 2,642 ns |
| **align total** | 1,877.8 s (**77.4%**) | 994.3 s (**69.7%**) | 10,958.5 s (**81.6%**) |
| — DP kernel | 1,741.1 s (71.7%) | 909.3 s (63.7%) | 10,454.5 s (77.9%) |
| — `al2subs` + quals | 136.8 s (5.6%) | 85.0 s (6.0%) | 504.0 s (3.8%) |
| — per aligned pair | 15,537 ns | 9,713 ns | 119,634 ns |
| other (`compute_lambda`, overhead) | 143.7 s (5.9%) | 105.1 s (7.4%) | 429.1 s (3.2%) |

**Both platforms are align-dominated, and PacBio is the *most* align-dominated
of the three.** That refutes the prediction recorded when #127 was opened: k=7
gives PacBio 16,384-entry k-mer vectors against MiSeq's 1,024, and 88% of its
pairs are shrouded, so PacBio looked like the screen-bound case. The screen does
cost 3.1× more per comparison there (2,642 vs 841 ns) — but the alignment it
guards costs **7.7×** more (119.6 vs 15.5 µs), and it fires on 2.1× fewer pairs.
The two effects move in opposite directions and the alignment wins.

So there is one lever, not two. Rolling the DP kernel up to whole-run shares:

| | MiSeq R1 | MiSeq R2 | PacBio |
|---|---|---|---|
| DP kernel, % of `run_dada` | **43.7%** | **29.7%** | **62.8%** |
| DP kernel, % of pooled wall | ~41% | ~28% | ~58% |

**The banded Needleman-Wunsch DP kernel is the single hottest thing in pooled
`dada` — by a wide margin, on both platforms.** Every phase #124 examined is a
rounding error beside it.

## The screen is not the problem, and it is worth keeping

Worth stating plainly, because "16.7% of compare" invites the wrong conclusion.
Without the screen, every comparison would align: MiSeq R1's 482.6 M
comparisons × 15.5 µs ≈ **7,450 s** against the 2,427 s actually spent. On
PacBio it is 771.7 M × 119.6 µs ≈ **25.6 hours** of busy time against 3.7. The
screen returns roughly **3× on MiSeq and 7× on PacBio**.

A *perfect, free* index — one that costs nothing and shrouds exactly what the
current screen shrouds — is therefore capped at 15–23% of compare. Worth having
eventually; not the thing to build first, and not the thing #127 was opened to
find.

## Where the DP kernel's time goes

The per-aligned-pair costs are large enough to be worth decomposing. The
vectorized aligner walks anti-diagonals of a banded matrix, so cells computed
≈ `nrow × (band+1)` with `nrow = len1 + len2 + 1`:

| | cells/align | ns/align | **ns/cell** |
|---|---|---|---|
| MiSeq R1 (240 bp, band 16) | 8,177 | 14,405 | **1.76** |
| MiSeq R2 (160 bp, band 16) | 5,457 | 8,882 | **1.63** |
| PacBio (~1,455 bp, band 32) | 96,063 | 114,132 | **1.19** |

Two things are reassuring here: cost scales with cells (MiSeq R2/R1 = 1.62×
measured vs 1.50× from the length ratio), and the vectorized path is confirmed
in use — HiFi reads sit under the 3,500 bp `i16` guard and take
`homo_gap_p == gap_p`, so neither fallback fires.

But **1.2–1.8 ns/cell is 5–10× off what an 8-lane `i16` SIMD kernel should
achieve**, and the two platforms appear to miss for different reasons.

### The memory traffic is fixed at 6 bytes per cell

`align_vectorized_with_buf` allocates the **full** compressed matrix — `d16`
(`i16` scores) *and* `p8` (traceback pointers), both `ncol × nrow` — and
`reset_buf` zero-fills both on every call. Each cell is therefore touched twice
(zeroed, then written) across 3 bytes: **6 bytes of write traffic per cell**,
independent of platform.

| | matrix (`d16` + `p8`) | traffic/align | GB/s/thread | × 24 threads |
|---|---|---|---|---|
| MiSeq R1 | 18 KB + 9 KB | 54 KB | 3.8 | 91 GB/s |
| MiSeq R2 | 12 KB + 6 KB | 36 KB | 4.1 | 99 GB/s |
| PacBio | 199 KB + 99 KB | 597 KB | 5.4 | **129 GB/s** |

All three land in the same 91–129 GB/s band despite a **12× spread in matrix
size**. On PacBio that is a bandwidth-bound signature: a 298 KB working set per
alignment, streamed, 24 threads deep. On MiSeq the 27 KB matrix is L2-resident,
so DRAM bandwidth cannot be binding there — and yet MiSeq's ns/cell is *worse*.
The likely mechanism there is per-anti-diagonal overhead: at band 16 the inner
loop runs ~17 cells, or two `i16` vectors plus a remainder, so loop
prologue/epilogue is amortised over very little work. PacBio's band-32
diagonals (~33 cells, four vectors) amortise it better, which is consistent with
its lower ns/cell despite the larger matrix.

Different mechanisms, same target, and one structural observation applies to
both: **`d16` is write-only outside a 3-row window.** The DP reads rows
`row`, `row-1` and `row-2`; the traceback reads `p8` only and never touches
`d16`. The full `d16` matrix — 199 KB per alignment on PacBio — is allocated,
zeroed and written for nothing. The zero-fill is redundant in the same way: the
sentinel columns are written explicitly and the active cells are all overwritten
before they are read.

That suggests a cheap experiment before any redesign: keep `p8` full (the
traceback genuinely needs it), reduce `d16` to a 3-row rolling buffer, and drop
both zero-fills. If the bandwidth reading is right, per-cell traffic falls from
6 bytes to ~1 and PacBio should move a lot. It is a pure-mechanics change with
no effect on the DP recurrence, so the byte-identical invariant should be
straightforward to hold — but it must still be *held*, not assumed: #124's
pruning arm hid a real tie-break bug that only exactness testing surfaced.

## Design constants

For costing any future change without re-running (24 threads, pinned node;
`ns/comp` and `ns/align` are per-thread busy time):

| constant | MiSeq (k=5, band 16) | PacBio (k=7, band 32) |
|---|---|---|
| k-mer screen | 818–841 ns/comp | 2,642 ns/comp |
| screen pass rate | 25.0–25.5% | 11.9% |
| DP kernel | 1.63–1.76 ns/cell | 1.19 ns/cell |
| `al2subs` + qual mapping | 830–1,132 ns/align | 5,502 ns/align |
| `compute_lambda` + overhead | 5.9–7.4% of compare | 3.2% of compare |
| DP write traffic | 6 bytes/cell | 6 bytes/cell |
| comparisons per unique | 1,350–1,771 | 1,410 |

Note the last row: `b_compare` is the `O(nraw × nclusters)` term, and at ~1,400
comparisons per unique on both platforms it is the reason pooled runs scale the
way they do. Nothing in this page changes that shape — it identifies which
constant in front of it is worth attacking.

## What this dictates

1. **`b_compare` is align-bound on both platforms.** One lever, not the two
   platform-specific ones #87 ended with.
2. **The k-mer screen is closed as a perf target.** It returns 3–7× on what it
   costs, and a perfect free replacement caps out at 15–23% of compare. Revisit
   only if the aligner gets fast enough to change the ratio.
3. **The DP kernel is the target** — 30–63% of `run_dada`, the largest single
   share this project has measured.
4. **The first thing to try is not a new algorithm.** The kernel moves 6 bytes
   of memory traffic per DP cell, ~5/6 of it for a matrix that is never read
   back. Test the rolling-buffer change before considering anything that
   perturbs the alignment itself; per
   [cost-model-before-build](shuffle-build-scan.md), model it against the
   ns/cell constants above and not against a laptop.
5. **`al2subs` is not the problem** (3.8–6.0%), which retires the hypothesis
   that post-processing explained the per-pair cost.
