# Inside `b_compare`: the store is a cache-line problem

**Verdict: confirmed and built.** The serial store in `b_compare` was not doing
serial work — it was walking a 160-byte struct to read 8 bytes of it. Moving
that one field into a dense array cuts the store by **71%** and `run_dada` by
**20.1–21.3%** on NovaSeq soil 16S, byte-identical.

## What was left to explain

[The serial remainder](compare-serial-fold.md) folded six reduction passes into
the store loop for −24 to −30% of `run_dada`, and ended by naming what remained:
the store itself, at a "suspiciously flat 18.4–19.6 ns/raw-visit across a 5×
range in pool size."

That flatness was the clue. A phase whose cost per unit does not move when the
working set moves five-fold is not bound by residency. It is paying a fixed
per-access price — which is what a cache miss on every access looks like.

## The commit path is negligible; it is all scan

Before measuring anything on a cluster, the arithmetic settles what the store
*is*. On soil 16S R1 the store walks 11,888,798,623 raw-visits and commits
12,082,545 of them: **0.102%**. Whatever the store costs, 99.9% of it is the
scan over the result vector and the per-raw state it touches, not the work done
on a hit. So the target is the scan, and the commit path can be ignored.

This is the cost-model-first discipline the project keeps re-learning: the
question "what is this phase made of" is answerable on paper here, and answering
it first avoided building against the wrong 0.1%.

## The microbenchmark: `examples/store_bandwidth`

Seven arms decompose the scan into the streams it actually performs, at the real
pool size (1,225,523 raws) and the real store rate (0.00102). Reads go through
`black_box` — the first draft reported 574 GB/s because the single-array arms
auto-vectorised while the paired arms could not, an artifact of the instrument,
not a result.

| arm | what it isolates |
|---|---|
| `both48` / `both16` | result vector + strided `Raw` read, at 48 B and 16 B elements |
| `dense48` / `dense16` | result vector + a *dense* `f64` array |
| `comps48` / `comps16` | the result vector alone |
| `stride` | the strided `Raw` read alone |

The answer was unambiguous and reproduced on two machines and under both the
`release` and `release-native` profiles: **the strided read is 83–87% of the
scan.** The result vector, which is a clean sequential stream, is nearly free.

`Raw` is 160 bytes. `e_minmax` is an `f64`. Reading it once per raw pulls a
64-byte line to use 8 bytes, and the raws are visited in an order that defeats
the prefetcher. Held instead in a `Vec<f64>` on `B` — the same parallel-array
pattern `B::raw_cluster` already uses — the identical reads become contiguous.

## Result

Two replicates per arm, NovaSeq soil 16S and soil ITS2, both reads, exclusive
node (compute-5-6, EPYC 7713), 64 threads, `numactl --interleave=all`.

| soil 16S | `main` | dense array | change |
|---|---|---|---|
| store R1 | 227.2 / 223.7 s | 63.8 / 63.6 s | **−71.7%** |
| store R2 | 275.4 / 277.9 s | 80.2 / 81.2 s | **−70.8%** |
| store ns/raw-visit R1 | 19.1 / 18.8 | 5.4 / 5.3 | 18.9 → 5.4 |
| store ns/raw-visit R2 | 18.7 / 18.9 | 5.5 / 5.5 | 18.8 → 5.5 |
| `run_dada` R1 | 834.2 / 820.5 s | 654.1 / 667.4 s | **−20.1%** |
| `run_dada` R2 | 926.4 / 960.5 s | 735.4 / 749.1 s | **−21.3%** |

Soil ITS2, the smaller and more variable pool, returns −69% of store time and
**−15.2%** of `run_dada` on both reads. All 300 per-sample outputs across the
six 16S runs are byte-identical once the embedded build hash is excluded, and
240/240 on ITS2.

The controls are flat. `b_shuffle` sits at 112.5–115.3 s across every 16S arm
and `reconcile` at 40.2–41.2 ns/comp — the untouched phases did not move, which
is what rules out a whole-run effect masquerading as a phase win.

The observable the user noticed first was core utilisation. Effective cores —
`(map thread-busy + serial wall) / run_dada` — rise from **31.4 to 38.4** on R1
and **29.4 to 36.5** on R2, of 64. Amdahl is the whole story of this phase, and
the store was the largest remaining serial block.

## What generalises

The synthetic under-predicted. `store_bandwidth` put the strided read at 83–87%
of the scan, which brackets a −6 to −17% `run_dada` win; the real deficit
multiplied out to −20%. Reported here because the direction of the miss matters:
the microbenchmark isolated the right *mechanism* while understating its
*magnitude*, so its ranking of levers was trustworthy even though its arithmetic
was not. Plan against a synthetic's ordering, not its numbers.

The reusable shape is the diagnostic, not the fix. **A per-unit cost that does
not move with working-set size is a cache-line problem, not a residency
problem** — and the fix is a layout change, not a smaller working set. The same
signature is worth looking for anywhere the profile shows a flat ns/unit across
pools that differ several-fold. `Raw` is 160 bytes and other phases read single
fields from it.

## What this leaves

The store is now 12–15% of compare, down from 33–37%, and 4.9–5.5 ns/raw-visit
is close to the sequential-stream floor the synthetic measured. There is no
large serial block left in `b_compare`: the parallel map is 66% of the phase and
bandwidth-bound above 48 threads, which is
[a different problem](compare-screen-vs-align.md) with a different lever.

Two smaller things are open and measured but not decided:

- **Collecting the map's result at 16 B instead of 48 B** in production runs
  (`--verbose` needs the 48 B form for its cost attribution). This is real but
  ambiguous: −5.9% / −1.5% wall on ITS2, `sys` time collapsing 93%, and `user`
  time *rising* 6–12% with no mechanism established for the rise. The leading
  account is glibc's dynamic `mmap` threshold, which caps at 32 MB: the 48 B
  vector is 39.6–62.6 MB and is `munmap`ped and re-faulted every call, while the
  16 B vector fits under the cap and is reused. That predicts the *allocation*,
  not the element width, is the thing to fix.
- **Reusing one scratch buffer** across calls rather than allocating per call,
  which would remove the allocation at either width.
