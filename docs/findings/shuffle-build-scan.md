# Pruning the `b_shuffle` build scan

**Verdict:** the build phase of `b_shuffle_converge` is **not layout-bound — it
is volume-bound**, and the volume is removable exactly. Data-layout rework
(narrower records, SoA comp streams, a fused per-raw map) is worth **0–9%** and
is usually a *regression*. Replacing the cluster-major full scan with a
**two-level bounded scan** examines only **9–24% of the candidates** and runs
**1.06–4.24× faster** across a 24-cell parameter sweep, never slower, with a
**bit-identical argmax** in every cell. This is the first direction from
[#124][124] to clear the gate.

The result is a model, not a measurement of the shipped code: it gates an
implementation, it does not substitute for one. See "What this does not show".

## Background

After [#86][86], `b_shuffle` is 46–55% of pooled wall and the Amdahl ceiling on
pooled runs. It splits into a **build** (cluster-major rebuild of `compmax`,
paid fresh each bud round, 27–41% of shuffle *time*) and a **reconcile**
(raw-major, scattered, the remainder). [#87][87] tried to remove the build by
carrying `compmax` across buds and was closed: it merely relocated the work into
the more expensive scattered path, and the sign of the result flipped between
two amplicons of one sequencing run (`docs/findings/shuffle-compmax-carry.md`).

That closure left the phase itself as the target. The structural fact #124
starts from is that during a shuffle **`lambda_ij` is fixed** — only the scalar
`reads_j` changes — so the whole phase is `argmax_j (lambda_ij × reads_j)` with
one factor static, evaluated ~9.6e6 times per build to discover ~550 moves.

## Method

`dev/shuffle-model/build_model.rs` reproduces the build loop's access pattern in
isolation and A/Bs seven designs on identical generated data. Every variant
prints a checksum of the resulting best-cluster assignment; **all variants agree
in all cells**, so the speed comparison is between designs that produce the same
partition, including the lowest-`ci` tie-break.

Two workload parameters turn out to carry the whole result:

- **`decay`** — the spread of `lambda`. Steeper decay (lambda falling fast with
  Hamming distance) makes a raw's top candidate more dominant.
- **`zipf`** — the skew of cluster read counts, `reads ∝ rank^-zipf`.

Measurements below are on Apple silicon. The *ratios* and the *examined
fraction* transfer; the absolute ns/comp does not (see below).

## Result 1: layout rework is a dead end

At the default cell (nraw 285k, 400 clusters, 9.6e6 comps):

| Variant | ns/comp | vs current |
| --- | ---: | ---: |
| A current (AoS 24 B comps, `emax` + `compmax` 24 B) | 1.65 | 1.00× |
| B fused 16 B per-raw map | 1.84 | 0.89× |
| C fused map + SoA 12 B comp streams | 1.78 | 0.93× |
| D `emax` + `best_ci` u32, SoA comps | 1.54 | **1.07×** |
| E control: sequential stream, no gather at all | 0.91 | 1.81× |

The current build is already within **1.8×** of a loop that does the same
arithmetic with *no gather whatsoever*. That is the ceiling on all layout work,
and the best real variant captures 7% of it. Shrinking the per-raw map from
32 B to 16 B makes things **worse**, not better. Repeating the sweep at
nraw = 1e6 and 3e6 does not change the ranking (D peaks at 1.09×).

So: the per-raw map spilling out of cache is *not* what the build costs. Its
cost is one multiply-and-compare per stored comparison, ~9.6e6 of them, and the
only way to make it cheaper is to **not do them**.

## Result 2: the naive bound works, then collapses under real abundance skew

Direction 1 in #124 — sort each raw's candidate list by lambda descending
(static, since lambda is fixed) and break once `lambda × max_reads ≤ best` — is
exact and, iterated raw-major over a **flat CSR** index rather than the current
per-raw `Vec`, is also a *sequential* stream.

With uniform cluster reads it is excellent: **5.3% of candidates examined,
3.05×**. Under a power-law abundance distribution it falls apart:

| decay | zipf | examined | vs current |
| ---: | ---: | ---: | ---: |
| 30 | 0 | 5.3% | 3.05× |
| 10 | 1.0 | 33.9% | 1.31× |
| 4 | 1.0 | 48.1% | **0.83×** |
| 10 | 2.0 | 45.9% | **0.87×** |

The bound is `lambda × max_reads`, so **a single very large cluster inflates the
bound for every raw in the pool** — and real cluster abundances are exactly that
skewed. Taken alone this direction would have been another workload-conditional
result with a sign flip, like #87.

## Result 3: a two-level bound is robust

The large clusters are *few*. Split the build accordingly:

1. Scan the **top `T` clusters by reads cluster-major** — contiguous, the access
   pattern the build already uses — to seed `emax`/`best_ci`.
2. Scan everything else **raw-major over the CSR index**, bounded by `reads_T`,
   the (T+1)-th largest cluster read count, which under skew is orders of
   magnitude below `max_reads`.

Selecting the top `T` costs O(nclusters) per build — hundreds of elements,
negligible. The composition is order-independent because the max uses
`(e > best) || (e == best && ci < best_ci)`, so both passes yield the same
lowest-`ci` argmax regardless of visit order.

`T = 32` over a 24-cell sweep (nraw ∈ {285k, 1e6} × density ∈ {0.03, 0.082,
0.2} × four `decay`/`zipf` combinations):

- candidates examined: **8.9%–23.5%** of the full scan
- speed: **1.06×–4.24×**, **never below 1.0×**
- checksum identical to the current build in **every cell**

The weak cells are the low-skew, flat-lambda corners (1.06–1.27×) — precisely
where the naive bound was strongest and where there is least to win. The strong
cells are the skewed ones, which is where real data lives. `T` is not delicate:
T = 8 and T = 32 are close everywhere, T = 128 is consistently worse (pass 1
starts scanning clusters that pass 2 would have pruned).

## What this dictates

- **Do not spend effort on the build's data layout.** Result 1 closes that:
  the ceiling is 1.8× and the achievable share of it is ~7%, with two of the
  three obvious "compaction" designs coming out negative. This is a fifth
  negative result for reasoning about this phase in bytes rather than volume.
- **Implement the two-level bounded build** as the #124 candidate, with the
  lambda-descending CSR candidate index. It needs the CSR index built once and
  extended per bud (the existing `CandIndex` is already appended per cluster;
  the change is flattening it and sorting each raw's slice by lambda descending,
  which is legal because lambda never changes after `b_compare`).
- **Gate it on the ASV concordance guardrail, not on wall time alone.** The
  harness checksum shows the argmax is exact by construction, but the real
  implementation must preserve `raws[i].comp` and the tie-break end to end.
- **Expect the real gain to exceed this model, not fall short of it.** The
  target machine's build runs at 4.3–6.3 ns/comp against this machine's 1.65,
  i.e. it is substantially more memory-limited. A design that cuts examined
  volume 4–10× should pay *more* there, not less. That is an argument for
  measuring on the target before sizing, not for assuming the upside.
- **The same machinery should then be pointed at the reconcile**, which is the
  larger half (59–73% of shuffle time, 11.4–13.4 ns/comp). It is already
  raw-major, so it inherits both the CSR flattening (removing a per-raw `Vec`
  pointer chase) and the bound. This is the follow-up, not part of the first
  change.

## What this does not show

- It is synthetic. `decay` and `zipf` are stand-ins for the real lambda spectrum
  and abundance distribution; the sweep is wide enough that no cell is a loss,
  but the *specific* speedup on 16S or ITS2 is not predicted here.
- It models the **build only**. Total shuffle time also contains the reconcile,
  so a 2× build is roughly a 13–20% shuffle reduction and ~6–11% of pooled wall
  on its own — worthwhile but not transformative until the reconcile follows.
- It ignores the cost of maintaining the sorted CSR index across buds. That is
  real work an implementation must pay and this harness does not charge for.

## Reproducing

```sh
cd dev/shuffle-model
rustc -O -C target-cpu=native build_model.rs -o build_model
./build_model 285000 400 0.082 3 10 1.0     # nraw nclusters density reps decay zipf
```

Calibrate `nraw`/`nclusters`/`density` against a real run's `--verbose`
`shuffle scan split` line (`comps/build`) before trusting a cell.

[86]: https://github.com/HPCBio/dada2-rs/pull/86
[87]: https://github.com/HPCBio/dada2-rs/issues/87
[124]: https://github.com/HPCBio/dada2-rs/issues/124
