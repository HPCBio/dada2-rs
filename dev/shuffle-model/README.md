# `b_shuffle` build-scan cost model

Standalone harness for [#124][124]. It reproduces the access pattern of
`b_shuffle_converge`'s **build** phase — the cluster-major walk of every
cluster's comp vec that rebuilds `compmax` once per bud round — in isolation,
and A/Bs candidate redesigns against it.

It is a *model*, not a benchmark of the real code: the point is to gate a design
before implementing it, per the discipline in
`docs/findings/shuffle-compmax-carry.md`. Every variant computes the same
argmax, and the printed checksum (Σ of each raw's best cluster index) must be
identical across all variants — that is the harness's exactness check.

## Build and run

```sh
rustc -O -C target-cpu=native build_model.rs -o build_model
./build_model [nraw] [nclusters] [density] [reps] [decay] [zipf]
```

| Arg | Default | Meaning |
| --- | --- | --- |
| `nraw` | 285000 | raws in the pool (pooled MiSeqSOP is ~273k–297k) |
| `nclusters` | 400 | clusters at the point being modelled |
| `density` | 0.082 | fraction of raws with a stored comp per cluster (cluster 0 always holds all of them) |
| `reps` | 5 | repetitions; the best time is reported |
| `decay` | 30 | lambda spread: `lambda = exp(-decay × U)`. Lower = flatter, harder to prune |
| `zipf` | 0 | cluster-abundance skew. 0 = uniform reads; >0 = `reads ∝ rank^-zipf` |

`decay` and `zipf` are the two parameters the result is actually sensitive to,
so sweep them rather than trusting a single cell. `nraw × nclusters × density`
sets the comparison volume; calibrate it against the `comps/build` figure in the
`--verbose` `shuffle scan split` line of a real run.

## Variants

| | Design |
| --- | --- |
| **A** | current: AoS 24 B `Comparison` stream, `emax: Vec<f64>` + `compmax: Vec<Comparison>` |
| **B** | fused 16 B per-raw map `(f64, u32)`, AoS comps |
| **C** | fused 16 B map + SoA comp streams (lambda `f64` + index `u32`) |
| **D** | split `emax` `f64` + `best_ci` `u32`, SoA comps |
| **E** | control: sequential stream, no gather (lower bound on the current shape) |
| **F** | raw-major CSR candidate index, lambda-descending, `lambda × max_reads` prune |
| **G** | two-level bound: top-`T` clusters by reads scanned cluster-major, the rest bounded by `reads_T` (run at T = 8, 32, 128) |

See `docs/findings/shuffle-build-scan.md` for the results and what they dictate.

[124]: https://github.com/HPCBio/dada2-rs/issues/124
