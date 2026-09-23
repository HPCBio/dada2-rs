# The port: what was translated, what was rewritten

DADA2 is two layers, and dada2-rs treats them completely differently.

- **The C++ core was translated**, file by file, structure by structure. All
  thirteen sources, in two days. This layer is the fidelity contract: where it
  behaves differently from R, that is a bug.
- **The R layer was rewritten.** `R/*.R` is orchestration, S4 classes and
  plotting — things that do not translate into a CLI so much as get re-expressed
  as one. Here dada2-rs is free to differ, and does.

Between them sits the decision that shaped everything after it: **Rcpp was cut
out rather than reimplemented.** `Rmain.cpp` became `src/dada.rs` with no R
runtime anywhere, which is what makes this a standalone tool rather than a faster
backend for the R package.

One more thing worth stating plainly, because it explains the shape of the
project's first months: **correctness against R was not the initial goal.**
Workflow completeness was — FASTQ in, sequence table out. The concordance work
began only once that pipeline existed end to end, and the first thing it found
was the [error model](findings/loess-error-model-correctness.md).

## The C++ layer, translated

| DADA2 C++ | dada2-rs | Notes |
|---|---|---|
| `misc.cpp` | `misc.rs` | nucleotide integer encoding/decoding |
| `containers.cpp` | `containers.rs` | `Raw`, `Bi`, `B`, `Sub`, `Comparison` |
| `kmers.cpp` | `kmers.rs` | k-mer frequency/order vectors and distances |
| `pval.cpp` | `pval.rs` | Poisson abundance p-values, lambda |
| `nwalign_endsfree.cpp`, `nwalign_vectorized.cpp` | `nwalign.rs` | both, merged |
| `cluster.cpp` | `cluster.rs` | the greedy divisive loop |
| `filter.cpp` | `filter.rs` | |
| `Rmain.cpp` | `dada.rs` | **decoupled from Rcpp** |
| `taxonomy.cpp` | `taxonomy.rs` | naive Bayes k-mer classifier |
| `evaluate.cpp` | `evaluate.rs` | alignment evaluation, k-mer utilities |
| `error.cpp` | `error.rs` | cluster stats, transition counts, post-hoc p-values |
| `chimera.cpp` | `chimera.rs` | bimera detection |

Nothing in the C++ layer was left behind. `RcppExports.cpp` and `dada.h` have no
counterpart by construction.

### What changed in translation, and why it is not cosmetic

- **Manual `malloc`/`realloc` became `Vec`.** The interesting case is
  `bi_pop_raw`, which removes a raw by overwriting it with the last element —
  Rust's `swap_remove` is the same operation, so the semantics carried across
  exactly rather than being approximated by a "tidier" removal that would have
  changed ordering.
- **SSE2 intrinsics became scalar loops.** `kmer_dist_SSEi`, `kord_dist_SSEi`,
  `nwalign_endsfree_SSEi` and `dploop_vec_SSE2` were hand-written SIMD; the port
  writes plain loops and lets LLVM auto-vectorise to NEON or AVX. This was a bet,
  and it paid — see [the DP kernel work](findings/compare-screen-vs-align.md) —
  but it is also why "does our kernel match theirs cell for cell?" was an open
  question for the first months.
- **RcppParallel became Rayon**, in `b_compare_parallel` and in
  `table_bimera2`.
- **R's `ppois` became `statrs`**, which uses the same regularised incomplete
  gamma, so the p-values agree rather than merely being close.
- **Indexing is 0-based throughout**, and callers add 1 if they need R-style
  output.

## The R layer, rewritten

R's `R/` directory is not an algorithm; it is a package. `allClasses.R`,
`show-methods.R` and `plot-methods.R` are S4 and ecosystem glue.
`dada.R`, `errorModels.R`, `filter.R`, `paired.R`, `chimeras.R`, `taxonomy.R`
and `sequenceIO.R` are the workflow — the part a user actually calls.

That workflow was re-expressed as subcommands rather than translated, and the
mapping is deliberately not one-to-one. `remove-primers`, for instance, fuses
R's `removePrimers` and `filterAndTrim` into a single parallel pass — which is
[the single biggest wall-clock win in the benchmarks](findings/threading-serial-steps.md),
and is a design difference rather than a porting decision.

The two layers are visible in the source tree, and the clearest case is the word
"filter" appearing twice. `src/filter.rs` is the translated `filter.cpp` — its
header says so. `src/filter_trim.rs` is the rewritten R workflow function,
mirroring `filterAndTrim` / `fastqFilter` / `fastqPairedFilter`, and its header
documents the *order* filters are applied in to match R. Same word, different
layers, different obligations.

Plotting did not come across. `scripts/plot_errors.R`,
`plot_quality_profile.R`, `plot_complexity.R` and `plot_expected_error.R` cover
the same ground from outside, reading the JSON outputs.

## Cutting Rcpp

Removed from `Rmain.cpp` in the port:

- SSE/x86-64 SIMD dispatch
- `Rcpp::checkUserInterrupt()` — no R event loop to poll
- the R-specific `b_make_*` output formatters
- RcppParallel

What remains is two entry points that owe nothing to R:

```text
dada_uniques(&[RawInput], &DadaParams) -> Result<DadaResult, String>
run_dada(Vec<Raw>, &DadaParams) -> B
```

This is the most consequential structural choice in the project. It means
dada2-rs can be a CLI, a library and a container image without an R install
anywhere near it — and it means the algorithm's inputs and outputs are plain
Rust types rather than `SEXP`. It also means **`dada_uniques` is the boundary**:
the idea of offering dada2-rs as a backend *for* R DADA2 later on is an option
this design leaves open, not a goal it was built toward.

## Nine days to a closed pipeline

| date | what landed |
|---|---|
| Apr 7 | `summary` — quality metrics, before any algorithm |
| Apr 8 | `derep`; then the whole C++ core: foundations → `nwalign` → `cluster`/`filter` → `Rmain` → `taxonomy` → `evaluate`/`error`/`chimera` |
| Apr 9 | `error_models`, subsampling |
| Apr 13 | `learn-errors` — the iterative self-consistency loop |
| Apr 14 | `dada`, `merge-pairs` |
| Apr 15 | `filter-and-trim`, `make-sequence-table`, `remove-bimera-denovo`, `seq-table-to-tsv`, `seq-table-to-fasta` — **the pipeline closes** |
| Apr 17 | `sample`, `errors-from-sample` |

The build order is the C++ dependency graph, bottom up: data structures and
k-mers before alignment, alignment before clustering, clustering before the
driver. Everything through `evaluate`/`error`/`chimera` landed on a single day,
which is only possible because each layer could be written against a fixed,
already-translated one below it.

### The pivot

A week after the pipeline closed, the emphasis changed:

- **Apr 22** — per-iteration cluster diagnostics and the first comparison
  scripts. Instrumentation before investigation.
- **Apr 22** — `align_vectorized` DP correctness fixed and banded alignment
  enabled ([#2](https://github.com/HPCBio/dada2-rs/issues/2)); `omega_c` default
  corrected to R's `1e-40`.
- **Apr 25** — [#4](https://github.com/HPCBio/dada2-rs/issues/4) filed: a 7%
  transition-count gap and a divergent error matrix against R.

That issue is where the project stopped being a port and started being a
reimplementation held to a standard. The habit every findings page now assumes —
score on ASV-level concordance, not on whether the output looks reasonable —
dates from there.

## What this dictates

- **Know which layer you are in.** A difference from R inside the translated C++
  layer is a bug. A difference in the rewritten workflow layer is a design
  choice, and needs documenting rather than fixing — `remove-primers` is the
  canonical example.
- **The C++ mapping is stable enough to navigate by.** Looking for the R
  counterpart of a Rust module is a table lookup, and the [algorithm
  pages](algorithm-core.md) follow the same structure.
- **`dada_uniques` is the contract.** Anything that wants to consume the
  denoising core — a library user, a future R binding — goes through it, which
  is why its signature has stayed stable while everything under it changed.
- **The auto-vectorisation bet is still being cashed.** Dropping hand-written
  SIMD is what allowed the kernel to be rewritten repeatedly for
  [banding](findings/band-size-platform-defaults.md),
  [rolling `d16`](findings/compare-screen-vs-align.md) and
  [an alternative backend](findings/wfa-viability.md) without maintaining four
  intrinsic paths.
