# The initial port: what was translated, what was rewritten

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
project's first months. The **short-term** goal was a rough but working port —
FASTQ in, sequence table out. Correctness against R was a **mid-term** goal, not
an afterthought, and the likely obstacle was named in advance: the error model,
precisely because it is a port of R's LOESS rather than of C++. That suspicion
was borne out, and working through it is the
[LOESS page](findings/loess-error-model-correctness.md) — a rare case in this
project of a predicted problem turning out to be the actual one.

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
  writes plain loops and lets LLVM auto-vectorise to NEON or AVX. This is
  [deferred, not decided](#hand-written-simd-deferred-not-rejected), and it is
  why "does our kernel match theirs cell for cell?" was an open question for the
  first months.
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
| Apr 7 | `summary` — the tutorial's first step, before any algorithm |
| Apr 8 | `derep`; then the whole C++ core: foundations → `nwalign` → `cluster`/`filter` → `Rmain` → `taxonomy` → `evaluate`/`error`/`chimera` |
| Apr 9 | `error_models`, subsampling |
| Apr 13 | `learn-errors` — the iterative self-consistency loop |
| Apr 14 | `dada`, `merge-pairs` |
| Apr 15 | `filter-and-trim`, `make-sequence-table`, `remove-bimera-denovo`, `seq-table-to-tsv`, `seq-table-to-fasta` — **the pipeline closes** |
| Apr 17 | `sample`, `errors-from-sample` |

Two different orders are at work here, which is worth separating.

**The core followed the C++ dependency graph, bottom up:** data structures and
k-mers before alignment, alignment before clustering, clustering before the
driver. Everything through `evaluate`/`error`/`chimera` landed on a single day,
which is only possible because each layer could be written against a fixed,
already-translated one below it.

**The subcommands followed the DADA2 SOP tutorial**, in the order a user meets
them. That is why `summary` came first, a day before any algorithm: the
tutorial's early step is `plotQualityProfile()` — survey raw read quality and
expected errors, *then* choose truncation and filtering parameters. A port that
cannot tell you what your reads look like cannot be followed along with the
tutorial, whatever else it can do. `scripts/plot_quality_profile.R` completes
that step by rendering the figure from `summary`'s JSON.

The same logic explains the rest of the sequence — `derep` next, then the
denoising core, then `merge-pairs`, `make-sequence-table` and
`remove-bimera-denovo` closing the pipeline on Apr 15 in exactly the order the
tutorial walks through them.

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

## Two early decisions still open

### JSON as the intermediate format

Every subcommand hands work to the next through JSON. That was chosen early and
for one reason — **simplicity while prototyping** — and it stuck, which is the
usual fate of a format decision made before anyone knows the workload.

It has since been pushed at from the efficiency side rather than replaced:
gzip output (`--gzip`), a single-parse reader that removed 63% of the pooled
derep load and 8.6% of pooled wall time
([#133](https://github.com/HPCBio/dada2-rs/issues/133)), and dropping resident
intermediates. The format itself has not changed.

Whether something denser — Parquet, Avro, bincode — would be better is
[#1](https://github.com/HPCBio/dada2-rs/issues/1), open since the prototyping
phase and deliberately unhurried. The case for revisiting is a measured one, not
a tidiness one: it needs a workload where serialisation is demonstrably the
constraint. Note that JSON's readability has been load-bearing more than once —
several findings on this site were possible because an intermediate artefact
could be opened and inspected without tooling.

### Hand-written SIMD: deferred, not rejected

Dropping the C++ intrinsic paths is not a permanent verdict, but it is a
well-supported default. The repeated finding has been that the hot kernels are
**not execution-bound where SIMD would help**: `b_compare`'s DP kernel is
[memory-bandwidth-bound above ~48 threads](findings/compare-screen-vs-align.md),
the pooled scans are
[bandwidth-bound rather than op-count-bound](findings/shuffle-build-scan.md), and
the serial store turned out to be
[a cache-line problem, not a compute one](findings/compare-store-scan.md).
Auto-vectorisation has been competitive everywhere it has been measured.

Reopening it would mean a **from-scratch survey against the current code**, not a
re-reading of the old ones — the kernel has changed substantially since those
surveys, and the conclusion is only as current as the code it was drawn on.

Worth watching as precedent: `wfa2lib-rs`, the crate behind
[our WFA backend](findings/wfa-viability.md), has
[added hand-rolled SIMD](https://github.com/COMBINE-lab/wfa2lib-rs/commit/bc838bba260d2f154c5bf045a1c863a5e47a8930)
— NEON processing 4 diagonals per iteration, AVX2 8, with runtime feature
detection and a scalar fallback. Their reported gains on Apple Silicon are
**1.7× on edit distance, 1.4× on gap-linear, ~1.06× on affine**, with affine-2p
reaching parity with the C reference. Two things to take from it: hand-rolled
intrinsics in a Rust bioinformatics port are tractable and bounded in size, and
their gains are largest on the simplest scoring model — which is not the regime
our banded NW spends its time in.

One detail there corroborates our own result. They leave the **extend kernel
scalar deliberately**, judging SIMD unlikely to help it on algorithmic grounds.
That is exactly the kernel our profiling blamed for the PacBio WFA slowdown, and
it is why the fix was an
[edit-budget cap rather than a faster extend](findings/wfa-viability.md).

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
- **Hand-written SIMD is deferred, not rejected** — see below. Its absence is
  what allowed the kernel to be rewritten repeatedly for
  [banding](findings/band-size-platform-defaults.md),
  [rolling `d16`](findings/compare-screen-vs-align.md) and
  [an alternative backend](findings/wfa-viability.md) without maintaining four
  intrinsic paths, which is a real benefit independent of the performance
  question.
- **JSON is a prototyping decision that stuck.** Revisit it as a measurement
  ([#1](https://github.com/HPCBio/dada2-rs/issues/1)), not as a cleanup.
