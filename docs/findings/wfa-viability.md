# WFA as an alignment backend: viable, but not a drop-in

**Verdict: WFA is a working experimental backend and stays experimental.** It is
competitive or faster than the banded Needleman-Wunsch path on every workload we
measured once two bugs were fixed — but it is **not byte-identical to NW**, and
the reason is an upstream paradigm mismatch we cannot fix from this repository.
Needleman-Wunsch remains the default and the error-model backend.

Two independent results, which are easy to conflate:

- **Speed — solved.** The "WFA is slow on PacBio" problem was two bugs, not the
  algorithm. Banding it and capping its edit budget turn a 2.3× slowdown into a
  result that *beats* NW (141.8 s → 59.1 s vs NW's 61.7 s). The residual slowdown
  was **k=5-only**; at the recommended PacBio k=7 WFA was already 1.8× faster
  uncapped.
- **Correctness — open, and it does not close by itself.** WFA under-credits free
  end-gaps when the match score is non-zero. This is invisible on the concordance
  fixtures and visible at scale (Jaccard ~0.92–0.95 vs NW, unbounded). The edit
  cap does not mask it, and no configuration change fixes it.

## Why WFA at all

The alignment paradigm is fixed to **ends-free global** alignment by the
error-model contract: `pval.rs` consumes *positional* substitutions, so any replacement has to solve the
same global problem and emit the same `Sub` contract, not merely produce a good
alignment. WFA qualifies — it is exact for its scoring model, and its cost is
O(n·s) in the *edit distance* rather than O(n·m) in the read lengths. Denoising
compares near-identical sequences, so s is small where n·m is not. That is the
whole thesis, and it survives: the WFA kernel is **~60× faster than NW** on a
near-identical 1.5 kb pair.

## The dependency

[COMBINE-lab/wfa2lib-rs](https://github.com/COMBINE-lab/wfa2lib-rs) — a
from-scratch **pure-Rust** WFA, BSD-3-Clause. Pure Rust matters more than it
sounds: it keeps `cc`/`bindgen` and a C toolchain out of the cross-compilation
path (musl, aarch64, CI containers).

It is **not on crates.io**, so it enters as a pinned git dependency on the HPCBio
fork, which is why the backend is behind an off-by-default `wfa` Cargo feature —
crates.io will not accept a crate with a git dependency ([#63](https://github.com/HPCBio/dada2-rs/issues/63)).
It also raises the effective MSRV to **1.91**.

**`libwfa` was considered and rejected** — do not resurface it. It is the one WFA
crate actually published on crates.io, and it is worse on every axis that
mattered: FFI bindings to the C WFA (~3,600 lines of vendored C, reintroducing
the toolchain problem), unmaintained since 2020, and it **pre-dates WFA2's
Eizenga non-zero-match support** — so DADA2's `match = +5` scoring cannot be
expressed on it at all. Being on crates.io did not outweigh that.

Three packaging fixes went upstream rather than into a divergent fork:
feature-gating the CLI deps, a public `set_alignment_scope` setter, and declaring
the license.

## Speed: two bugs, in order

### Bug 1 — WFA was running unbanded

Our adapter ran WFA globally optimal while NW was band-limited. Since WFA's cost
grows with edit distance, cost crosses NW's at **~130 edits / 1450 bp** and keeps
going (1.97× at 200 edits) — and PacBio mock-community cross-species pairs sit
past that crossover. This was also a latent *semantic* difference: a globally
optimal alignment is not the same object as a band-limited one.

Passing the band through as `BandedStatic ±band` ([#52](https://github.com/HPCBio/dada2-rs/pull/52))
took PacBio `learn-errors` from **124 s → 68 s** (4.3× → 2.4× NW).

### Bug 2 — the residual 2.4× is a k=5 artifact, not a WFA property

Profiling attributed the remainder to genuine WFA `extend` compute, not per-call
overhead. The fix is an **edit-budget cap**: `set_max_alignment_steps` aborts a
pair that exceeds the budget, and it re-aligns on the banded NW path —
NW-identical for exactly those pairs. PacBio concordance fixture, 6740 filtered
reads, band 32, 1 thread, `learn-errors`:

| k-mer size | NW | WFA uncapped | WFA cap ≈ 40 edits |
|---|---|---|---|
| **k=5** (R's `KMER_SIZE`) | 61.7 s | **141.8 s** (2.3× slower) | **59.1 s** (beats NW) |
| **k=7** (recommended PacBio) | 13.1 s | 7.3 s (1.8× faster) | 8.7 s |

The k=7 column is the interesting one: **the slowdown never existed at the
setting we actually recommend.** A looser k=5 screen admits divergent
*non-error-copy* pairs, and WFA pays full O(n·s) on them; k=7 removes them
before they reach the aligner. Banding alone does not bound this — k=5 uncapped
was 141.8 s *with* band 32.

The budget is an **absolute edit count, not a fraction of read length**, because
error copies stay a few edits apart regardless of how long the read is. Default
50 edits, exposed as `--wfa-max-edits` on every aligning subcommand and recorded
in the error-model JSON alongside `--align-backend`.

### Illumina is the opposite case

MiSeq SOP, 20 samples, paired V4, ~250 bp, per-sample denoising, 4 threads:

| step | NW | WFA |
|---|---|---|
| learn-errors (F+R) | 35.5 s | 33.1 s |
| dada (F+R) | 8.1 s | 7.9 s |
| **total** | **43.6 s** | **41.1 s** (~6% faster) |

Short, similar reads are WFA's sweet spot; the per-call overhead that hurt at
1500 bp is a net win at 250 bp.

## Correctness: the ends-free divergence

**Mechanism.** DADA2 *maximises* a score with `match = +5`; WFA *minimises* a
penalty where match costs 0. Under free end-gaps the number of scored columns
changes, so the `+5`-per-match reward moves the true optimum — and WFA's
cost-model pruning cannot see that. It under-credits free end-gaps. This is
upstream [WFA2-lib #102](https://github.com/smarco/WFA2-lib/issues/102), not a
bug in the adapter.

Two failure modes: a free leading/trailing end-gap is not credited (Δ ≈ one gap
penalty), or a free end terminates prematurely (a mismatch plus a free end-gap
is chosen over a higher-scoring interior gap). Minimal reproducer, kept as the
deterministic test `nwalign.rs::wfa_endsfree_known_divergence`: with `s2 = s1`
minus its leading base, `align_endsfree` credits a free leading end-gap and
scores 255, while WFA places a penalised internal gap and scores 247.

**The gap model is not implicated.** Running WFA *global* against the scalar
global `align_standard` on the same 10k random pairs isolates it:

| comparison | disagree |
|---|---|
| WFA global vs `align_standard` | 9 / 10 000 (0.09%) |
| WFA ends-free vs `align_endsfree` | 550 / 10 000 (5.5%) |

~98% of the divergence is free-end-gap crediting. This is what rules out affine
gap scoring as a mitigation: affine makes single-base gaps *costlier*
(`open + extend` ≥ linear), which pushes WFA further toward the
mismatch-plus-free-end-gap resolutions that already lose — and it diverges from
DADA2's linear-gap reference by construction, so it cannot improve parity.

### What it costs in ASVs

This is where the finding is easy to misread, because the answer depends on the
size of what you measured.

- **Small fixtures: nothing.** Illumina forward reads, byte-identical ASV set
  (360 = 360, churn 0). PacBio HiFi full-length 16S, a **full** swap since HiFi
  leaves `use_homo` false so every alignment goes through WFA: 11 = 11, churn 0.
  Full pipeline at k=5 with the cap, against the static R reference: **93 = 93
  ASVs, recall 1.000, precision 1.000, count-correlation 0.994** — identical to
  the NW baseline.
- **One divergent call on a real sample.** Illumina reverse reads: 304/305
  shared, WFA splitting one extra abundance-4 ASV in F3D144. The mechanism was
  traced rather than assumed — both backends find the same Hamming-11 nearest
  neighbour, but WFA's different ends-free alignment of those 4 reads shifted the
  substitution counts, dropping `birth_pval` to 6.95e-42, below `OMEGA_A` in WFA
  and not in NW, flipping one partition decision.
- **At scale: Jaccard ~0.92–0.95 against NW**, unbounded, as low-abundance ASV
  over-calls.

The divergent pairs are *low-edit* pairs. That is why the concordance gates pass
and the effect still exists: small fixtures do not contain enough of the
distribution to show it, the same trap documented in
[the minimizer work](minimizer-screening.md).

## Do not retry

- **The edit cap as a correctness fix.** It bounds cost, not correctness. The
  diverging pairs finish far under budget, so they never take the NW fallback.
- **BiWFA.** An O(s) memory optimisation of the same paradigm; it inherits the
  identical suboptimality.
- **Swapping pattern and text**, as suggested upstream. Tested on 20k random
  pairs: "pattern = shorter" fixed 8 and broke 53. Each argument order is wrong
  in a *different* place, so it does not generalise.
- **Affine gap scoring** — see above.
- **The `libwfa` crate** — see above.

## What this dictates

1. **NW stays the default and the error-model backend.** WFA stays behind
   `--align-backend wfa2`, a source build with `--features wfa`, and an
   EXPERIMENTAL label. It is not a candidate for promotion while it changes
   low-abundance calls.
2. **A real fix is a match-aware wavefront cost** — fork work on `wfa2lib-rs`,
   not configuration and not adapter work. Everything reachable from this side
   has been tried.
3. **Both WFA gates are in CI.** Illumina and PacBio k=5 (at the default cap)
   both run the WFA backend through the concordance harness, so a regression
   surfaces without anyone remembering to check.
4. **The screen governs WFA's worst case, not the aligner.** The k=5/k=7 table is
   a statement about which pairs reach the aligner at all. Note the boundary:
   [`b_compare` is align-dominated](compare-screen-vs-align.md) on both platforms
   under the *default NW* path, and the k-mer screen closes as a perf target
   there. These are compatible claims about different things — the screen sets
   WFA's cost *distribution*; it is not itself where the default path's time
   goes.
5. **Packaging is the other blocker, and it is not technical.** Until
   `wfa2lib-rs` reaches crates.io, the backend cannot ship in a released crate
   and the fork needs maintaining against upstream C++ WFA2. That is a standing
   cost to weigh before any promotion decision
   ([#51](https://github.com/HPCBio/dada2-rs/issues/51)).

## Tracking

[#49](https://github.com/HPCBio/dada2-rs/issues/49) (paradigm evaluation, open —
also holds the separate pair-HMM avenue),
[#50](https://github.com/HPCBio/dada2-rs/pull/50) (backend, merged),
[#51](https://github.com/HPCBio/dada2-rs/issues/51) (long-read optimisation),
[#52](https://github.com/HPCBio/dada2-rs/pull/52) (banding, merged),
[#63](https://github.com/HPCBio/dada2-rs/issues/63) (crates.io feature gate).
Flags and build instructions: [Parameters](../parameters.md#experimental-wfa-alignment-backend)
and [Installation](../installation.md#developer-build-experimental-wfa-backend).
