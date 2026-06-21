# wfa-ffi-oracle

A **throwaway oracle** for the WFA ends-free + nonzero-match-score divergence
(upstream WFA2-lib [#102](https://github.com/smarco/WFA2-lib/issues/102)). It
exists to answer one question:

> Does the **latest C++ WFA2-lib** return the *optimal* ends-free alignment when
> the match score is nonzero (DADA2 uses match = +5), or does it reproduce the
> same suboptimality we see in the pure-Rust port?

## Results (WFA2-lib `bcf473a`, 2026-06) — #102 is **partially** fixed

Run via `./run.sh` (C harness `oracle.c` linked against the C++ static lib).
With `match=-5` the WFA penalty score equals the DADA2 score, so any result
below the known optimum is a genuine optimality miss, not a penalty-space tie:

| case | mode | optimum | C++ | verdict |
|------|------|---------|-----|---------|
| leading-gap 52nt | Mode 1 (free end-gap crediting) | 255 | 255 | ✅ fixed |
| `GCGG`/`CG` (the #102 thread reproducer) | Mode 1 | 10 | 10 | ✅ fixed |
| leading-gap 55nt | Mode 1 | 275 | 275 | ✅ fixed |
| trailing interior gap | Mode 2 (mismatch + free-end indel over interior gap) | 272 | **271** | ❌ still off by 1 |
| near-end interior gap | Mode 2 | 347 | **346** | ❌ still off by 1 |

**3/5.** Current C++ landed the leading/trailing free-end-gap fix (**Mode 1**,
the common case — including the issue's own `GCGG`/`CG`) but still exhibits
**Mode 2**: its CIGARs end `…X​M​D` / `…X​M​I`, i.e. it prefers `[mismatch +
free-end indel]` over the higher-scoring `[interior indel + trailing matches]`.

**Implication.** Porting the C++ ends-free fix into `wfa2lib-rs::termination.rs`
would remove the *bulk* (Mode 1) of our divergence, but **not** make WFA a clean
drop-in: Mode 2 survives even upstream, so WFA stays non-identical to NW and can
still flip the occasional low-abundance ASV. Option (c) stands — NW remains the
error-model backend; this porting is worthwhile only if WFA's *speed* ever
justifies it (currently it does not; see `project_wfa_edit_budget_cap_issue51`).

## Two harnesses in this dir

- **`oracle.c` + `run.sh`** — the **working** oracle: links the C++ WFA2-lib
  static lib directly (simplest path for an FFI oracle). `run.sh` clones+builds
  WFA2-lib at a pinned commit, compiles, and runs. This produced the results
  above.
- **`Cargo.toml` + `src/main.rs`** — an optional pure-Rust-binding stub (its own
  `[workspace]`, never built by the root). Kept as the skeleton if we later want
  the oracle as a Rust FFI crate instead of a C harness.

## Why this is isolated under `dev/` and never ships

- The production aligner is the **pure-Rust** `wfa2lib-rs` (a re-implementation,
  not an FFI binding). We deliberately avoided a C/C++ build dependency
  (cc/bindgen, cross-compile friction). See memory
  `project_wfa_dependency_wfa2lib_rs`.
- This crate uses **FFI into the C++ library on purpose**, but *only as a test
  oracle*. Nothing here is a dependency of dada2-rs. It is a standalone Cargo
  project (its own `[workspace]`), so the root `cargo build`/`cargo test` never
  touches it and it cannot leak a C++ toolchain requirement into CI or releases.

## The decision it gates

The #102 bug is **not** the match-score setting (the Rust port already exposes
`match_` with Eizenga reward-conversion, and our global score already matches
NW's optimum). The bug is that the ends-free **termination** ignores the match
reward — `wfa2lib-rs/src/termination.rs::termination_endsfree` is purely
length-based. So:

1. **If the latest C++ lib gets these reproducers right** → the fix is real and
   we port the match-aware termination logic into the pure-Rust
   `termination.rs` (shipping crate stays 100% Rust at runtime).
2. **If the C++ lib is *also* wrong** → no point porting; #102 is unfixed
   upstream and WFA stays experimental with NW as the error-model backend.

**Do not port any termination code until this oracle confirms (1).**

## The reproducers

Both come from `dada2-rs` `src/nwalign.rs` (DADA2 scoring: match +5, mismatch
-4, gap -8, ends-free on both sequences):

| case | pattern / text | optimal (NW `align_endsfree`) | pure-Rust WFA |
|------|----------------|-------------------------------|---------------|
| leading-gap (Mode 1) | `s2 = s1` minus leading base, 52nt | **255** (free leading end-gap, all 51 cols match) | 247 (penalized internal gap) |
| #102 minimal | pattern `GCGG`, text `CG` | optimal `GCGG / -CG-` | suboptimal `-GCGG / CG--` |

`src/main.rs` hard-codes these and prints, per case: the C++ WFA score + CIGAR
vs the known optimum, and PASS/FAIL (PASS = C++ matches the optimum = the fix is
real).

## Wiring the FFI binding (TODO)

`src/main.rs` currently runs as a **stub** (no external dep) so it compiles and
documents the expected values today. To make it a live oracle, add a binding to
the C++ WFA2-lib and replace the `wfa_endsfree_score()` stub. Options, in order
of preference:

1. An existing Rust binding crate to **WFA2-lib (C++)** that exposes ends-free +
   settable match score (verify it tracks a recent enough WFA2-lib that includes
   any #102 fix; the older crates.io `libwfa` binds the *original* C WFA and
   predates Eizenga match support — not usable here).
2. A minimal hand-written `bindgen`/`cc` shim against a local WFA2-lib checkout
   (pin the exact commit/tag tested in this README when results are recorded).

Record the WFA2-lib version/commit actually tested next to the results.

## Run

```
cd dev/wfa-ffi-oracle
cargo run            # once a binding is wired; prints PASS/FAIL per reproducer
```
