# Parameters: `setDadaOpt()` parity

dada2-rs aims to mirror R DADA2's algorithm parameters faithfully. Every
**result-affecting** option from R's `setDadaOpt()` (the `dada_opts` defaults in
[`dada.R`](https://github.com/benjjneb/dada2/blob/master/R/dada.R)) is exposed as
a CLI flag with the **same default value**, so a stock dada2-rs run reproduces a
stock R run.

## Equivalency table

| R `dada_opts` | R default | dada2-rs flag | dada2-rs default | match |
|---|---|---|---|:---:|
| `OMEGA_A` | 1e-40 | `--omega-a` | 1e-40 | ✅ |
| `OMEGA_P` | 1e-4 | `--omega-p` | 1e-4 | ✅ |
| `OMEGA_C` | 1e-40 | `--omega-c` | 1e-40 (dada) / 0 (learn) | ✅ † |
| `DETECT_SINGLETONS` | FALSE | `--detect-singletons` | false | ✅ |
| `USE_KMERS` | TRUE | `--no-kmer-screen` | on | ✅ |
| `KDIST_CUTOFF` | 0.42 | `--kdist-cutoff` | 0.42 | ✅ |
| `MAX_CONSIST` | 10 | `--max-consist` (learn) | 10 | ✅ |
| `MATCH` | 5 | `--match` | 5 | ✅ |
| `MISMATCH` | -4 | `--mismatch` | -4 | ✅ |
| `GAP_PENALTY` | -8 | `--gap-p` | -8 | ✅ |
| `BAND_SIZE` | 16 | `--band` | 16 | ✅ |
| `MAX_CLUST` | 0 | `--max-clust` | 0 | ✅ |
| `MIN_FOLD` | 1 | `--min-fold` | 1 | ✅ |
| `MIN_HAMMING` | 1 | `--min-hamming` | 1 | ✅ |
| `MIN_ABUNDANCE` | 1 | `--min-abund` | 1 | ✅ |
| `USE_QUALS` | TRUE | `--use-quals` | true | ✅ |
| `HOMOPOLYMER_GAP_PENALTY` | **NULL** | `--homo-gap-p` | → `--gap-p` | ✅ ‡ |
| `GREEDY` | TRUE | `--greedy` | true | ✅ |
| `PSEUDO_PREVALENCE` | 2 | `--pseudo-prevalence` | 2 | ✅ |
| `PSEUDO_ABUNDANCE` | Inf | `--pseudo-min-abundance` | off (= Inf) | ✅ |
| `VECTORIZED_ALIGNMENT` | TRUE | — (always on) | — | ◻︎ perf |
| `GAPLESS` | TRUE | — (always on) | — | ◻︎ perf |
| `SSE` | 2 | — (Rust autovectorizes) | — | ◻︎ N/A |

These flags apply to the denoising and error-learning subcommands (`dada`,
`dada-pooled`, `dada-pseudo`, `learn-errors`, `errors-from-sample`).
`remove-bimera-denovo` additionally honors `--match` / `--mismatch` / `--gap-p`,
matching R's `removeBimeraDenovo` (which respects the global `GAP_PENALTY`).

## Notes

- **† `OMEGA_C`** — the base `dada_opts` default is `1e-40`, but `learnErrors()`
  overrides it to `0` internally (error inference shouldn't merge by abundance
  p-value). dada2-rs encodes the same split: `1e-40` for the `dada*` commands,
  `0` for `learn-errors`. It is **not** inherited via `--inherit-err-params`.

- **‡ `HOMOPOLYMER_GAP_PENALTY = NULL`** — in R, `NULL` means "score homopolymer
  gaps with the normal `GAP_PENALTY`." dada2-rs reproduces this: when
  `--homo-gap-p` is unset it **falls back to whatever `--gap-p` resolves to**
  (so `--gap-p -4` alone yields `homo_gap_p = -4`), rather than pinning a literal
  default.

- **Scalar `MATCH` / `MISMATCH`** — R's full 4×4 `SCORE_MATRIX` is commented out
  in `dada.R`; the active path uses the scalar `MATCH`/`MISMATCH` (uniform +5
  diagonal / −4 off-diagonal). dada2-rs matches that active path.

- **`VECTORIZED_ALIGNMENT` / `GAPLESS` / `SSE`** — deliberately **not** exposed.
  These are implementation/performance toggles, not algorithm parameters:
  dada2-rs always uses its vectorized banded aligner, and `SSE = 2` (R's SSE2
  selection) is moot since the Rust kernel autovectorizes (with AVX2/AVX-512
  available via a `release-native` build). They are result-equivalent.

## Beyond `setDadaOpt()`: `--kmer-size`

`KMER_SIZE` is **not** a `setDadaOpt()` option — R fixes the k-mer screen at
`k = 5` in its C layer. dada2-rs exposes it as a tunable flag, **`--kmer-size`
(3–8, default 5)**. The default matches R exactly, but the knob lets PacBio HiFi
runs trade memory for speed (e.g. `k = 7`); see the
[PacBio walkthrough](walkthrough-pacbio.md). This is the one place dada2-rs goes
*beyond* R's parameter surface.

## Experimental: WFA alignment backend

!!! warning "Experimental — developer build only"
    These flags require a source build with `--features wfa` (see
    [Installation](installation.md#developer-build-experimental-wfa-backend)).
    A default build — and the published crate — is Needleman-Wunsch only and
    **errors** if `--align-backend wfa2` is selected. Tracked in
    [issue #63](https://github.com/HPCBio/dada2-rs/issues/63).

dada2-rs can optionally route the ends-free alignment through a
[wavefront-alignment](https://github.com/smarco/WFA2-lib) (WFA) backend instead
of Needleman-Wunsch. It is **ASV-equivalent** to NW on the tested Illumina and
PacBio HiFi datasets, but the per-pair alignments are *not* byte-identical to NW
(a known upstream ends-free crediting difference).

| Flag | Default | Description |
| --- | --- | --- |
| `--align-backend` | `nw` | Alignment backend: `nw` (Needleman-Wunsch) or `wfa2` (experimental WFA). |
| `--wfa-max-edits` | `50` | WFA edit-budget cap (issue #51). A pair needing more than this many edits aborts WFA and falls back to the banded NW path (NW-identical for that pair), bounding WFA's cost on divergent pairs that slip past the k-mer screen. `0` = unbounded; ignored for the `nw` backend. |

The edit budget is an absolute edit count, not a fraction of read length:
denoising only aligns near-identical reads (~99.9% identity), so real
error-copies stay a few edits apart regardless of read length. See
[issues #49 and #51](https://github.com/HPCBio/dada2-rs/issues/51) for the
performance rationale and parity analysis.
