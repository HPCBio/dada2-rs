# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html) (pre-1.0, so
minor versions may carry breaking changes).


## [Unreleased]

### Added
- `--indel-events` and `--indel-opportunities` (Diagnostics, experimental) on
  `dada`, `dada-pooled` and `dada-pseudo` (#291): indel events of every unique
  against its final centre, and read-weighted substitution, insertion and
  deletion rate tables, from the same alignments as `trans`. Denoising output
  is unchanged. `dev/indel_spikein.py` spikes known indel variants into a real
  run for evaluating indel-aware denoising.
- An indel-aware λ (Experimental) on `dada`, `dada-pooled` and `dada-pseudo`
  (#291): `--indel-rates` scores each indel between a unique and its centre
  from a prior run's `--indel-opportunities` table, with `--hp-ceiling`,
  `--indel-pair-window` and `--indel-end-zone`. Off by default; R-parity output
  is unchanged. `--indel-opportunities` gains `del_fq`, `ins_len` and `del_len`.
- `learn-errors --indel-rates-out` (Experimental) learns the indel rates for
  `dada --indel-rates` in a second self-consistency loop after the
  substitution model, which stays R's; `--indel-history` keeps every
  iteration's tables (#291). `--indel-opportunities` gains context, joint and
  per-base tables, and `scripts/plot_indel_errors.R` plots observed indel rates.

### Fixed
- `remove-primers` `--rc-primer-rev`, `--trim-fwd`, `--trim-rev`, `--orient`
  and `--compress`, and `filter-and-trim --compress`, can now be turned off
  with `=false` (#296). They defaulted to true but took no value, so none could
  be disabled, and the documented `--rc-primer-rev false` was a parse error.
  A reverse primer given in R's read orientation was therefore
  reverse-complemented twice and every read discarded. The bare flag still
  means true; defaults and output are unchanged.
- `assign-species --help` and its docs page say the match is an exact
  *substring* match, as it has been since #200; "exact match" read as
  full-sequence equality, the behaviour #200 removed.

## [0.5.0] - 2026-10-19

The headline of this release: fidelity to R DADA2 results
* `dada`, `dada-pseudo` and `dada-pooled` now reproduce R's ASV tables cell for cell 
  on the Illumina and PacBio (k=5) concordance fixtures
* CI now holds as an **exact gate** on all fixtures directly comparable to R. 
* Pooled MiSeq (362 samples), NovaSeq ITS2 (30), and NovaSeq 16S () give R's 
  ASV set in both read directions. 

Several of the fixes that got there change results; changes to results, flags or 
output are marked **BREAKING** below.

### Added
- `--end-model ends-free|free-3p|global` (Experimental) on `learn-errors`,
  `errors-from-sample`, `dada`, `dada-pooled` and `dada-pseudo` (#289). The
  default, `ends-free`, is R DADA2's model; `free-3p` charges leading gaps
  (reads anchored after a trimmed forward primer), `global` charges both ends.
  The model is recorded in the error model (`params.end_model`) and `dada`
  warns on a mismatch. Chimera detection keeps R's model.
- `reference-eval` subcommand: aligns an ASV set (FASTA or `dada` /
  `dada-pooled` JSON) against a truth set, classifies each ASV TP / near / FP,
  reports recall and precision, joins each ASV's birth p-value, and with
  `--non-chimeric` cross-tabulates chimera calls against reference class (#91).
  This command is still considered a work in progress.
- `summary` detects binned quality scores: an always-on `quality_metrics` block
  (distinct Q values and counts, `is_binned`, `binned_scores`), plus
  `--binned-threshold` and `--report`. New `summary-merge` subcommand unions
  the per-sample results into a run-level bin set, flags samples missing a bin,
  and validates against a declared scheme with `--expected-bins`.
- `learn-errors` accepts derep / sample JSON (`.json`, `.json.gz`) as well as
  FASTQ, as `dada` does (#68).
- `--metrics-json PATH` on `dada`, `dada-pooled` and `dada-pseudo`: phase times,
  shuffle / bud / p-update counters, resident footprint, screen occupancy,
  derep split, peak RSS and, for `dada-pooled`, per-sample `pipeline.inputs`
  counts (#162, #178). `--metrics-attribution` adds per-comparison timings
  (slows the run slightly). A field the run did not measure is omitted, never `0`.
- `source_fastq` field in `dada`, `dada-pooled` and `dada-pseudo` output: the
  input FASTQ, or the FASTQ a derep JSON records (#111). `merge-pairs` checks
  provenance against it.
- `dada-pooled --pooled-record PATH` writes the merged pool (pooled abundances,
  global `map`, ASVs) as one record; `kdist-calibrate --from-dada-pooled` reads
  it and screens the pool as one population.
- `dada-pooled --cluster-trace` (with `--trace-no-members`,
  `--trace-min-abund`). Traces record `nraw`, `omega_a` and `omega_p`, so a
  prior birth's abundance-scale p-value (`birth_pval * nraw`) can be recovered.
- `dada-pseudo --reestimate-err-between-rounds` (opt-in) re-fits the error model
  from round 1's transitions before round 2, emulating what R's
  `pool = "pseudo"` does (#100). Recorded as `reestimated_err_between_rounds`
  in the output params. Off by default.
- `dada-pooled --gzip` writes the per-sample JSON files gzip-compressed as
  `{sample}.json.gz` (#70). `dada` and `dada-pseudo` gain the same `--gzip`
  flag, for consistency across the per-sample denoisers.
- Experimental minimizer pre-alignment screen, `--screen-backend minimizer`
  (`--minimizer-k`, `--minimizer-w`, `--screen-audit`); off by default, and the
  k-mer path is bit-identical. On PacBio HiFi it matches the k-mer screen
  exactly at cutoff >= 0.45 and is 18.7% faster; pooled ITS2 is 28-30% faster
  with 18 of 3,028 ASVs churned. `kdist-calibrate --derive-cutoff` derives the
  cutoff. The index is chosen by a run-time probe, which warns when its sample
  was unrepresentative.
- Opt-in `DADA2RS_GAPLESS_X86=1` reproduces x86-64 R's gapless test, which also
  fires for sequences of unequal length, for parity runs against an x86-64 R
  reference (#277). With it, pooled PacBio matches x86-64 R exactly (0 of
  62,585 cells). It changes results, so it warns on every run.
- `DADA2RS_MEMBER_ORDER=insertion|sorted|shuffle:<seed>` diagnostic gate for
  measuring the member-order floor (#157). Non-default arms change results and
  always warn; an unrecognised value is fatal.
- `--verbose` now reports the **resolved** values of the `DADA2RS_*` development
  tuning gates, and any `DADA2RS_*` variable that is set but not recognised is
  warned about **regardless of `--verbose`** (#145). Retired spellings carry an
  alias naming the replacement. Motivated by two full soil-pool A/B runs spent
  comparing a configuration against itself after `DADA2RS_SHUFFLE_CARRY` was
  renamed to `DADA2RS_SHUFFLE_NO_CARRY`: both arms ran carry-on, and came back
  byte-identical and within 1%, which is what a *correct* null looks like.
- `--verbose` progress lines for pooled `dada`: per-window phase deltas,
  effective cores, map parallel efficiency and align fraction, every
  `DADA2RS_PROGRESS_SECS` seconds (default 30, `0` disables) (#150). Totals hid
  that occupancy climbs 22 → 42 effective cores of 64 across a single run.
- `scripts/learnerrors_to_dada2rs.R` now writes R's `trans` transition counts
  when converting a full `learnErrors()` object, so `dev/compare_error_models.py`
  can check a dada2-rs model against R's in one step (identical transition
  totals, then rates in log10) (#271). It writes 17 significant digits (was 15),
  so converted models round-trip exactly (#277). Inference ignores `trans`;
  converted models load and denoise exactly as before. A bare error matrix
  still converts without it.
- `scripts/plot_errors.R` gains an `optparse` CLI (`--nti`/`--ntj`, layer
  toggles, `--mass`, `--size-trans`, `--mass-shading`, `--min-mass`) and colours
  points by share of observations, for binned-quality data. Zero rates are
  blanked instead of plotting as `-Inf`.

### Changed
- **BREAKING (results): `dada-pooled` now breaks abundance ties in the merged
  table by first appearance across the inputs, as R's `combineDereps2` does**
  (#260). It previously sorted ties by sequence. This changes pooled results,
  and they now depend on input order, as R's do. Against R `dada(pool=TRUE)`,
  the pooled ASV set is now identical on MiSeq (362 samples) and NovaSeq ITS2
  (30), both read directions, where it differed by up to 12 ASVs. The PacBio
  HiFi residual is covered under Fixed (#219, #204, #277). `--pool-tiebreak
  lexical` restores the old, order-independent rule. Outputs record the rule
  (`params.pool_tiebreak`) and each sample's `pool_input_index`;
  `--pooled-record` lists the inputs in order.
- **BREAKING (results): the default LOESS fitting surface is now `interpolate`**
  (#205), matching the surface R's `loess()` — and therefore `loessErrfun` —
  uses. This changes error models and can shift ASVs and counts; a stock
  `learn-errors` run now agrees with R without extra flags. `--loess-surface
  direct` restores the old behaviour. The previous default was exact against
  `loess(surface = "direct")`, but R DADA2 never calls that: on a 95-sample
  pooled PacBio run the surface was the largest remaining term in the error
  model (median 1.088e-03 between surfaces, against 2.599e-04 for the k-mer
  screen), and on a 362-sample MiSeq run `interpolate` landed +14 reads from R
  against `direct`'s -1,699. Cost at ASV level is nil on both platforms.
- **BREAKING:** `--sample-names` on `make-sequence-table` and `merge-pairs` is
  now comma-separated, like `dada-pooled` and `dada-pseudo`. It was
  space-separated and swallowed any positional inputs that followed it.
- **BREAKING (non-default methods):** `remove-bimera-denovo`'s parent-abundance
  defaults now depend on `--method`, as in R (#211): `consensus` keeps 1.5 / 2,
  while `pooled` and `per-sample` use R's 2 / 8. On the MiSeq SOP, `per-sample`
  keeps 1,877 ASVs instead of 1,299 and `pooled` 402 instead of 399. Passing
  either flag explicitly still overrides.
- **BREAKING (stderr):** `--verbose` now reports run shape only. The compare /
  shuffle / bud / p-update attribution tables and projections moved to
  `--metrics-json` (#162); a pooled fixture goes from 64 to 35 stderr lines.
  `dada-pooled` no longer prints unnamed `[derep]` lines.
- `--errfun binned-qual` now errors when an anchor inside the observed quality
  range has no observations; rates across the gap were silently `1e-7` (#263).
  It warns when the observed minimum or maximum quality is not an anchor
  (#208), or when an anchor holds less mass than a nearby non-anchor quality.
  Each warning prints once per run.
- `assign-taxonomy --seed` now has a **fixed default**, so taxonomic assignment
  is reproducible out of the box (#187). There is no unseeded mode: entropy only
  makes a run unrepeatable, and varying the seed gives the same spread while
  staying reproducible. This follows `kdist-calibrate`, which already seeds by
  default; R DADA2 currently cannot do this (benjjneb/dada2#1115).
- `plot_quality_profile.R`, `plot_complexity.R` and `plot_expected_error.R` now
  parse flags with `optparse`, like `plot_errors.R` (#198). All four answer
  `--help`, accept `--out FILE` as well as `--out=FILE`, and reject a
  non-numeric `--width`/`--height`, which previously wrote a default-sized PDF.
  Existing flags and output are unchanged; `optparse` is now needed by all four.
  `learnerrors_to_dada2rs.R` also uses it now, for `--help` and a usage error on
  the wrong number of arguments; its two positionals are unchanged.
- `--version` and every output's `dada2_rs_version` now tag non-default builds:
  `+wfa`, `+rmath-ppois`, `+native` / `+cpu=<name>`, and the Cargo profile when
  it is not `release` (#289). A default `cargo build --release` is unchanged.
- In-line help is a one-line-per-flag reference card under fixed headings
  (help text −66%), and each subcommand links its page under `docs/commands/`
  (#168). Top-level `--help` groups subcommands into Workflow, Quality control,
  Export and Diagnostics (#83).
- Errors print as `dada2-rs: error[<ErrorKind>]: <message>` (exit 1) instead of
  the `io::Error` Debug dump, so tooling can match on the kind.
- `kdist-calibrate` is no longer hidden. An explicit `--max-pairs` with
  `--nearest-parent` is now rejected, since that mode ignores it.
- With `dada --sample-jobs` > 1, each bud round's progress is emitted as one
  line tagged with its sample (#172). Concurrent samples previously spliced
  partial lines together, and log filters dropped them. Serial and pooled
  output is byte-identical.
- `merge-pairs` and `make-sequence-table` order abundance ties as R does (first
  occurrence, stable sort). Output is now byte-reproducible across runs; it
  previously depended on per-process `HashMap` order (#240).
- The `DADA2_RS_PAR_GRAIN` tuning gate is renamed `DADA2RS_PAR_GRAIN` for
  consistency with every other gate (#145). The old spelling is still honoured
  and warns; it was the lone `DADA2_RS_` name, so anyone setting it from memory
  got silence and the default grain.
- Experimental WFA backend (`--features wfa`) pinned to wfa2lib-rs `96c1f23`.
  Its scores now equal NW's optimum everywhere tested. `DADA2RS_WFA_LINEAR` is
  a recognised gate (#51).
- The published crate also excludes the documentation sources (`docs/`,
  `mkdocs.yml`); they are published on Read the Docs.


### Deprecated
- `--loess-preset` on `learn-errors` and `errors-from-sample`. It bundled
  `--loess-surface`, `--loess-cell`, `--loess-max-rate` and `--loess-min-rate`,
  with `default` selecting `direct` and `r-dada2` selecting `interpolate`. Now
  that `interpolate` is the default, `r-dada2` is redundant and `default` names
  the non-default surface, which is actively misleading. The flag still works
  and maps onto `--loess-surface`, but warns; `--loess-surface` wins if both
  are given.


### Fixed
- **BREAKING (results): pooled `dada` matched R's abundance-sorted pool only by
  accident** (#219). The merged pool was built in first-seen order, which can
  leave cluster 0's centre where `b_bud` can never reach it. On 95-sample
  PacBio a 66,937-read organism R calls with 106,853 reads was absent. It is
  now sorted by abundance, as `combineDereps2` does.
- **BREAKING (results):** `calc_pA` returned 1.0 instead of R's 0 at zero
  expected reads, so reads R leaves unplaced were absorbed (#204). With
  the #219 fix, pooled PacBio pre-chimera ASVs go 2,810 → 2,791, equal to R. MiSeq
  moves by one read.
- **BREAKING (output):** reported ASV abundance now counts only the reads `map`
  assigns, as R's `clustering$abundance` does (#204). Reads failing `OMEGA_C`
  were counted before. On the MiSeq SOP, per-sample and pseudo counts now equal
  R's in 20/20 samples, both directions (was 0/20, ~2,500 reads over).
  `dada-pooled` per-sample outputs were already correct; its `_pooled.json` ASV
  list changes.
- `assign_center` breaks read-count ties toward the first raw, as R does
  (#239). It fires only when the top of a derep is tied: 4 of 724 MiSeq
  inputs, 2 with misassigned reads. R agrees with the fix on all 4.
- The final-substitutions pass no longer takes the gapless shortcut, matching
  R, which runs it without k-mers (#277). This affects `dada --aux-outputs`
  transitions and cluster quality only; learned error models are unchanged.
- Error models round-trip f64 exactly through JSON; 269 of 1,504 cells read
  back 1 ulp off before (#277).
- LOESS now matches R's `loess()` to ~1e-14 in every regime (#213, #215). The
  kd-tree bounding box is padded as R's `ehg126` pads it, the kd-tree leaf
  rules are R's (vertex sets now equal R's `kd$xi`), the local design is
  centred at the query, and `span > 1` widens the bandwidth. Sparse-quality
  fits (3-8 populated Q columns) were up to 0.5 log10 off; dense-data ASVs are
  unchanged.
- An error model that does not cover a read's highest quality is extended by
  repeating its last column, as R does, instead of panicking (#102). This
  warns on every run.
- Non-ACGT bases (`N`) are rejected with an error naming the sequence, base,
  position and fix (`filter-and-trim --max-n 0`) instead of panicking a worker
  thread (#101). `learn-errors` surfaces the first per-sample failure instead
  of only "All DADA runs failed".
- Multi-file inputs are checked for existence before the list lengths are
  compared, with a hint when a missing path looks like an unexpanded shell
  glob (#117). `merge-pairs` used to blame the wrong flag.
- `merge-pairs` no longer warns about provenance on every sample when `dada`
  ran from derep JSON (#111).
- Single-input `dada --metrics-json` and `--cluster-trace` now label the run
  with the output's sample name, which respects `--sample-name` and strips
  `.fastq.gz` (#67).
- `merge-pairs` output with two pairs merging to one sequence now sums their
  counts, as R's `getUniques` does; the last count used to win (#240).
- `build.rs` no longer forces a rebuild on every cargo command in a git
  worktree or outside a repository; a build with no `.git` embeds the plain
  version.
- **BREAKING (results): `merge-pairs` now finds and judges the overlap as R's
  `mergePairs` does** (#272). It aligned with the denoising scores (5 / −4 /
  −8), so on length-variable amplicons a longer, imperfect overlap at a wrong
  offset could outscore a short perfect one and the pair was rejected; it now
  uses R's merge scores (1 / −64 / −64 at `--max-mismatch 0`, 1 / −8 / −8
  otherwise). The overlap is counted as R's `C_eval_pair` counts it, from the
  column where both reads have begun, and accepted on R's rule (`--min-overlap`
  counts matching bases; `--max-mismatch` counts mismatches plus indels). On
  NovaSeq ITS2 this restores the 10 merged cells (93 reads) R kept and we
  dropped; merged MiSeq V4 output is unchanged.
- `assign-species` matched a query against a reference only when the two
  sequences were **equal in full**, so it returned nothing on any full-length
  reference — the ordinary case, since species FASTAs hold ~1.4 kb 16S and
  queries are sub-region amplicons (#200). It now matches by **containment**,
  as R does (`vcountPDict(...) > 0`). On the MiSeq SOP, 232 ASVs against Silva
  138.2: 0 species before, 16 after — per-ASV identical to R's `addSpecies`.
  The failure was silent: exit 0, valid JSON, a present-but-null Species column.
- **BREAKING (seeded output):** `assign-taxonomy` is now reproducible regardless
  of **input order and input set** (#187). Each sequence's bootstrap RNG is
  derived from the sequence (pinned MD5) instead of its position, so a
  sequence's classification depends only on itself, the reference and the seed.
  Shuffling 3,994 queries previously changed 7.6% of assignments at a fixed
  seed; it now changes none. Thread-count independence is unchanged. This is the
  dada2-rs counterpart of the long-standing R issue benjjneb/dada2#1115, from a
  different cause. **Seeded output moves once with this change**, since every
  sequence gets a new stream; unseeded runs were never reproducible and are
  unaffected.
- The `nalign` and `nshroud` counters no longer overflow on large pooled runs
  (#143). Both count pairwise comparisons (`nraw × nclusters`) and were `u32`,
  so any pool exceeding ~4.3 billion comparisons wrapped silently: a 1.23 M-
  unique NovaSeq soil 16S pool reported 3,015,581,030 alignments against an
  actual 11,605,515,622. The wrapped values reached the `ALIGN:` verbose line
  *and* the `nalign` / `nshroud` fields of `dada`, `dada-pseudo` and
  `dada-pooled` output JSON, so affected artifacts carry wrong values in those
  two fields. Nothing reads either counter — no clustering, error-model,
  p-value or table output was affected, and no re-run is needed to correct an
  archived result. Smaller pools (under ~4.3e9 comparisons) were never
  affected.
- LOESS error models no longer collapse to a uniform `1e-7` floor on input with
  five or fewer populated quality columns (#95). The tricube kernel zeroes the
  farthest neighbour in the local window, which left a `span = 0.75`,
  `degree = 2` fit short of the coefficients it needed; the local fit now falls
  back to the highest degree the surviving points can identify. Affected
  `--errfun loess` (and the `pacbio` loess path) on binned-quality data, where a
  handful of distinct Q values is normal. Fits with six or more populated
  columns are byte-identical to before. `loess_errfun` and `pacbio_errfun` now
  return `Result` and report an error rather than emitting an unusable uniform
  matrix when every transition fails to fit.
- `learn-errors` and `errors-from-sample` now gzip-compress their output when
  the `-o` path ends in `.gz`, for consistency with the other JSON-emitting
  subcommands (previously the file was written uncompressed) (#71).
- `kdist-calibrate --from-dada` now resolves derep inputs whose filenames were
  renamed by the pipeline (e.g. `{sample}.derep.R1.json.gz`): it falls back from
  the exact `{sample}.json[.gz]` match to a `{sample}.*.json[.gz]` scan,
  disambiguating an ambiguous prefix by the derep JSON's own `sample` field
  (#72).
- `kdist-calibrate --from-dada` now reads gzip-compressed `dada` output JSONs
  (`*.json.gz`, e.g. from `dada-pooled --gzip`); previously the gzip bytes were
  parsed as raw JSON and failed with `expected value, line 1 column 1`.
- `make-sequence-table` now reads gzip-compressed `dada` / `merge-pairs` inputs
  (`*.json.gz`, e.g. from `dada --gzip` / `dada-pooled --gzip`); previously the
  gzip bytes were parsed as raw JSON and failed with `expected value, line 1
  column 1`.
- Single-file JSON `-o` outputs now honor a `.gz` output path, routing through
  the shared gzip-aware writer instead of a raw write (matching
  `learn-errors`, #71): `summary`, `summary-merge`, `merge-pairs`, `remove-primers`,
  `filter-and-trim`, `make-sequence-table`, `remove-bimera-denovo`,
  `assign-taxonomy`, `assign-species`, and `dada`'s single-sample output.


### Performance
- `b_p_update` prefetches the `Raw` it will read 16 iterations ahead (#154),
  cutting the phase **34-38%** on four arms across two NovaSeq soil pools and
  **−1.0 to −6.8% of `run_dada`** depending on how large a share the phase is
  (5.1% on 16S read 1, 14.3% on ITS2 read 2). Byte-identical by construction —
  a prefetch has no architectural effect. The phase is memory-bound on a random
  gather, not compute-bound: over 91% of repricings take an early exit doing no
  arithmetic, so the cost is a cache miss on a 160-byte `Raw` read for ~20 bytes
  of it. Three layout fixes were priced first and all were worse; the only
  competitive one wins by fitting a 32 MB per-CCD L3 and therefore *loses* on
  the larger pool. x86_64 only; `DADA2RS_PUPDATE_PREFETCH=0` disables it.
- `Raw::e_minmax` moves off `Raw` into a dense `Vec<f64>` on `B`, parallel to
  `B::raw_cluster` (#147). `b_compare`'s serial store reads this value once per
  raw per call; `Raw` is 160 bytes, so as a field each read pulled a full
  64-byte cache line to use 8 bytes of it, and that strided read measured as
  **83–87% of the store scan**. Contiguous, the store falls from 18.8–19.1 to
  4.9–5.5 ns per raw-visit — **−71% of store time** — worth **−20.1 to −21.3%
  of `run_dada`** on NovaSeq soil 16S and −15.2% on ITS2, byte-identical on
  both pools across two replicates. Effective cores rise 31.4 → 38.4 (R1) and
  29.4 → 36.5 (R2) of 64.
- `b_compare`'s serial reduction is folded into the store loop it already
  shared a result vector with (#143): seven passes over that vector become
  one. Worth **−24 to −30% of `run_dada`** on two NovaSeq soil pools
  (16S and ITS2, both reads), byte-identical. The reduction had been six
  separate passes costing 337 s on soil 16S — 31% of `run_dada` — on a phase
  whose parallel map is already bandwidth-saturated and cannot absorb it.
- `b_shuffle` carries `compmax` across bud rounds instead of rebuilding it
  per bud (#139), deleting the build scan: 9,700 rebuilds collapse to 1.
  Measured **−45 to −48% of `run_dada`** on soil 16S and −15.4 to −15.9% on
  ITS2 (64 threads, post-#143), byte-identical across arms.
  `DADA2RS_SHUFFLE_NO_CARRY=1` restores the per-bud rebuild.
- Pooled `dada` shuffle and bud are incremental: a persistent candidate index
  and a per-cluster bud cache. The 384-sample MiSeq pool goes from 427 s to
  359 s (−16%), byte-identical; `dada` RSS rises ~5% (#85).
- The shuffle's move pass visits only dirty clusters (#132): −4.3 to −5.6% of
  `run_dada` on MiSeq, −1.8% on PacBio. The reconcile tests only changed
  candidates (#136): −5.2 to −15.3%. Both are byte-identical.
- Vectorized banded NW keeps a 3-row rolling score buffer, cutting DP memory
  traffic ~6×: −3.5 to −6% on the DP kernel, ~−3% on `run_dada` (#127).
- Pooled derep load reads JSON in a single pass: −24% on the derep phase and
  −36% peak RSS on a 339 MB PacBio set (#133). The memory half applies to every
  tagged-JSON reader.

### CI
- The R concordance gate is now two-level: exact (every sample × ASV count
  equals R's) for Illumina NW, both pseudo arms and PacBio k=5 NW; near (same
  ASV set, ≤ max(10 reads, 0.1%) moved) for PacBio k=7 and WFA. Tests run on a
  `test-opt` profile (release with debug assertions).

## [0.2.0] - 2026-07-06

### Added
- `summary` subcommand: per-position quality metrics, cumulative expected-error
  metrics, and an optional per-read sequence-complexity histogram (#8).
- `sample` and `errors-from-sample` subcommands for subsampling input FASTQs and
  learning an error model from the subsample (usable for bootstrapping).
- `--failed-uniques` diagnostic on `dada` / `dada-pseudo` / `dada-pooled`:
  emits a TSV of uniques that failed to seed a partition (#60).
- `chimera-diagnostics` for higher-order (trimera) screening, with
  nearest-parent and flank gates (#54).
- Experimental **WFA alignment backend** (`--align-backend wfa2`) with a
  `--wfa-max-edits` edit-budget cap (#49, #51). ASV-equivalent to Needleman-
  Wunsch on tested Illumina and PacBio HiFi data, but not byte-identical.
- `merge-pairs` records input-file provenance and warns on mismatch (#10).
- `just` / `make` task runners for build, install, test, and docs (#46).

### Changed
- The experimental WFA backend is now gated behind an **off-by-default `wfa`
  Cargo feature** (#63). Default builds — and the published crate — are
  Needleman-Wunsch only; selecting `--align-backend wfa2` without the feature
  errors rather than silently falling back to NW. Building WFA requires a source
  checkout (`cargo build --features wfa`) because it depends on a git crate that
  cannot ship on crates.io.
- Denoising/error-learning verbose logs now echo the active alignment backend
  (#51).
- Trimmed the published crate: development, benchmarking, concordance data
  (multi-MB PacBio FASTQs), examples, notes, and CI/infra files are excluded.

### Performance
- Sparse k-mer-8 screen, gated to k ≥ 8, cutting resident memory in the
  high-k regime (#43).
- `dada-pooled` streams dereplication into merging instead of holding all
  dereps in memory, and uses an integer merge-quality accumulator (#39, #41).
- Various k-mer-vector memory reductions in the alignment screen (#32).

### CI
- R DADA2 concordance guardrail for Illumina and PacBio, with abundance/recall/
  precision gates (#35).

## [0.1.0] - 2026-06-03

- Initial crates.io release: Rust ports of the core DADA2 pipeline steps
  (`filter-and-trim`, `derep`, `learn-errors`, `dada`, `merge-pairs`,
  `remove-bimera-denovo`, `assign-taxonomy` / `assign-species`, sequence-table
  helpers).

[Unreleased]: https://github.com/HPCBio/dada2-rs/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/HPCBio/dada2-rs/compare/v0.2.0...v0.5.0
[0.2.0]: https://github.com/HPCBio/dada2-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HPCBio/dada2-rs/releases/tag/v0.1.0
