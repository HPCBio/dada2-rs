# `dada`

Denoise one or more samples independently — the equivalent of R DADA2's
`dada(..., pool = FALSE)`. For pooled and pseudo-pooled inference see
`dada-pooled` and `dada-pseudo`.

```bash
# single sample -> stdout or -o
dada2-rs dada sample.derep.json.gz --error-model err.json -o sample.dada.json

# many samples, processed independently and serially
dada2-rs dada derep/*.json.gz --error-model err.json --output-dir dada/ --threads 24
```

Each input may be a FASTQ file (dereplicated in memory) or a JSON file produced
by `derep` or `sample` (`.json` / `.json.gz`). Pre-dereplicated input avoids
re-reading the FASTQ when iterating on parameters. Output is a JSON object
describing the inferred ASVs.

With a single input the result goes to `--output` / `-o` (or stdout). With more
than one input the samples are processed independently and serially — **not**
pooled — and one `{sample}.json` per sample is written to `--output-dir`, which
is then required.

## Input

**`--sample-name`** — sample identifier written to the output JSON. Defaults to
the filename stem of the input.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

**`--prior`** — FASTA file of prior sequences (uncompressed or gzipped). Each
sequence in the file that matches a dereplicated unique exactly is flagged as a
prior, making it immune to the abundance p-value filter. Prior-based splitting
uses `--omega-p` instead of `--omega-a`.

## Error model

**`--error-model`** (required) — JSON error model produced by `learn-errors`.

**`--use-err-in`** — use the model's `err_in` matrix instead of `err_out`.
`err_out` (the rates estimated from the final transition counts) is the default
and is what R DADA2 uses downstream.

**`--inherit-err-params`** — inherit any unspecified algorithm parameters
(`omega_*`, `min_*`, `detect_singletons`, `band`, `homo_gap_p`, `kdist_cutoff`,
`kmer_size`, `no_kmer_screen`) from the error model JSON's `params` block. Any
flag passed explicitly on the CLI still wins.

Without this flag the built-in CLI defaults apply, and a warning is emitted for
each CLI value that disagrees with the error model's value. That warning is the
point: denoising with parameters that differ from the ones the model was fit
under is usually a mistake, and this makes it visible rather than silent.

Note that `OMEGA_C` is **not** inherited — see the
[parameters page](../parameters.md).

## Denoising

These are R's `setDadaOpt()` parameters and all default to the R values; the
[parameters page](../parameters.md) has the full equivalency table. Omitting a
flag means "inherit from the error model if `--inherit-err-params` is set,
otherwise use the R default".

| Flag | R option | Default |
|---|---|---|
| `--omega-a` | `OMEGA_A` | 1e-40 |
| `--omega-c` | `OMEGA_C` | 1e-40 |
| `--omega-p` | `OMEGA_P` | 1e-4 |
| `--min-fold` | `MIN_FOLD` | 1 |
| `--min-hamming` | `MIN_HAMMING` | 1 |
| `--min-abund` | `MIN_ABUNDANCE` | 1 |
| `--detect-singletons` | `DETECT_SINGLETONS` | false |
| `--max-clust` | `MAX_CLUST` | 0 (unlimited) |
| `--greedy` | `GREEDY` | true |
| `--use-quals` | `USE_QUALS` | true |

`--detect-singletons`, `--greedy`, `--use-quals` and `--no-kmer-screen` are
tri-state: omit to inherit or use the default, or pass `true` / `false` to set
them explicitly.

## Alignment

**`--band`** — alignment band radius, R's `BAND_SIZE`. 16 is the Illumina
default; 32 is recommended for PacBio HiFi 16S amplicons (per the DADA2 LRAS
manuscript); `-1` is unbanded (O(n²), rarely needed).

The platform split is not cosmetic. On MiSeq data the band can drop to 8 safely,
but 4 breaks; on PacBio HiFi it cannot drop from 32 at all — 32 → 16 already
changes the ASV set. See
[Band size & platform defaults](../findings/band-size-platform-defaults.md).

**`--gap-p`** — gap penalty, R's `GAP_PENALTY` (default −8).

**`--homo-gap-p`** — homopolymer-run gap penalty, R's
`HOMOPOLYMER_GAP_PENALTY`. When unset it falls back to whatever `--gap-p`
resolves to, reproducing R's `NULL` default. PacBio pipelines sometimes set this
closer to 0 (e.g. `-1`) on the theory that homopolymer indels dominate — but
note that the `-1` recommendation originates with 454, not PacBio, and HiFi
data does not have that error mode.

**`--match`** / **`--mismatch`** — R's `MATCH` (+5) and `MISMATCH` (−4). R's
full 4×4 score matrix is commented out upstream; the active path is these
scalars, and dada2-rs matches that path.

**`--align-backend`** — `nw` (default) is Needleman-Wunsch. `wfa2` is the
experimental WFA backend (wfa2lib-rs): ASV-equivalent on the Illumina and PacBio
HiFi data tested, but alignments are not byte-identical. `wfa2` requires a build
with `--features wfa`; a default build — and the published crate — errors if it
is selected.

## Screening

The screen decides which pairs are worth aligning. It does **not** define the
clusters. Set it too aggressively and you lose real comparisons; set it loosely
and you pay in alignment time.

**`--kdist-cutoff`** — k-mer distance cutoff, R's `KDIST_CUTOFF` (default 0.42).
Pairs above this distance are not aligned. Lower screens more aggressively
(faster, more false negatives); raise it for divergent sequences.

Tightening this is not free: it acts through partition and error-model bias, not
just by screening out error copies, and `n_asv` is a misleading metric for
judging it — use set identity. The clean lever is to decouple this cutoff from
the one `learn-errors` runs under, leaving error learning at 0.42 while
denoising tightens to 0.30; see
[KDIST cutoff decoupling](../findings/kdist-cutoff-decoupling.md).

**`--kmer-size`** — k-mer size for the screen and the Raw k-mer vectors, R's
`KMER_SIZE` (default 5, range 3–8; 8 is a hard ceiling because k-mer indices
must fit in a `u16`).

For PacBio HiFi, do **not** leave this at 5: on ~1.4 kb reads the screen becomes
a no-op — nearly every pair is aligned, roughly 4–5× slower at scale. Use k=7
for speed or k=6 to cap memory; both give effectively identical ASVs. Memory
scales as `4^k` per Raw (k=5 → 1 KB, k=6 → 4 KB, k=7 → 16 KB, k=8 → 64 KB).
See [K-mer screen size](../findings/kmer-size-screening.md).

**`--no-kmer-screen`** — disable the screen and align every pair. Much slower;
use only when you suspect the screen is filtering valid comparisons.

## Performance

Neither flag affects results.

**`--threads`** (default 1) — threads for both dereplication and the DADA2
comparison map.

**`--sample-jobs`** — multi-input only: how many samples to denoise
concurrently, each on its own `threads / sample-jobs` sub-pool. A single
sample's comparison map is often too small to feed many threads, so fanning
samples across smaller sub-pools keeps every core fed — roughly 4 threads per
sample is the sweet spot — and bounds memory to this many samples in flight.
Defaults to `round(threads / 4)`, i.e. serial at ≤4 threads. Dial it down for
very large or complex samples if memory is tight.

## Output

**`--output` / `-o`** — single input only; write JSON here instead of stdout.

**`--output-dir`** — required for more than one input; one `{sample}.json` per
sample, the same convention as `dada-pooled`. Created if absent.

**`--compact`** — minified JSON instead of pretty-printed.

**`--gzip`** — gzip the per-sample files as `{sample}.json.gz`.

## Diagnostics

None of these change the ASVs.

**`--verbose`** — progress to stderr.

**`--aux-outputs`** — emit R-DADA2-parity per-cluster diagnostics in the output
JSON: `cluster_stats` (n0/n1/nunq/birth_qave/post-hoc p-value),
`cluster_quality` (mean quality at each reference position), `birth_subs` (the
substitutions that drove each cluster split), and `transitions` (a 16 × nq
transition-by-quality matrix). This costs one extra alignment per Raw against
its cluster center, so enable it only when you need the diagnostics.

**`--cluster-trace`** — write a single full cluster trace (`clusters.json`)
here: cluster centers, members with hamming distance / λ / p-value, and birth
metadata. Useful for ASV-calling QC and one-off plots; see
`examples/cluster_trace/`.

**`--trace-no-members`** — omit the per-cluster `members` array, emitting only
centers and birth metadata.

**`--trace-min-abund`** (default 1) — only include trace members at or above
this abundance.

**`--metrics-json`** — write a structured record of the run's phase times,
compare attribution, shuffle split, footprint and index decision to this file.
Free to leave on: it collects only what the run already pays for.

**`--metrics-attribution`** — add the per-comparison timings (`busy`, map
parallel efficiency, the screen / DP / `al2subs` split). This **slows the run**,
so do not combine it with a timing measurement. Requires `--metrics-json`.

See [Run metrics](metrics-json.md) for the schema and the two levels.

**`--failed-uniques`** — write a TSV of uniques that failed to denoise (final
p-value < `omega_c`, `map == null`). Tidy long format with a header —
`sequence<TAB>sample<TAB>reads` — one row per failed unique per sample it
appears in. This answers "what failed to denoise, and how many reads did it
cost?" without hand-joining `map` back to the derep uniques.

## Experimental

These flags are unstable and off by default. They are not drop-in replacements
for the defaults.

**`--screen-backend`** — an alternative to the k-mer screen. `kmer` (default) is
the ESPRIT-style `4^k` frequency vector that R/C++ DADA2 uses; `minimizer` is a
winnowed-minimizer sketch. Neither defines the clusters.

*Expect mostly-concordant but not identical results.* On every dataset tested
the ASV sets and counts agree closely, but there are real differences: typically
a fraction of a percent of reads and a small number of low-abundance ASVs (all
under ~15 reads in the datasets measured).

It pays off where a lot of screening happens — diverse pools, where most pairs
are dissimilar, so the screen runs on everything and the aligner on almost
nothing. There the k-mer screen can reach 44–77% of denoising time and this
replaces it with under 10%. On low-diversity data the screen is ~1% of runtime
and there is nothing to win.

*It requires tuning.* The default `--kdist-cutoff 0.42` is wrong for this
backend — it over-screens by roughly 3×. Use `--minimizer-k 8` with a cutoff
derived for your data (~0.64 on diverse Illumina pools, ~0.50 on PacBio HiFi)
via `kdist-calibrate --screen-backend minimizer`. See
[Minimizers as the screen](../findings/minimizer-screening.md).

**`--minimizer-k`** (default 8, range 5–31) — k-mer size for the minimizer
sketch. Independent of `--kmer-size`, which still governs the frequency screen.
k=8 was chosen on Illumina via a pair audit (k=11 misses real neighbours) and
has **not** been re-derived per platform the way `--kmer-size` had to be.

**`--minimizer-w`** (default 5, range 1–64) — winnowing window in k-mers. The
sketch retains ~`2/(w+1)` of positions, so larger `w` is smaller and faster but
less sensitive. This was never varied in any published measurement — every
result on the findings page is a single w=5 point.

**`--screen-audit`** — evaluate *both* screens on every comparison and align the
union, reporting how often they disagree and, for each disagreement, how many
substitutions the alignment actually found. Answers whether the minimizer screen
passes a superset, a subset, or a different set of the pairs `--kdist-cutoff`
passes. ASVs are unaffected (audit-only alignments are discarded), but the run
is substantially slower and its timings are not comparable to a normal run.

**`--wfa-max-edits`** (default 50, `0` = unbounded) — WFA edit-budget cap, in
edit operations, used with `--align-backend wfa2`. WFA aborts a pair once it
needs more than this many edits and falls back to the NW path for that pair
(NW-identical there). The budget is an absolute edit count, **not** a fraction
of read length: denoising only aligns near-identical reads (~99.9% identity), so
real error copies stay a few edits apart regardless of read length, while
divergent non-error-copy pairs that slip past the screen stay bounded. Ignored
by the `nw` backend.

Internally an edit budget `E` maps to a WFA cost of `E·|gap_p|` — e.g. 50·8 =
400 with default scoring. The `DADA2RS_WFA_MAX_STEPS` environment override is
specified in those raw cost units, not edits.

## See also

- [Per-sample denoising algorithm](../algorithm-dada.md)
- [Parameters: `setDadaOpt()` parity](../parameters.md)
- [Tuning for your data](../tuning-for-your-data.md)
