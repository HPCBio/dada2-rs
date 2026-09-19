# `kdist-calibrate`

Calibrate the k-mer screen: emit k-mer distance against true alignment
divergence.

For sampled pairs of unique sequences, reports the k-mer distance (the
`KDIST_CUTOFF` screen metric) alongside the true **unbanded** ends-free
alignment divergence. This lets the 0.42 default — nominally ~10% divergence,
calibrated on Illumina 16S — be checked per dataset, platform, `k`, and pooling
regime.

Default output is CSV:
`sample,kdist,edits,core_len,pct_div,screened_in,ab_i,ab_j`. Several modes below
change the columns.

## Picking a pooling regime

The default pools all input uniques into one set — a full-pool geometry
baseline. For the abundance-aware modes, use a biologically meaningful
population instead:

- **`--per-sample`** for the independent (per-sample denoising) regime.
- **`--from-dada-pooled`** for a true pooled run. Pass the `_pooled.json[.gz]`
  record from [`dada-pooled`](dada-pooled.md), **not** the raw derep union, so
  the pool is scored once with pooled abundances and no cross-sample
  double-counting.

```bash
# all-vs-all geometry baseline, per sample (pairs capped by --max-pairs)
dada2-rs kdist-calibrate derep/*.json.gz --k 5 --per-sample --threads 24 -o kdist.allvall.csv

# abundance mode (nearest more-abundant parent) — ALWAYS pair with --per-sample
dada2-rs kdist-calibrate derep/*.json.gz --k 5 --nearest-parent --per-sample --threads 24 -o kdist.abundance.csv

# pooled regime, post-inference: score the merged pool as one population
dada2-rs kdist-calibrate run_pooled.json.gz --k 5 --from-dada-pooled --threads 24 -o kdist.pooled.csv
```

## Input

**`<INPUTS>...`** — derep JSON files, or `dada` output with
`--from-dada[-pooled]`.

**`--derep-dir`** — with `--from-dada`: the directory holding the derep JSONs
that fed `dada`. Matched to each output by sample name: an exact
`{sample}.json[.gz]` first, else a `{sample}.*.json[.gz]` file (e.g. a
pipeline-renamed `{sample}.derep.R1.json.gz`). Ambiguous prefixes are resolved
by the derep JSON's own `sample` field.

## Pooling regime

**`--per-sample`** — compute pairs *within* each sample instead of pooling all
uniques into one set.

**`--nearest-parent`** — abundance-aware mode. Instead of random pairs, link
each unique to its nearest **more-abundant** neighbour — its candidate
error-copy parent — and report the screen's headroom above real error-copy
distances. Output columns become
`sample,ab,parent_ab,ab_ratio,kdist,edits,core_len,pct_div,screened_in`.

!!! warning "`--nearest-parent` is O(n²) and ignores `--max-pairs`"
    It scans every unique against its more-abundant prefix. Pair it with
    `--per-sample`, where n is small per sample and each parent link is a
    genuine within-sample error copy. For a pooled run use
    `--from-dada-pooled`. Do **not** run it over a pooled raw-derep union: that
    is both enormous and cross-sample noise, since a unique's "parent" may live
    in another sample. Bound it with `--max-uniques`.

**`--from-dada`** — post-inference mode: treat the inputs as `dada` output JSONs
rather than derep JSONs, and label every input unique by what denoising actually
did to it — center (survived as an ASV), member (absorbed as an error copy), or
failed (shed by the abundance test). Requires `--derep-dir`. Columns become
`sample,class,cluster,ab,center_ab,ab_ratio,birth_type,birth_pval,kdist,edits,core_len,pct_div,band_req,screened_in`.

**`--from-dada-pooled`** — pooled post-inference mode: treat the inputs as the
`_pooled.json[.gz]` record(s) from `dada-pooled` and screen the merged unique
table against the single global partition. Self-contained, so no `--derep-dir`,
and the pool is assessed as one population with pooled abundances rather than
re-aggregated per-sample splits. Same columns as `--from-dada`, with
`sample = __pooled__`.

**`--derive-cutoff`** — derive-only: report the minimizer cutoff that reproduces
the k-mer screen's **pass rate** on this data, then stop. Skips alignment
entirely, so it runs in seconds where the full curve takes hours — the curve's
cost is aligning every sampled pair unbanded to get true divergence, which the
matched-pass rule never consults.

It targets the k-mer screen's selectivity, which is the safe target, not the
cheapest cutoff that still agrees with it. On PacBio HiFi the ASV table is
identical from 0.45 to 0.60 and this picks 0.50. Sweep if you can afford to.

**`--derive-uniform-pairs`** — with `--derive-cutoff`: sample pairs uniformly at
random instead of abundance-weighted. Uniform is what a calibration *curve*
wants, since it describes the metric — but it is the wrong population for a
*pass rate*, because `b_compare` compares every raw against each cluster
**centre**, and centres are the abundant uniques. On pooled PacBio the
minimizer/k-mer pass ratio is 0.744 on the pairs actually screened and 0.911 on
uniform pairs, so uniform sampling makes the minimizer look 23% less selective
than it is and the derived cutoff overshoots. Kept because the published curves
are uniform.

**`--max-pairs`** (default 200,000) — maximum pairs computed per population,
random-subsampled above this to bound the O(n²) cost. Applies to all-pairs mode
and to the center pairs of `--from-dada[-pooled]`. It has no effect under
`--nearest-parent`, so passing both is rejected.

**`--max-uniques`** (default 0 = keep all) — randomly subsample each sample to at
most this many uniques before pairing. Unlike `--max-pairs` this *does* apply
under `--nearest-parent`, so it is the lever for bounding that mode.

**`--seed`** — RNG seed for reproducible subsampling.

## Screening

**`--k`** (default 5) — k-mer size. R's default is 5; PacBio full-length wants 7.

**`--cutoff`** (default 0.42) — the screen cutoff used for the `screened_in`
flag and the leakage summary.

**`--leak-pct`** (default 5.0) — divergence above which a screened-in pair counts
as *leaked*, i.e. too far apart to be an error copy. This is a crude threshold;
the true ceiling is abundance-dependent.

## Alignment

**`--band`** (default −1) — alignment band radius; negative means unbanded.
Unbanded is the correct default here, because a band truncates the divergence of
distant pairs and so corrupts the very quantity being calibrated.

## Performance

**`--threads`** (default 1) — threads for the parallel alignment.

## Output

**`--output` / `-o`** — write the CSV here instead of stdout.

## Diagnostics

**`--verbose`** — per-population progress and the leakage summary to stderr.

## Experimental

**`--screen-backend`** (default `kmer`) — which screen to calibrate. `kmer` is
the ESPRIT frequency vector; `minimizer` is the winnowed sketch.

**A cutoff does not transfer between them.** On the MiSeq SOP, 0.42 passes 27.6%
of pairs on the frequency vector and 9.0% on the sketch, so the minimizer
backend needs its own curve. See
[Minimizers as the screen](../findings/minimizer-screening.md).

**`--minimizer-k`** / **`--minimizer-w`** — sketch k-mer size and winnowing
window, used with `--screen-backend minimizer`.

## See also

- [KDIST cutoff decoupling](../findings/kdist-cutoff-decoupling.md)
- [K-mer screen size](../findings/kmer-size-screening.md)
