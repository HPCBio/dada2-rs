# Minimizers as the pre-alignment screen

**Verdict so far:** a winnowed-minimizer sketch reproduces the k-mer screen's
decisions on both platforms — **93=93 ASVs on PacBio HiFi and 8=8 on Illumina,
churn=0 at every pipeline stage** — while being a **strict subset** of what the
k-mer screen passes, with **zero disagreements below 10 substitutions**. Its win
is **not speed against a tuned baseline**: at PacBio's production `k=7` it passes
15.57% of pairs against the k-mer screen's 15.78%, a 0.2-point difference. The
win is that it reaches that selectivity with a **74% smaller resident screen**
and **without platform-specific tuning of `k`** — it decouples the two knobs the
frequency vector welds together.

Experimental and opt-in (`--screen-backend minimizer`), on the same footing as
the [WFA alignment backend](https://github.com/HPCBio/dada2-rs/issues/51). The
default is unchanged.

> **Status: early.** Every number below comes from the small concordance
> fixtures (2 samples each). They establish a *mechanism*, not a production
> verdict. The gate before promotion is a deep, diverse run — see
> [What would change the verdict](#what-would-change-the-verdict).

## What the screen is, and what "replacing" it means

The pre-alignment screen decides which candidate pairs are worth aligning. It is
a *speed* optimisation inherited from the ESPRIT lineage — **not** the cluster
definition, which is the partition plus the abundance test downstream. This is
the same framing as [k-mer screen size](kmer-size-screening.md), and it is what
makes a screen swap a legitimate thing to try at all: the screen is allowed to
change its mind about which pairs to align, as long as the pairs it drops were
never going to matter.

Both backends feed the *same* `--kdist-cutoff` gate. The minimizer distance is
deliberately [`kmer_dist8`]'s exact formula

```text
dist = 1 - Σ min(count_a[m], count_b[m]) / min(total_a, total_b)
```

evaluated on a winnowed subsample of the k-mer space rather than on all of it.
That is a choice, not a coincidence: it keeps `KDIST_CUTOFF = 0.42`
interpretable across both backends, so an A/B changes *which k-mers are
consulted* and nothing else about the gate.

## The problem it attacks

The frequency-vector screen welds two knobs together. `k` sets **both**
specificity **and** memory, because the vector is `4^k` entries regardless of
read length. That coupling is the direct cause of two standing problems:

- **`k` must be tuned per platform.** [k-mer screen size](kmer-size-screening.md)
  established that Illumina wants `k=5` and PacBio wants `k=7` *for speed*, and
  that PacBio at `k=5` is "effectively a no-op". This branch measured that claim
  directly for the first time, and it is worse than "effectively": on PacBio HiFi
  at `k=5` the screen passes **100.00%** of comparisons. Not approximately — every
  single one. The default configuration does zero screening on long reads.
- **Memory scales as `4^k`**, which forced the sparse representation in
  [#43](https://github.com/HPCBio/dada2-rs/issues/43) and gates raising the cap
  in [#44](https://github.com/HPCBio/dada2-rs/issues/44).

A winnowed sketch is `O(len / w)` entries *regardless of `k`*. So `k` becomes
free to raise for specificity while memory falls, and one setting can serve both
platforms.

## Design

- **k = 11, w = 5** by default. `k` is chosen purely for discrimination now that
  memory does not depend on it: large enough that unrelated amplicons share
  almost nothing, small enough that a substitution — which destroys the `k`
  k-mers spanning it — removes only ~2–3 of ~48 minimizers in a 250 bp read.
- **Hashed before winnowing** (SplitMix64). On the raw 2-bit packing `A → 00`, so
  a poly-A run is numerically minimal in every window it touches and would be
  selected everywhere, concentrating sketches on exactly the low-complexity
  content that carries no signal.
- **Consecutive duplicate hashes collapse.** A deliberate divergence from the
  mapping-oriented convention (minimap2 et al. key on position, because they need
  to know *where* a seed matched). We only need to know *whether* content is
  shared, and per-occurrence emission actively hurts that: inside a homopolymer
  every window has the same minimal k-mer at a new offset, so the run emits at
  density 1.0 against winnowing's ~2/(w+1). That inflates its count, and
  `Σ min(count)` would then read two sequences as similar merely because both
  contain a homopolymer. This was caught by a unit test asserting the hash
  avalanche was sufficient — it was not, and the fix was in the emission rule,
  not the hash.
- **Fails open, never closed.** When either sketch is absent or empty, the screen
  returns distance 0 ("align it"). `minimizer_dist` would report 1.0 for an empty
  sketch — maximal distance — which would drop the pair on the basis of *no
  evidence at all*. Failing open costs an alignment; failing closed silently
  loses ASVs. Emptiness, not sequence length, is the test: a sequence can be long
  enough for a window and still sketch to nothing if it is dense with non-ACGT.
- **Struct-of-arrays.** `Vec<(u64, u8)>` aligns to 8 and occupies **16** bytes per
  entry, 7 of them padding — 44% of a structure whose entire justification is
  being small. Split into `Vec<u64>` + `Vec<u8>`, an entry costs 9.

## Evidence

### Screen agreement (the mechanism gate)

`--screen-audit` evaluates **both** screens on every comparison, aligns the
**union**, and buckets each disagreement by the substitution count the aligner
actually found. Counts alone cannot distinguish a lost error copy from a
correctly-rejected stranger; the `nsubs` histogram can. Audit-only alignments are
discarded, so the run's ASVs remain the active backend's.

PacBio HiFi, 1 sample, 38,079 comparisons, band 32:

| `--kmer-size` | k-mer passes | minimizer passes | kmer-only | minimizer-only | kmer-only at nsubs ≤ 10 |
|---|---|---|---|---|---|
| 5 (R-parity default) | **100.00%** | 15.57% | 84.44% | **0** | **0** |
| 6 | 19.20% | 15.57% | 3.63% | **0** | **0** |
| 7 (PacBio production) | **15.78%** | 15.57% | 0.21% | **0** | **0** |

Illumina MiSeq, 1 sample, 4,088 comparisons, `k=5` (which *is* the Illumina
production setting):

| | k-mer | minimizer |
|---|---|---|
| passes | 62.62% | **44.10%** |
| kmer-only | 18.52% | — |
| minimizer-only | — | **0** |
| kmer-only at nsubs ≤ 10 | **0** | — |

Two facts hold across every configuration tested:

1. **The minimizer screen is a strict subset.** `minimizer-only = 0` everywhere —
   it never passes a pair the k-mer screen rejects, so it cannot introduce
   spurious comparisons.
2. **Every pair it additionally rejects is far.** All of them sit at >10
   substitutions; the count at nsubs ≤ 10 is exactly zero in all four
   configurations. There is no observed pair in the range where a shrouded
   comparison could cost an ASV.

### ASV concordance

Set identity and abundance, never `n_asv` — per the
[concordance guardrail](../benchmarking.md#5-concordance-validation-tooling) and
the lesson in [k-mer screen size](kmer-size-screening.md).

| Stage | PacBio HiFi | Illumina (paired) |
|---|---|---|
| intermediate `dada` | 93 = 93, **churn 0** | fwd 13 = 13, rev 9 = 9, **churn 0** |
| post-chimera final table | 93 = 93, **churn 0** | 8 = 8, **churn 0** |
| vs static R DADA2 reference | PASS, metrics **identical** to the k-mer arm | PASS |

On a single-sample Illumina run the *partition* is identical too: all 449 uniques
map to the same ASV sequence, with identical abundances. What differs is cluster
creation **order** and the birth diagnostics (`birth_e`, `birth_fold`,
`birth_pval`) — the screen changed which comparisons contributed expected
abundance to which cluster, without changing the outcome.

### Memory

PacBio HiFi, ~1,495 bp reads. Screen structure only (the `kord` vector, 2,982
B/raw, is paid by both backends and excluded):

| screen | bytes/raw | vs k=7 |
|---|---|---|
| k-mer, k=5 (dense) | 1,024 | — *(but a 100%-pass no-op here)* |
| k-mer, k=6 (dense) | 4,096 | −75% |
| k-mer, k=7 (dense) | **16,384** | baseline |
| k-mer, k=8 (sparse, [#43](https://github.com/HPCBio/dada2-rs/issues/43)) | 5,695 | −65% |
| **minimizer, k=11/w=5** | **4,306** | **−74%** |

The honest comparison is against **k=7**, because that is the configuration whose
selectivity the sketch matches (15.57% vs 15.78%). Against `k=5` the sketch costs
4× *more* memory — but `k=5` is the arm that screens nothing on long reads, so
that comparison is not meaningful.

Measured sketch density is **0.333** on PacBio and 0.349 on Illumina, against
winnowing's theoretical `2/(w+1) = 0.333`. The sketch is behaving exactly as the
theory predicts, which is a useful sanity signal that nothing is degenerate.

### The inverted index

`b_compare` scans all raws against one cluster center, so the raws are the stable
side and the center is the query. `MinimizerIndex` builds posting lists over the
raws once and derives every raw's shared-minimizer count in a single scatter,
replacing `nraw` merge-joins with an array read.

It is an **exact** acceleration, and that is deliberate: it changes cost, never
output, so index-on and index-off must agree byte-for-byte. **Verified** — the
two arms produce identical JSON on real data. Keeping the metric change and the
acceleration on separate axes is what lets the ASV results above be attributed to
the metric alone. Disable for A/B with `DADA2RS_MINIMIZER_INDEX=0`.

## What this does *not* show

Stated plainly, because the first framing of this result was wrong:

- **It is not a 6.4× speedup.** That figure compares against PacBio `k=5`, a
  baseline we already knew was a no-op. Against the tuned `k=7` the alignment
  count is a wash (0.2 points).
- **No wall-clock number is reported here at all.** The fixtures are seconds
  long, and per [A/B measurement hygiene](measuring-on-numa.md) a timing claim
  needs a pinned NUMA node, replicated arms, and an untouched control channel.
  None of that has been done.
- **The fixtures are small** — 8 and 93 ASVs. `churn=0` on 8 ASVs is close to
  vacuous. The PacBio 93=93 with churn 0 at *both* stages is the stronger of the
  two, and still only 2 samples.
- **Error-model provenance is not wired.** `screen_backend` is deliberately not
  resolved through the error model's `params` block, because released models have
  no such field and inheritance would silently override an explicit
  `--screen-backend`. Consequence: a model learned under one screen and applied
  under the other is **not** flagged the way a `kdist_cutoff` mismatch is.

## Incidental finding: ASV ordering is not deterministic

While checking byte-identity, two runs of the **same binary on the default k-mer
backend** produced the same ASV set with the same counts in a **different
order**. This is pre-existing and unrelated to this branch, but it matters
methodologically: **byte-identity checks on `seqtab.nochim.json` are unreliable**,
and any comparison must be set- and count-based. `dev/compare_asvs.py` already
is; ad-hoc `diff`/JSON-equality checks are not. Worth its own issue.

## What would change the verdict

In rough order of decisiveness:

1. **A deep, diverse run.** The MiSeq SOP and the 95-sample PacBio set, through to
   the chimera-filtered table. Small fixtures agreeing proves the screen is not
   grossly wrong; it does not prove the tails agree, and the tails are where a
   screen change shows up.
2. **A wall-clock A/B done properly** — pinned, replicated, with a control
   channel. The prize is bounded: the screen is 15–23% of `b_compare`
   ([screen vs align](compare-screen-vs-align.md)), and against `k=7` the
   minimizer screen does not reduce alignments. The *plausible* win is on runs
   left at the `k=5` default, where it converts a 100%-pass no-op into a real
   screen without asking the user to know about `--kmer-size`.
3. **`k`/`w` sensitivity.** Every number here is one point (k=11, w=5). Whether
   the "zero disagreements below 10 substitutions" result is robust across the
   parameter space, or a lucky point, is untested.
4. **Low-complexity and short-amplicon data.** The homopolymer emission rule and
   the fail-open path are both reasoned-from-first-principles and unit-tested,
   but neither has met real ITS2 or a short amplicon where a meaningful fraction
   of reads falls below one window.

## Reproducing

```bash
# End-to-end A/B through chimera removal (both harnesses take SCREEN_BACKEND).
bash dev/concordance/run_illumina.sh ./target/release/dada2-rs \
    dev/concordance/data/illumina /tmp/ab_kmer 4
SCREEN_BACKEND=minimizer bash dev/concordance/run_illumina.sh \
    ./target/release/dada2-rs dev/concordance/data/illumina /tmp/ab_mini 4
python3 dev/compare_asvs.py --baseline kmer=/tmp/ab_kmer/seqtab.nochim.json \
    --compare minimizer=/tmp/ab_mini/seqtab.nochim.json

# Pair-level screen agreement, the mechanism gate.
./target/release/dada2-rs dada <filtered.fastq.gz> --error-model <err.json> \
    --screen-backend minimizer --screen-audit --kmer-size 7 -o /dev/null

# Exactness of the index (must be byte-identical).
DADA2RS_MINIMIZER_INDEX=0 ./target/release/dada2-rs dada ... --screen-backend minimizer
```

[`kmer_dist8`]: https://github.com/HPCBio/dada2-rs/blob/main/src/kmers.rs
