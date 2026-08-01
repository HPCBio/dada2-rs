# Reading the prep before the result

**A run can be internally consistent and still be entirely wrong.** This page
records a case where public metadata named the wrong primer pair, the library
retained both primers behind variable-length spacers, and the resulting analysis
produced a confident, reproducible, **completely misleading** result — three
times in a row, each retraction moving in the same direction.

Nothing in any artifact flagged it. The spurious ASVs were real sequences at real
abundances, the error models fit, the pipeline reported no warnings. It surfaced
only because someone asked what the merged read lengths looked like.

The finding it corrupted is the [NovaSeq 6000 errfun
A/B](binned-quality-illumina-novaseq.md); this page is the methodological
residue, and it applies to **any dataset whose prep you did not perform
yourself** — SRA/ENA deposits above all.

## What happened

The dataset is a public NovaSeq 6000 soil 16S run. Its BioProject metadata
states:

```
Design: 515F-806R (V4)
```

That is wrong. The data are **V3–V4 (341F/805R)**. Direct primer search across
merged reads:

| motif | found in | position |
|---|---|---|
| **341F** `CCTACGGGNGGCWGCAG` | **99.96%** | 5' end, offsets 0–4 |
| **805R** (reverse complement) `GGATTAGATACCCBDGTAGTC` | **99.99%** | 1 nt from the 3' end |
| 515F `GTGYCAGCMGCCGCGGTAA` | 98.08% | **position ~152 — internal** |

If these were V4, 515F would sit at position 0 and 341F would be absent. Instead
515F is an internal site of the *longer* V3–V4 product, exactly as expected.

Two further properties, neither documented:

1. **Both primers were retained** — no primer removal was performed at any point.
2. **Both are IUPAC-degenerate** (341F has `N` and `W`; 805R has `H` and `V`),
   and each sits behind a **heterogeneity spacer** — Fadrosh-style 0–4 nt padding
   used to raise base diversity on patterned flow cells.

## Why this is worse than it sounds

A degenerate primer behind a variable spacer is a **synthetic variation
generator**. The same biological molecule enters denoising as several distinct
sequences, differing at the spacer offset and at the degenerate positions. Those
differences are genuine sequence differences, so `dada` is right to keep them
apart — it is being fed a lie, not making an error.

With this happening at **both ends of every read**, the two effects multiply.
Collapsing the final ASVs to the **insert only** — the span strictly between 341F
and 805R-revcomp — removed **65.7%** of one arm's ASVs and **74.5%** of the
other's as duplicates of sequences already present.

Roughly **two-thirds to three-quarters of the ASV table was an artifact of
untrimmed prep.**

## The symptoms all pointed somewhere else

This is the part worth internalising. Every downstream symptom had a plausible,
wrong explanation ready to hand:

| symptom | the tempting diagnosis | the actual cause |
|---|---|---|
| ~26% ASV difference between error models | the error model matters enormously here | primer/spacer variation, amplified differently by each model |
| 44% merge rate | poor libraries, primer dimer, over-cycling | a 443 bp product on 2×250 leaves only ~37 nt of overlap |
| `loess` failed to converge (3 of 4 runs) | binned quality breaks the LOESS fit | the fit was chasing artifactual variation |
| ASV counts in the tens of thousands | soil is diverse | inflated 3–4× |

Each of these is a *real observation*. Each supports a conclusion that is false.
The merge rate in particular reads as a library-quality problem and was initially
written up as one; the length distribution disconfirmed it (unimodal 443–449 bp,
**0.00%** of merges below 200 bp, **100.0%** of accepted merges carrying zero
mismatches in the overlap). Primer dimer and mis-amplification produce short
fragments and a multimodal profile. Neither was present.

## After correction

The reads were re-trimmed with `cutadapt` against the actual V3–V4 primers
(`dev/cutadapt_v3v4.sh`) and the entire analysis was re-run from `learn-errors`
forward. Every symptom above resolved:

| | untrimmed | primers removed |
|---|---|---|
| Merge rate | 44.1% | **83.8%** |
| `learn-errors` convergence | 1 of 4 runs | **4 of 4** |
| Post-chimera ASVs (`loess`) | 26,097 | **22,498** |
| ASV Δ between errfun arms | −26.1% | **+3.4%** — *and the sign flipped* |

The sign flip is the sharpest indictment. On untrimmed data `loess` produced
*fewer* ASVs despite fitting *lower* error rates — a contradiction recorded at
the time as an unexplained anomaly. On trimmed data it produces *more*, which is
what the model predicts. **The artifact was strong enough to reverse the
direction of the effect**, not merely inflate its magnitude.

## What this dictates

**1. Verify the amplicon from the reads, never from the metadata.** One command,
before anything else. Search for the expected forward primer and the **reverse
complement** of the reverse primer, and report both the hit rate and the
distribution of start offsets:

```bash
zcat R1.fastq.gz | awk 'NR%4==2 {
    n++
    if (match($0, /CCTACGGG.GGC[AT]GCAG/) && RSTART <= 9) { f++; off[RSTART-1]++ }
} END {
    printf "fwd primer: %d/%d (%.1f%%)\n", f, n, 100*f/n
    for (i=0; i<10; i++) if (off[i]) printf "  offset %d: %.1f%%\n", i, 100*off[i]/f
}'
```

Three outcomes, three different responses:

- **Not found** → the metadata names the wrong primers. Find the real ones before
  proceeding.
- **Found at a single offset** → fixed-offset trimming is safe.
- **Found across several offsets** → heterogeneity spacers. Fixed-offset trimming
  **cannot** work; see below.

!!! danger "`filter-and-trim` cannot remove a variable-offset primer"
    `filter-and-trim` removes leading bases by **fixed offset**
    (`--trim-left`); it does not match primer sequence. With the primer at five
    different offsets there is no correct value — too small leaves spacer bases,
    too large eats real sequence, and either way one biological sequence becomes
    several ASVs. Use `cutadapt` (or R's `removePrimers()`) upstream, and reserve
    `filter-and-trim` for quality filtering. Tracked in [#113].

**2. Search the reverse complement at the 3' end.** A merged read contains the
reverse primer *reverse-complemented*. Searching for it in forward orientation
returns 0.00% and reads as "the reverse primer was already trimmed" — a wrong
conclusion that, in this investigation, produced an incomplete correction and a
retraction that had to be issued twice. The forward-orientation hit rate was
0.00%; the reverse complement was present in **99.99%** of reads.

**3. Screen read overlap against the expected amplicon size, up front.** A
per-sample fragment-length/overlap heatmap run *before* denoising surfaces all of
this at once: wrong primer pair, retained spacers, thin overlap budget, and
genuine background such as primer dimer. Everything on this page was recovered
retrospectively from merge output, which is both harder and later than it needs
to be. This applies to every paired run regardless of chemistry.

**4. Distrust a confident result on data you did not prepare.** The corrupted
analysis was self-consistent at every internal check. What broke it open was an
external expectation — *"these should be ~240 nt"* — colliding with an observed
median of 443. Absolute magnitudes on unfamiliar data need an outside anchor
before they mean anything.

**5. A within-dataset A/B survives what absolute counts do not.** Both arms here
carried identical artifacts, so the *comparison* remained structurally valid even
while every *number* in it was wrong. That is a real property worth relying on —
but note it did not save the direction of the effect, which flipped. A/B design
buys less protection than it appears to.

## Scope

This is a methods page, not a result. It carries no claim about error models,
binned quality, or this dataset's biology — for those see the [NovaSeq 6000
finding](binned-quality-illumina-novaseq.md), which is the corrected analysis
that emerged once the prep was understood.

[#113]: https://github.com/HPCBio/dada2-rs/issues/113
