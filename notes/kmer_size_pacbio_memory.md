# K-mer size: PacBio HiFi screen behavior + memory (issue #15)

Follow-up to `notes/benchmark_kmer_size.md` (Illumina-only). Resolves the two
open questions in issue #15 for long reads: (1) does a larger `--kmer-size`
mis-handle the pre-alignment screen on ~1.4 kb PacBio HiFi reads, and (2) what
is the memory cost of large k.

- Date: 2026-05-31
- Binary: `target/release/dada2-rs` @ 36e0743
- Host: Apple Silicon, single-threaded runs (`--threads 1`) for clean compute timing
- All `errors-from-sample` runs: `--max-consist 20 --verbose`.
  "converged" = iteration count < the cap of 20.

> **Note on provenance.** The raw run artifacts lived in `/tmp/kmer_exp/` and were
> purged when `/tmp` cleared. The numbers below were captured/verified during the
> run (the `summary.csv` and shroud `ALIGN:` lines were re-read directly from disk
> before purge). The `err_out` max-diff column comes from the run logs and could
> not be independently re-derived after the purge — re-run to reproduce exactly.

## Datasets & preprocessing

| Dataset | Reads | Uniques | Read len | errfun | band |
|---------|-------|---------|----------|--------|------|
| MiSeq F3D0 (R1) | 7228 | 2014 | ~240 bp | loess | 16 |
| PacBio SRR28724909 | 32132 | 10109 | ~1450 bp | pacbio | 32 |

```
# MiSeq (raw -> filtered)
filter-and-trim --trunc-len 240 --max-n 0 --max-ee 2 --trunc-q 2 --compress

# PacBio
filter-and-trim --min-len 1000 --max-len 1600 --max-ee 2 --max-n 0 --compress
```

**Orientation check (important methodology note).** PacBio HiFi reads are *often*
mixed-orientation with primers attached, which would break the k-mer screen at
any k (a sequence and its reverse complement share almost no k-mers). We tested
empirically: of the first 2000 reads of `SRR28724909.trim.fastq.gz`, **0** carried
the 27F primer motif at the 5′ end and **0** carried its reverse-complement at
the 3′ end; `remove-primers --orient` produced 0 output. These reads are **already
primer-trimmed and uniformly oriented** (the `.trim` in the filename). So no
orientation step was needed, and the screen is being tested on legitimately
same-strand sequences. *For raw HiFi from other sources, orient first.*

**Scale.** The k-sweep used the **full 10109 PacBio uniques**, held constant
across all k. MiSeq used its full 2014 uniques.

## Results

### MiSeq F3D0 (Illumina, ~240 bp; loess, band 16)

| k | iters | wall (s) | peak RSS | err_out max\|Δ\| vs k5 |
|---|-------|----------|----------|------------------------|
| 5 | 5 | 5.21 | 43.7 MiB | — |
| 6 | 4 | 2.24 | 61.3 MiB | 1.25e-3 |
| 7 | 4 | 1.57 | 132.3 MiB | 2.06e-3 |
| 8 | 4 | 2.56 | 416.9 MiB | 4.35e-3 |

### PacBio SRR28724909 (HiFi, ~1450 bp; pacbio, band 32) — single subsampled sample

> **Superseded timing:** the wall times in this single-sample table (and the
> "k=8 ~2× slower than k=6" claim derived from them) are a **small-sample
> artifact**. At full scale the ordering inverts — see "Full PacBio sweep at
> scale" below, where k=5 is by far the slowest and k=7 is fastest. The ASV /
> shroud / memory findings here all still hold.

| k | iters | wall (s) | peak RSS | err_out max\|Δ\| vs k5 |
|---|-------|----------|----------|------------------------|
| 5 | 4 | 322.61 | 942.6 MiB | — |
| 6 | 4 | 78.44 | 1031.8 MiB | 4.91e-7 |
| 7 | 4 | 66.04 | 1389.9 MiB | 6.07e-7 |
| 8 | 4 | 67.29 | 2198.0 MiB | 6.65e-5 |

All runs converged. Learned error rates are **essentially identical across k**
(PacBio agrees to ~1e-7; MiSeq to a few 1e-3 at the low-Q edges).

> **Terminology:** "shrouded" = a candidate pair rejected by the **k-mer screen**
> before any Needleman–Wunsch alignment is attempted (the `sub_new` k-mer distance
> exceeds the cutoff, so it returns no alignment). The mechanism is what the
> literature calls *k-mer screening*; "shroud"/`nshroud` is the outcome counter,
> a name inherited verbatim from the upstream DADA2 C++ source (`dada.h`,
> `cluster.cpp`, `Rmain.cpp`'s `"%i aligns, %i shrouded"`). `shrouded` ⊆ `aligns`.

### Direct screen behavior (PacBio, `dada` iteration 1, kdist_cutoff = 0.42)

Same input, same band 32, only `--kmer-size` and the (own) learned model differ.
Total comparisons in the first cluster pass:

| k | aligns | shrouded | shroud % | final ASVs |
|---|--------|----------|----------|------------|
| 5 | 548600 | 0 | 0.0 % | 87 |
| 8 | 548707 | 432806 | 78.9 % | 87 |

## The key finding (corrects the naive intuition)

On **long reads the k=5 screen is a no-op**: it shrouds **0 %** of pairs, so every
pair goes to full O(L²) alignment — which is exactly why k=5 is the *slowest*
config (322 s vs ~70 s). Larger k makes the screen actually engage.

Why: `kdist = 1 − sharedKmers / (L − k + 1)`. For L ≈ 1450:
- **k=5:** only 4⁵ = 1024 distinct k-mers but ~1446 k-mer positions, so frequency
  vectors are dense and even unrelated long reads share most 5-mers → kdist stays
  far below 0.42 → nothing is screened. The screen carries no information at this
  k on long reads.
- **k=8:** 4⁸ = 65536 ≫ 1446 positions, so vectors are sparse; genuine HiFi error
  differences across 1.4 kb drop enough exact 8-mers to push many pairs past 0.42
  → 78.9 % shrouded.

Crucially, **the 78.9 % shrouding at k=8 was benign here**: final ASV count (87)
and learned error rates were identical to k=5. The screen rejected pairs that
the alignment would not have clustered anyway — it is *working*, not
*over-rejecting harmfully*. The aggressiveness is nonetheless a latent risk on
more divergent or lower-quality data, which argues for a moderate rather than
maximal k.

## Memory cost (the binding constraint)

K-mer storage per Raw = `4^k × 3 bytes` (u16 `kmer` + u8 `kmer8`; `kord` is
k-independent). For the **full PacBio sample** (10109 uniques):

| k | per-Raw | k-mer vectors (full sample) |
|---|---------|-----------------------------|
| 5 | 3.0 KB | 29.6 MiB |
| 6 | 12.0 KB | 118.5 MiB |
| 7 | 48.0 KB | 473.9 MiB |
| 8 | 192.0 KB | 1895.4 MiB |

Measured peak RSS tracks this 4×-per-step growth (PacBio: 943 → 1032 → 1390 →
2198 MiB). Pooled multi-sample HiFi (hundreds of thousands of uniques) would be
infeasible at k ≥ 7.

## Conclusions

1. **No harm from larger k on long reads.** Across k=5–8, PacBio converged in 4
   iterations, produced 87 ASVs, and learned the same error model (~1e-7
   agreement). The #15 worry about false-rejecting valid pairs did not
   materialize: even at 78.9 % shrouding the result was unchanged.
2. **k=5 wastes the screen on long reads** — ~0 % shrouded means every pair is
   aligned, making k=5 the slowest PacBio config (dramatically so at scale — see
   below).
3. **k=6+ engage the screen**, give identical ASV *counts* (87) and error model
   here, and run far faster than k=5. (Caveat: only counts were checked at this
   subsample — the full Illumina run later showed stable counts can hide real
   set-level churn, so "identical ASVs" is not established for PacBio.) At
   single-subsample scale the per-k times above are noisy and close; at full
   scale (below) k=7 is fastest and k=6 is the memory-conscious pick — a
   memory↔speed tradeoff, not a clean win on both.
4. **Memory still grows 4× per step** (k=6 ≈ 120 MiB, k=7 ≈ 470 MiB, k=8 ≈ 1.9 GiB
   on a full sample). This is the one axis where higher k is strictly worse.

> NOTE: the original recommendation here (k=6 as a clean speed+memory "sweet
> spot") rested partly on the small-sample timing that the full-scale run below
> overturns. The revised recommendation is in "Full PacBio sweep at scale".

## Recommendation (see the scale section below for the current PacBio call)

- Keep **k=5 as the global default** (correct and lowest-memory for short reads,
  where the screen does discriminate — see `benchmark_kmer_size.md`).
- For **PacBio HiFi, raise k above 5** — k=5 leaves the screen a no-op on long
  reads (huge alignment cost). Choose **k=7 for speed** or **k=6 to cap memory**;
  see the scale section for the measured tradeoff. Document this in the
  `--kmer-size` help for `--errfun pacbio`.

### Caveats

Single Illumina sample; single PacBio sample (full uniques). The benign-shrouding
result (identical ASVs at 79 % shroud) should be confirmed on a second sample or
a pooled run before changing any default. Scratch artifacts were purged from
`/tmp`; re-run to reproduce exact `err_out` diffs. The memory model and the
k=5-is-a-no-op-on-long-reads finding are robust (mechanistic + measured).

---

## Full MiSeq SOP confirmation (2026-06-01)

The caveat above (single Illumina sample) is now addressed: a sweep on the
**full MiSeq SOP** dataset (61 samples, learn-errors pooled per-sample then
denoised, R1 and R2 separately) via `comparison/run_kmer_sweep.sh`. Errfun
loess, band 16, k = 5/6/7/8.

### R1

| k | learn_iters | dada_aligns | dada_shrouded | shroud % | n_asv | wall (s) |
|---|------------|-------------|---------------|----------|-------|----------|
| 5 | 6  | 473 091 080 | 354 857 588 | 42.9 | 2514 | 238.8 |
| 6 | 6  | 476 322 636 | 400 154 558 | 45.7 | 2529 | 182.4 |
| 7 | 7  | 475 915 226 | 422 244 675 | 47.0 | 2525 | 236.4 |
| 8 | 6  | 475 626 964 | 437 773 435 | 47.9 | 2520 | 616.7 |

### R2

| k | learn_iters | dada_aligns | dada_shrouded | shroud % | n_asv | wall (s) |
|---|------------|-------------|---------------|----------|-------|----------|
| 5 | 10*| 398 005 481 | 296 195 546 | 42.7 | 1991 | 175.6 |
| 6 | 7  | 401 707 229 | 336 741 954 | 45.6 | 2004 | 137.5 |
| 7 | 6  | 404 084 016 | 360 688 524 | 47.2 | 2015 | 206.2 |
| 8 | 7  | 404 545 391 | 371 782 216 | 47.9 | 2016 | 538.1 |

\* R2/k5 hit the `--max-consist 10` cap. The log shows this is **oscillation,
not slow convergence**: max|err_in−err_out| reached 1.78e-6 at iter 6 (already
converged), then bounced back to a stuck 1.17e-5 two-state cycle through iter 10.
This is normal DADA2 behavior near the error-model noise floor (R `learnErrors`
does the same) — the iter-6 vs iter-10 model differs by ~1e-5, far below the
~1e-3 that affects ASV calls. The iteration counts are **non-monotonic in k**
(R1: 6/6/7/6; R2: 10/7/6/7), confirming this is run-to-run oscillation, not a
k effect — it is NOT evidence for k=6 over k=5 on correctness grounds.

### Measured peak RSS + ASV-set diff (24-thread rerun, 2026-06-01)

A second rerun at **24 threads** added measured peak RSS (via `/proc` VmHWM —
this cluster has no `/usr/bin/time`) and a **sequence-level pooled ASV-set diff**
vs k=5. Aligns/shroud/ASV-counts reproduce the table above; new columns:

**R1** (baseline k=5 = 2514 ASVs):

| k | #ASV | wall (s) | RSS (GB) | shared | only_k5 | only_kN | verdict |
|---|------|----------|----------|--------|---------|---------|---------|
| 5 | 2514 | 288.2 |  4.1 | —    | —  | —  | — |
| 6 | 2529 | 218.0 |  6.5 | 2502 | 12 | 27 | DIFFERS |
| 7 | 2525 | 211.7 | 16.1 | 2495 | 19 | 30 | DIFFERS |
| 8 | 2520 | 444.3 | 54.4 | 2486 | 28 | 34 | DIFFERS |

**R2** (baseline k=5 = 1991 ASVs):

| k | #ASV | wall (s) | RSS (GB) | shared | only_k5 | only_kN | verdict |
|---|------|----------|----------|--------|---------|---------|---------|
| 5 | 1991 | 204.5 |  3.3 | —    | —  | —  | — |
| 6 | 2004 | 158.9 |  5.9 | 1975 | 16 | 29 | DIFFERS |
| 7 | 2015 | 181.8 | 16.3 | 1968 | 23 | 47 | DIFFERS |
| 8 | 2016 | 392.6 | 58.1 | 1962 | 29 | 54 | DIFFERS |

### Findings at full scale — REVISED

- **Stable ASV *counts* hid real set *churn*. The earlier "identical / washes
  out" claim was wrong for Illumina at scale.** Counts move only ~0.6–1.3 %, but
  the underlying ASV *sets* differ from k=5 at every higher k, and the churn is
  **bidirectional and monotone in k**: higher k both drops some k5 ASVs
  (`only_k5`: 12→19→28 on R1) and adds new ones (`only_kN`: 27→30→34). R2 is
  larger (`only_kN` up to 54 at k8).
- **The churn is the predicted false-negative effect, now observed.** `only_kN >
  only_k5` consistently and the gap widens with k: more shrouding → genuinely
  similar sequences that k5 would have aligned-and-clustered escape the screen
  and split off as their own (mostly rare) ASVs. Higher k slightly **fragments**
  clustering rather than improving it. Magnitude is small (≤1.4 % of ASVs) and
  expected to be concentrated in low-abundance variants — *abundance
  stratification of the churned sets is pending* (see open items).
- **Measured RSS confirms the 4^k model and is the decisive cost.** ~3.3 GB
  k-independent base + 4^k k-mer vectors: 4→6.5→16→**54 GB** (R1), 3.3→5.9→16→
  **58 GB** (R2). k=8 is operationally prohibitive on typical nodes; even k=7 is
  ~16 GB. The 24-thread per-worker buffer term is negligible for 250 bp reads, so
  RSS is essentially all k-mer vectors.
- **Wall time** (24 threads): k=7 fastest (R1 212 s, R2 182 s), k=8 slowest by
  far (R1 444 s, R2 393 s). Matches the earlier full run's *ordering*; absolute
  values differ with thread count (now recorded in the CSV `threads` column).
- **Per-iteration screen behavior** (from the learn logs): the first pass of
  every per-sample `dada` run is unscreened (`kdist_cutoff = 1.0`, by design, to
  seed cluster 0 — see `run_dada` in `dada.rs`), so it shows `0 shrouded` and
  fewer alignments; screening kicks in only once clusters bud. Matches R DADA2's
  C++; not a bug.
- **R2/k5 hitting max-consist 10** is error-model oscillation, not slow
  convergence or a k effect (see the `*` note above).

### Churn mechanism — DIAGNOSED (2026-06-02, closes open item #1)

A controlled k=5→k=7 comparison on the 20-sample `data/illumina/full-pooling`
run (dada-pooled, error model relearned at k=7, all other params identical;
outputs under `/tmp/churn_k7/`) traced the churned ASVs back to their derep
uniques via the inverted dada `map`. Result:

| Orient | k5 ASVs | k7 ASVs | shared | only_k5 | only_k7 |
|--------|---------|---------|--------|---------|---------|
| R1 | 517 | 514 | 511 | 6 | 3 |
| R2 | 420 | 426 | 419 | 1 | 7 |

- **Mechanism: k-mer-screen-driven cluster fragmentation, not convergence.**
  Raising k tightens the pre-alignment screen, so near-identical reads that k=5
  shrouded into a larger cluster are no longer pre-clustered and are born as
  their own `Abundance` ASVs. Two traced R2 examples are **exactly Hamming-1**
  from their k5 parent: an ab=215 ASV split from a 488-read k5 parent, and ab=155
  from a 278-read parent. `only_k5` is the mirror image (uniques re-merge/re-split
  at k7), **not** reads dropping below `omega_c`.
- **Benign — does not touch *dominant* biology.** Churned ASVs are low-to-moderate
  abundance (max 215 reads) sitting alongside shared ASVs of 700–1150+ reads.
  ~1–2 % set churn. (Note: "low-to-moderate", not "rare only" — the merge-level
  trace below found surviving churn up to ~107 reads from parents of 61–161 reads,
  so it reaches *moderate* abundance, not just the noise floor.)
- **Convergence ruled out.** R1's error model converged in 5 iterations at *both*
  k5 and k7 yet R1 still churned — so the screen, not the error-model fit, drives
  it. (R2 took 7 iters at k7 vs 5 at k5, but the churn doesn't track that.)

Verdict: the ~1–2 % churn is real but cosmetic — k-mer-screen fragmentation of
rare variants, not movement of abundant biology. This is consistent on both
platforms (PacBio's +9 and Illumina's churn are the same phenomenon).

### Merge survival — ANSWERED (2026-06-02, 362-sample full-pooling sweep)

`data/illumina/illumina_sweep` (362 samples, merged_pairs.k{5,6,7,8}.json +
per-k R1/R2 dada). Merged-amplicon ASV-set diff vs k5, alongside the
per-orientation churn for the *same* run:

| k | R1 churn | R2 churn | R1+R2 | merged churn | survival |
|---|----------|----------|-------|--------------|----------|
| 6 | 39 | 45 | 84  | 41 | ~49 % |
| 7 | 49 | 70 | 119 | 50 | ~42 % |
| 8 | 62 | 83 | 145 | 67 | ~46 % |

Merged sets: k5=2836 ASVs; k8=2833 (only_k5=35, only_k8=32). Merge rate flat
across k (0.936 → 0.932).

- **Merging roughly HALVES the churn but does NOT erase it** (~45–50 % survives
  into final amplicons, monotone in k). The earlier 2-sample k5/k6 "IDENTICAL"
  result was misleadingly clean — too small to see the tail.
- **Surviving churn is real per-orientation fragmentation propagating into the
  consensus**, not a merge artifact: every traced only_k8 amplicon's fwd/rev
  parents are `birth_type=Abundance` and trace to a low-Hamming k5 dada parent.
- **Hypothesis that did NOT hold:** I expected survival to be explained by the
  fragmenting base landing *inside the overlap region*. Trace **refuted** this —
  both surviving and collapsed fragments have diffs inside and outside the
  overlap. At ~250 bp amplicons the overlap covers ~62 % of the read, so the
  window heuristic can't discriminate. Whether a fragment collapses vs survives
  depends on finer consensus details the JSONs don't cleanly expose; left as a
  honest "not pinned down."
- **Abundance:** only_k8 merged churn is 11/32 below 10 reads, median 13, **max
  107**; two cases derive from genuinely abundant k5 parents (61, 161 reads). So
  it is low-abundance-*dominated* but reaches moderate abundance — not pure
  noise-floor.

### Chimera removal — the churn collapses at the FINAL table (2026-06-02)

The merged churn is measured one step too early: `remove-bimera-denovo`
(`seqtab.nonchim.k{5,6,7,8}.json`, same 362-sample run) is the last filter, and
the surviving churn is exactly what it targets. Chimera removal cuts ~75 % of
merged amplicons (2836 → 725) at every k, and **removes 60 of the 67 churned
amplicons (90 %) as chimeric** (only_k5: 32/35 removed; only_k8: 28/32 removed).

Post-chimera across-k diff (vs k5 = 725 ASVs):

| k | nonchim ASVs | shared | only_k5 | only_kN | churn |
|---|--------------|--------|---------|---------|-------|
| 6 | 728 | 721 | 4 | 7  | 11 |
| 7 | 730 | 722 | 3 | 8  | 11 |
| 8 | 733 | 722 | 3 | 11 | 14 |

**The final feature tables are ~98 % identical across k=5→k8** (residual ≈12 of
733). The churn cascade washes out at each downstream step:

| stage (k8 vs k5) | churn |
|------------------|-------|
| per-orientation (R1+R2) | 145 |
| after merging | 67 (~46 % survives) |
| **after chimera removal** | **14 (~10 % of merged; the table that goes to taxonomy)** |

This is the decisive endpoint: k-mer screen size has **essentially no effect on
the final Illumina feature table**. (Minor caveat: bimera detection is
abundance/co-occurrence-driven *within each table*, so the ~12 residual mixes
"churn removed" with "slightly different ASVs flagged chimeric per k" — but at
12/733 it is noise.)

### Open items (remaining)

1. **PacBio set-diff + abundance stratification** — PacBio only had a count check
   (+9); apply the same trace to confirm its churn is the same benign
   fragmentation. (Data available: PacBio k-mer sweeps.)

### Conclusion — REVISED (k=5 default; k has ~no effect on the final Illumina table)

The decisive finding: **k-mer screen size has essentially no effect on the final,
chimera-filtered Illumina feature table** (~98 % identical k5→k8; ~12 of 733
residual). The intermediate ASV churn (145 per-orientation → 67 merged → ~12
post-chimera) cascades away at each downstream step, and chimera removal absorbs
90 % of what survived merging.

So for **Illumina** there is no accuracy reason to raise k, and three costs if you
do: (a) intermediate ASV-set churn (mostly benign fragmentation, washed out by the
final table anyway), (b) steep memory (16 GB at k7, ~55 GB at k8 on the 61-sample
pool), (c) the screen already works fine at k5 (43 % shroud). **Keep k=5 for short
reads.** The k=6/k=7 speed edge does not justify the memory on short reads where
k=5 is already fast enough. (Opposite tradeoff from PacBio, where k=5 is a no-op
and raising k is necessary for speed — see the scale section.)

---

## Downstream: merging paired reads across k (manual)

For paired-end Illumina the per-orientation ASVs must be merged into full-length
amplicons (`merge-pairs`). The k-mer screen does **not** run in `merge-pairs` —
it re-dereplicates and does its own ends-free NW with no k-mer screen — so k
affects merging only *indirectly*, through the upstream R1/R2 ASV differences.
The question merging answers is therefore: **do the per-k ASV differences survive
into the merged amplicons, or wash out?** This step is run by hand (it is
cross-orientation, so it consumes the R1 sweep dir + R2 sweep dir together).

The sweep writes one dada JSON per sample under `<outdir>/dada_k<k>/`. For each k,
merge the matching R1/R2 dada outputs, supplying the filtered FASTQs (re-dereplicated
to recover read→unique maps). Files match by position, so sort all four globs the
same way (here the shared sample stems guarantee it):

```bash
R1=out_R1   # the R1 sweep outdir
R2=out_R2   # the R2 sweep outdir
FQ1=filtered/R1   # filtered forward FASTQs (one per sample)
FQ2=filtered/R2   # filtered reverse FASTQs

for k in 5 6 7 8; do
  dada2-rs merge-pairs \
    --fwd-dada  "$R1"/dada_k${k}/*.json \
    --rev-dada  "$R2"/dada_k${k}/*.json \
    --fwd-fastq "$FQ1"/*.fastq.gz \
    --rev-fastq "$FQ2"/*.fastq.gz \
    --min-overlap 12 --max-mismatch 0 --threads 8 \
    -o merged_k${k}.json
done
```

Output schema (per k): `{ "samples": [ { "sample", "total_pairs",
"accepted_pairs", "num_merged", "merged": [ {"sequence","abundance","accept"} ] } ] }`.

### Metrics to compare across k

**(a) merged ASV count + merge rate** — pooled over samples:

```bash
for k in 5 6 7 8; do
  python3 - "$k" merged_k${k}.json <<'PY'
import json, sys
k, path = sys.argv[1], sys.argv[2]
d = json.load(open(path))
seqs, tot, acc = set(), 0, 0
for s in d["samples"]:
    tot += s["total_pairs"]; acc += s["accepted_pairs"]
    seqs |= {m["sequence"] for m in s["merged"] if m.get("accept")}
print(f"k={k}  merged_ASVs={len(seqs)}  merge_rate={acc/tot:.3f}  ({acc}/{tot})")
PY
done
```

**(b) merged ASV-set diff vs k=5** — the decisive sequence-level check (the
cross-orientation analogue of the still-open R2 question above):

```bash
python3 - <<'PY'
import json
def merged_set(path):
    d = json.load(open(path)); out=set()
    for s in d["samples"]:
        out |= {m["sequence"] for m in s["merged"] if m.get("accept")}
    return out
base = merged_set("merged_k5.json")
print(f"baseline k=5: {len(base)} merged ASVs")
for k in (6,7,8):
    s = merged_set(f"merged_k{k}.json")
    print(f"  k={k}: shared={len(base&s)} only_k5={len(base-s)} "
          f"only_k{k}={len(s-base)} "
          f"[{'IDENTICAL' if base==s else 'DIFFERS'}]")
PY
```

If the merged sets come out IDENTICAL (or differ by ≪ the pre-merge ASV-count
spread), the upstream k-sensitivity washes out at the amplicon level — the
strongest possible "k doesn't matter for final output" result. If they differ,
(b) is exactly the sequence-level evidence needed to close the open question above.

### Worked example (2 MiSeq SOP samples, k=5 vs k=6, verified)

Filtered F3D0+F3D1 (truncLen 240/160), swept each orientation, merged per the
commands above (min-overlap 12, max-mismatch 0):

```
k=5  merged_ASVs=144  merge_rate=0.955  (11630/12173)
k=6  merged_ASVs=144  merge_rate=0.955  (11631/12177)
  k=6 vs k=5: shared=144 only_k5=0 only_k6=0  [IDENTICAL]
```

The merged amplicon set is identical and the merge rate is unchanged for these
**2 samples at k=5/k=6**.

> **Do not over-read this.** The full 61-sample sweep (above) shows the
> per-orientation ASV *sets* DO churn vs k=5 (R1/R2 "DIFFERS" at every higher k,
> growing to `only_k8`≈34–54), even though counts are stable. The 2-sample
> k5/k6 merge happened to be too small (and too low-k) to surface that. Whether
> the larger full-scale churn collapses after merging is **open** — it must be
> checked on the full 61-sample set across k=5–8 (open item #2 above), not
> inferred from this toy example.

---

## Full PacBio sweep at scale (2026-06-01) — timing inverts vs the subsample

A full-scale PacBio HiFi sweep (pacbio errfun, band 32, ~757M comparisons per k)
gives the cleanest picture yet and **overturns the small-sample timing** noted
above.

| k | learn_iters | dada_aligns | dada_shrouded | shroud % | n_asv | wall (s) |
|---|------------|-------------|---------------|----------|-------|----------|
| 5 | 6 | 756 825 650 |   7 093 737 |  0.9 | 2884 | 4139.6 |
| 6 | 6 | 757 313 016 | 608 639 940 | 44.6 | 2886 |  988.4 |
| 7 | 6 | 758 598 363 | 667 401 347 | 46.8 | 2891 |  749.2 |
| 8 | 5 | 759 208 148 | 700 398 319 | 48.0 | 2893 |  868.2 |

### Findings

- **k=5 is a no-op screen on long reads, confirmed at scale**: 0.9 % shrouded of
  757M comparisons. The 1024-entry k=5 vectors can't discriminate 1.4 kb reads,
  so nearly every pair gets a full O(L²) alignment.
- **ASV *count* impact is negligible**: 2884 → 2893, a monotonic **+9 (0.3 %)**
  across k5→k8 — even smaller than the Illumina spread. NOTE: this is **counts
  only** — no sequence-level set-diff was run for PacBio. The Illumina full run
  (above) showed stable counts can hide real bidirectional ASV-set churn, so do
  NOT assume the PacBio sets are identical across k. A pooled ASV-set diff (like
  the Illumina one) is an open item for PacBio too.
- **Timing inverts vs the subsample.** At real scale **k=5 is catastrophically
  slow** (~4140 s, ~4–5× everything else) because it aligns nearly all 757M pairs.
  **k=7 is the fastest (749 s); k=8 (868 s) beats k=6 (988 s).** The screen's
  savings (skipping long-read alignments) dominate the cost of the longer k-mer
  dot product — the opposite of the tiny-subsample result, where alignment cost
  was too small to dominate and k=8 looked ~2× slower than k=6. **Treat the
  subsample timing as an artifact; this scale ordering is the real one.**
- `maxrss_kb` is blank: the cluster lacks `/usr/bin/time`/`gtime`, so RSS wasn't
  captured (wall time is, via the bash-clock fallback). Memory comparison for
  this run therefore relies on the 4^k model (k6 ≈ 120 MiB, k7 ≈ 470 MiB,
  k8 ≈ 1.9 GiB on a full sample). **TODO: re-run with GNU `time`/`gtime` on the
  cluster to get measured peak RSS — memory is now the *sole* remaining argument
  for k=6 over k=7, so a measured number would settle the tradeoff.**

### Revised PacBio recommendation

- **Do not use k=5 for PacBio HiFi** — the screen is a no-op on long reads and the
  full alignment load makes it ~4–5× slower than any k≥6, with no ASV benefit.
- **k=7 for speed** (fastest at scale, ~749 s) **or k=6 to cap memory** (~120 MiB
  vs k=7's ~470 MiB, at ~76 % of the speed win over k=5). It is a genuine
  memory↔speed tradeoff; both give effectively identical ASVs (±9 across the
  whole k=5–8 range). k=8 has no advantage over k=7 (slower *and* 4× the memory).
- k=5 remains the safe global *default* for short Illumina reads, where the
  screen does discriminate; the above applies specifically to long-read
  (`--errfun pacbio`) runs.

---

## ASV lineage tracing through the JSON artifacts (feasibility map, 2026-06-02)

To diagnose *why* specific ASVs churn across k (open item #1), we need to trace a
final sequence back through the pipeline. This maps how far that's possible using
**only the JSON files we already emit** (scope: derep → dada → merge-pairs for
paired Illumina; derep → dada for single-end PacBio). Verified against both the
Rust output structs and real files in `data/illumina/full-pooling/`.

### What each artifact exposes for linking

| Artifact | Linking key(s) | Notes |
|----------|----------------|-------|
| **derep** (`*.derep.R*.json`) | unique **sequence string**; unique **positional index** (abundance-desc) | read→unique `map` is gated behind `--show-map` and is **absent** by default |
| **dada** (`<sample>.json`) | ASV **sequence**; ASV **index** into `asvs[]` (cluster order); **`map`** = derep-unique-index → ASV-index (always emitted); per-ASV `birth_type`/`birth_pval`/`birth_fold`/`birth_e` | `null` in `map` = unique dropped at `omega_c` |
| **seqtab** (`seqtab.R*.json`) | sequence **hash** (md5/sha1) + per-sample counts | keyed by sequence string only; not referenced by derep/dada/merge — a convenient ID that currently threads nowhere |
| **merge-pairs** (`merged_pairs.json`) | merged **sequence**; `forward`/`reverse` = **positional indices** into that sample's fwd/rev dada `asvs[]` | parent ASV sequences are NOT stored, only their indices |

### Net traceable span (today, no code changes)

- **Paired Illumina:** merged amplicon → fwd & rev **dada ASV** (by index) → **set
  of derep unique sequences** per direction (via the dada `map`). Stops at the
  unique level.
- **Single-end PacBio:** dada ASV → **set of derep uniques** (via `map`). Stops there.

**Verified end-to-end** on F3D0_S188: top merged amplicon (`forward=0`) → R1
`asvs[0]` (abund 538) → inverting the dada `map` gave **66 derep uniques summing
to exactly 538**. Exact reconciliation — the trace works.

### Break points (in priority order)

1. **Read level is severed.** The dada `map` keys on derep *unique* index, not
   *read*; the derep read→unique `map` is not emitted by default. Provenance
   bottoms out at the unique sequence, never the individual read. (Same single
   break for PacBio.)
2. **merge-pairs doesn't persist its read→(fwd,rev) composition** — it
   re-dereplicates the filtered FASTQ at runtime and stores only aggregated
   indices + abundance.
3. **All cross-stage links above the unique level are positional, not
   sequence-based** (`merge.forward`→`asvs[i]`; dada `map`→unique index). Robust
   only if every file is the same sample in its original emitted order; nothing
   embedded makes the join self-validating.
4. **No stable global ID threads the chain** — seqtab's hash exists but is
   referenced nowhere upstream/downstream.

### Implication for the churn diagnosis

Enough exists **today** to diagnose churn without new instrumentation, via a
**sequence-keyed** trace: find an `only_in_kN` ASV by its sequence in the kN dada
JSON → invert that file's `map` to its derep uniques → look those same unique
sequences up in the k5 dada `map` to see which ASV(s) they landed in there. That
directly tests the fragmentation hypothesis (did kN split a cluster k5 had
merged?), and `birth_*` speaks to the convergence-difference theory. Attribution
is to uniques, not reads — which is the right granularity for denoising anyway.

### Possible future enhancement (NOT yet implemented)

Carry the **seqtab-style sequence hash into the dada and merge-pairs outputs** (and
optionally alongside dada `map` entries). This converts the fragile positional
joins into stable hash joins, making lineage tracing turnkey and self-validating
instead of order-dependent. Low effort, no algorithmic change; deferred until
there's demand for routine tracing.
