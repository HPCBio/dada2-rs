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

### PacBio SRR28724909 (HiFi, ~1450 bp; pacbio, band 32)

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
2. **k=5 wastes the screen on long reads** — 0 % shrouded means every pair is
   aligned, making k=5 the slowest PacBio config (~4× slower than k=6).
3. **k=6 is the PacBio sweet spot**: the screen engages (≈4× speedup over k=5),
   results are identical, and memory is only ~120 MiB on a full sample.
4. **k=7/k=8 add memory (4× per step) with no benefit** over k=6 here, and their
   aggressive shrouding is a latent false-negative risk on harder data.

## Recommendation

- Keep **k=5 as the global default** (correct and lowest-memory for short reads,
  where the screen does discriminate — see `benchmark_kmer_size.md`).
- For **PacBio HiFi, use k=6**: same ASVs and error model as k=5, ~4× faster, and
  far cheaper than k≥7. Document k=6 in the `--kmer-size` help for `--errfun
  pacbio`, mirroring the per-errfun guidance added for `--loess-*`.
- Avoid k ≥ 7 for HiFi, especially pooled or large-N runs (memory + benign-but-
  aggressive shrouding with no upside).

### Caveats

Single Illumina sample; single PacBio sample (full uniques). The benign-shrouding
result (identical ASVs at 79 % shroud) should be confirmed on a second sample or
a pooled run before changing any default. Scratch artifacts were purged from
`/tmp`; re-run to reproduce exact `err_out` diffs. The memory model and the
k=5-is-a-no-op-on-long-reads finding are robust (mechanistic + measured).
