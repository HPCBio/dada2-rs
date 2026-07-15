# K-mer screen size

**Verdict:** the k-mer screen size (`--kmer-size`) has **~no effect on the final
chimera-filtered feature table** on either platform. Intermediate ASV churn from
changing `k` is benign cluster fragmentation that cascades away by the end of the
pipeline. So `k` is a **speed/memory knob, not an accuracy knob** — pick it for
throughput, and the platform-aware defaults follow directly.

This resolves [issue #15](https://github.com/HPCBio/dada2-rs/issues/15) ("optimal
`--kmer-size` for PacBio HiFi"), closed 2026-06-02.

## What the screen is (and isn't)

The k-mer screen is a pre-alignment prefilter: candidate uniques whose k-mer
distance to a cluster center exceeds the KDIST cutoff are skipped before
Needleman–Wunsch. It is a *speed* optimization that avoids alignments, inherited
from the ESPRIT lineage — **not** the cluster definition itself. Changing `k`
changes which pairs get aligned, which can shuffle intermediate clusters, but the
partition + abundance test downstream re-decide the actual ASVs.

## Evidence

A/B sweeps (k5 vs k8) with `compare_asvs.py`, measured at each pipeline stage —
the key methodological point being to measure at the **end**, not on intermediate
ASV counts:

| Platform | intermediate churn | after merge | post-chimera (final) |
|---|---|---|---|
| Illumina | 145 per-orientation | 67 merged | **~12 of 733 (~98% identical)** |
| PacBio | 37 raw | (single-end, no merge) | **~6 of 2132 (~99.7% identical)** |

- **Chimera removal absorbs 84–90% of the churn** — the fragmentation products
  look like bimeras and are filtered out, because a k-mer-screen fragment is
  typically an Abundance-birth, Hamming-1 from its k5 parent.
- The churn mechanism is **benign cluster fragmentation, not error-model
  convergence** — the error model is unchanged; only which reads seed which
  transient cluster moves.
- **One honest exception:** a single 375-read PacBio ASV (Hamming-3 from a
  909-read parent) survives chimera removal only at k8 — one moderately-abundant
  call (~0.3–0.5% of the table) that genuinely depends on screen size.

**Methodological lesson (reusable):** ASV *counts* are misleading for parameter
sensitivity — they stayed flat while the underlying sets churned. Only the
chimera-filtered final table is decisive. Always measure parameter sensitivity at
the end of the pipeline. This lesson recurs across the [KDIST cutoff
work](kdist-cutoff-decoupling.md) and the [band-size study](band-size-platform-defaults.md)
(use ASV set-identity, never `n_asv`).

## What this dictates

Because `k` doesn't move the final table, choose it purely for cost:

- **Illumina / short reads: keep `k=5`.** The final table is k-invariant; raising
  `k` only costs memory and churn for no accuracy gain.
- **PacBio / long reads: do *not* use `k=5`.** On full-length 16S the k=5 screen
  is effectively a no-op (almost every unique proceeds to alignment) → ~4–5×
  slower. Use **`k=7` for speed** or **`k=6` to cap memory**. (The default stays
  5 to match R DADA2's fixed `KMER_SIZE`; the recommendation lives in
  `--kmer-size` help.)
- **Memory scales as ~`4^k`.** Illumina k5/6/7/8 ≈ 4 / 6.5 / 16 / 55 GB; PacBio
  ≈ 31 / 37 / 56 / 135 GB. `KMER_SIZE_MAX = 8` is a hard u16 ceiling on the dense
  screen array.

## Path forward

Since accuracy is k-invariant, the only reason to reach for larger `k` is **screen
specificity** — a tighter prefilter lets fewer candidate alignments survive,
potentially speeding `dada`, not changing the answer. That reframes "raise the
cap" as a throughput question gated on memory, which is exactly what the
sparse-k-mer work addresses (it decouples the resident screen size from `4^k`,
making k=8 feasible). Any move past k=8, or any change to the global default,
needs an ITS/18S/non-16S validation set — all evidence here is 16S bacterial
amplicon on two chemistries.
