# Band size & platform-aware defaults

**Verdict:** the platform-aware alignment band default — **16 for Illumina, 32 for
PacBio HiFi** — is *vindicated*, but not for the reason you'd expect. Band 32 on
HiFi earns its keep through **error-model stability**, not alignment reach. A
single global `BAND_SIZE` would be wrong for one platform or the other, because
the two platforms fail band-tightening through **opposite** mechanisms.

## Setup

`--band` / `BAND_SIZE` is the width of the diagonal band the Needleman–Wunsch
aligner explores. Tighter = faster/leaner but risks clipping alignments that need
to wander off the diagonal (i.e. indels). The question: can the defaults be
lowered? Tooling: `dev/run_band_sweep.sh` + `dev/compare_bands.py` (learn-errors
at fixed `--nbases`, then `dada-pooled` per band; gate = ASV set-identity **and**
reads-moved-% ≤ 0.01%, with center-aware churn labels).

## Evidence — the two platforms diverge

| Platform | k | Result | err_out max Δ |
|---|---|---|---|
| MiSeq (384-sample) | 5 | band **8 SAFE** both orientations; band **4 breaks** (R2 loses 1 ASV) | ~1e-4 |
| PacBio HiFi (74-sample) | 7 | band **cannot drop from 32**; 32→16 already fails (net −7 ASVs) | 8e-4 (b16), 1.9e-3 (b8) |

- **MiSeq**: band 8 is a safe floor (0.0013% / 0.0003% drift, zero multi-read
  churn). Band 4 is the boundary — R2 loses one knife-edge abundance call
  (birth p-value 9.8e-41 vs OMEGA_A 1e-40) that scatters into a Hamming-1/2
  neighbor.
- **PacBio HiFi**: the *opposite*. 32→16 already loses ASVs — and the lost ones
  are Abundance-born, often strongly significant (p≈0), merging into a Hamming-1
  neighbor already in the set. That's genuine loss of single-nucleotide-variant
  resolution, not a center swap.

## The two mechanisms (they move in opposite directions)

1. **Band clipping a real error-copy alignment** — driven by *indels*, not
   substitution rate (substitutions align on the diagonal at band 0; only indels
   need band ≥ gap). Illumina indels are rare and get rarer with cleaner
   chemistry, so band 8 gets *safer* on newer runs. Genuine error copies need
   only band ≤ 8 on **both** platforms (confirmed at scale — the HiFi large-band
   tail is truncation/off-target singletons, **not** homopolymer indels, which
   CCS consensus reads don't have).
2. **Borderline-ASV flip** — what actually broke both the MiSeq band-4 case and
   the whole PacBio result. The band perturbs the *learned error model*, tipping a
   knife-edge OMEGA_A abundance p-value. On 1.5 kb HiFi reads this perturbation is
   **~10–20× larger** than on MiSeq (err_out max Δ 8e-4 vs ~1e-4), because a
   long-read error model is more sensitive to how alignments resolve.

Mechanism 1 says "band ≤ 8 is fine everywhere." Mechanism 2 is why that's *not*
enough on HiFi: band 32 buys error-model stability on long reads, and that's the
binding constraint. The `band_req` proxy (how much band an alignment demands) is
**necessary but not sufficient** — only the end-to-end ASV A/B is ground truth.

## What this dictates

1. **Keep the platform-aware default (16 Illumina / 32 HiFi).** It is not
   over-provisioning; band 32 on HiFi is doing real work via mechanism 2. A single
   global band would churn one platform.
2. **Don't lower the default on chemistry-conditional evidence.** The safe floor
   moves per dataset: which ASVs sit near the abundance threshold depends on the
   error rates, so mechanism 2's boundary shifts run to run. The direct
   band-safety claim (error copies need ≤ 8) is robust; the end-to-end "is
   lowering the default safe" claim is chemistry-conditional.
3. This **refuted** an earlier `diagnostics.md` claim that "band ≤ 8 is viable on
   both platforms" — true for mechanism 1, false end-to-end on HiFi. Corrected on
   main.

## Path forward — the binned-quality wildcard

The MiSeq result is on *old* chemistry (older basecaller, more errors — near
worst-case for mechanism 1). The open validation matrix: MiSeq (done) → PacBio
(done) → **MiSeq i100 → NovaSeq**. The wildcard is **binned qualities** (i100 is
binned {12, 24, 38}; NovaSeq ~4 bins): a binned error model is a coarse step
function, so a band-induced shift moves a whole bin's rate at once — chunkier
mechanism-2 perturbations that may flip *more* borderline ASVs despite lower raw
error. Binned datasets should be swept under both `--errfun loess` and
`binned-qual` (errfun held constant across bands keeps the A/B valid; the verdict
is conditional on errfun). Watch `reads_moved_pct` and any `center->failed` label
— that's a borderline flip / lost ASV.
