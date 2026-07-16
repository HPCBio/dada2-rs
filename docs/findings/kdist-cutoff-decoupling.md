# KDIST cutoff decoupling

**Verdict:** the k-mer-distance (KDIST) cutoff used by `learn-errors` and the one
used by `dada` can be set independently, and they *should* be. Most of the
speedup from a tighter cutoff buys in the **dada stage**, and that portion is
safe — it leaves the final ASV table essentially unchanged. The learning stage
speedup is the risky part: pushing the same cutoff into `learn-errors`
perturbs the error model enough to churn real-abundance ASVs. **How safe the
dada-stage tightening is, and how far you can push it, is platform-dependent** —
free on Illumina MiSeq, but subject to an irreducible borderline-ASV churn floor on
PacBio HiFi (see the PacBio section). 

In the benchmarking script, `--learn-kdist-cutoff`
lets you set the error model for `learn-error` independently (e.g., to 0.42) while still tightening `dada`.  When running the workflow steps independently, 
you can independently set `--kdist-cutoff` for `learn-errors` and for 
`dada-pooled`/`dada-pseudo`/`dada` (thus leaving `learn-errors` at the default and the denoising/inference step using a tighter cutoff).

## Background

The KDIST cutoff is the k-mer-distance threshold above which a candidate unique
is screened out before Needleman–Wunsch alignment. A tighter cutoff screens more
aggressively → fewer alignments → faster. Prior work
(*kdist cutoff + error-model stability*) established that the cutoff does **not**
act by screening out error copies (the final sub-alignment re-includes them);
it acts by altering the **partition / error model**. That predicted a clean 
alternative: decoupling the cutoff `learn-errors` sees from the one `dada` sees.
This benchmark tests that prediction.

## Setup

Illumina MiSeq, F1000 `MiSeqSOP` full data set (384 samples), full pooled mode, 
k=5, 24 threads on a node-exclusive run. Three arms, one run each, same binary
(`dada2-rs 0.2.0-ed578190` — this is a *settings* A/B, not a code A/B):

| Arm | Flags | learn cutoff | dada cutoff |
|---|---|---|---|
| `main` | (defaults) | 0.42 | 0.42 |
| `fixed_kdist` | `--kdist-cutoff 0.30` | **0.30** | **0.30** |
| `decoupled_kdist` | `--kdist-cutoff 0.30 --learn-kdist-cutoff 0.42` | 0.42 | **0.30** |

`fixed_kdist` tightens *both* stages; `decoupled_kdist` tightens **only** `dada`,
holding `learn-errors` at the default 0.42.

## Performance

Wall time, Δ% vs `main` (↓ better):

| Step | main | fixed_kdist | decoupled_kdist |
|---|---:|---:|---:|
| learn_fwd | 18.4s | **−52.3%** | −0.8% |
| learn_rev | 21.7s | **−57.5%** | −0.0% |
| dada_fwd | 163.3s | −30.1% | **−32.1%** |
| dada_rev | 130.2s | −24.8% | **−25.6%** |
| **TOTAL** | 354.8s | **−29.2%** | −24.3% |

- `decoupled` leaves `learn-errors` untouched (learn timings ≈ `main`), so its
  entire win comes from the **dada stage** (−32% / −26%).
- `fixed` adds the learn-stage win on top (−52% / −57%), reaching the larger
  −29.2% total — but that extra ~5 points is exactly what moves the error model
  (see below).

Memory is flat across all three (pipeline peak lives in `merge`, not `dada`;
within ±3%).

## Correctness (ASV concordance)

Using `dev/compare_asvs.py` on the final chimera-filtered table, baseline `main`
(725 ASVs):

| Arm | ASVs | churn | dropped from main | novel | novel ≥10 abund |
|---|---:|---:|---:|---:|---:|
| `fixed_kdist` | 727 | 6 | 2 | 4 | **2** (abund 31, 23) |
| `decoupled_kdist` | 726 | **1** | **0** | 1 | 0 |

- **`decoupled` loses nothing.** All 725 of `main`'s ASVs survive
  (`only_main = 0`). It adds exactly one new ASV (abundance 4, Hamming-25 from its
  nearest neighbor — genuinely distinct, not cluster fragmentation). Effectively
  the same table.
- **`fixed` churns real signal.** Pushing 0.30 into `learn-errors` drops 2 of
  `main`'s ASVs and adds 4, including **two meaningful-abundance novel ASVs**
  (31, 23) sitting Hamming-4/7 from high-abundance neighbors — the signature of
  an error-model-driven partition shift, not new biology.

## What this suggests:

1. **The dada-stage KDIST win is safe to take.** Tightening the cutoff `dada`
   sees buys ~−32%/−26% on the dada steps with a near-identical ASV table. This is
   the portion worth exposing as a recommended tuning.
2. **Do not tighten `learn-errors` for speed.** The learn-stage speedup is the
   ~5-point tail that perturbs the error model and churns real-abundance ASVs.
   Keep `--learn-kdist-cutoff` at the default.
3. **`--learn-kdist-cutoff` is the clean lever.** Decoupling the two cutoffs —
   exactly as the error-model-stability work predicted — lets a user take the
   dada speedup without paying the concordance cost. This validates keeping the
   two knobs separate rather than collapsing them into one.

## The one novel ASV, characterized

The single sequence `decoupled` adds is worth chasing down, because it decides
whether the "churn=1" is a benign gain or a subtle artifact. It appears to be 
a **genuinely distinct, real-looking rare taxon** — additive promotion, 
exactly the mode the mechanism decomposition predicts for a tight-inference / frozen-error-model run:

- **253 bp**, a clean modal-length V4 amplicon — no length anomaly, no
  homopolymer/indel signature.
- **4 reads across 3 independent samples** (1 / 2 / 1). A localized PCR/sequencing
  artifact would concentrate in one sample; spread across three unrelated samples
  argues for a real community member.
- **Hamming-25 (~90% identity) from its nearest same-length `main` ASV**
  (abundance 3684) — a *different organism*, categorically not a 1–2 base
  fragmentation of an existing center.
- **Not a chimera** — it survives `remove-bimera-denovo` in the decoupled arm.
- **Never born in `main` at all** — absent even from `main`'s *pre*-chimera table
  (2836 ASVs vs decoupled's 2837). `main` did not infer-then-absorb it; the looser
  0.42 dada cutoff simply never split these 4 reads into their own partition. The
  tighter 0.30 cutoff let them form a cluster that then cleared the abundance
  p-value.

So on this dataset the tighter dada cutoff didn't merely preserve `main`'s table —
it recovered one additional plausible rare taxon at the cost of nothing (zero ASVs
lost). The only caveat is depth: at 4 reads it sits right at the abundance-p-value
borderline, a real-but-marginal call that would firm up or drop out with more
sequencing.

## PacBio HiFi: the same lever, a different floor

The MiSeq result above says the dada-stage tightening to 0.30 is essentially free.
On PacBio HiFi it is **not** free — but the decoupling machinery is still the right
tool, and the experiment surfaces a sharper sub-finding about *where* the
learn-stage stability edge sits.

Same 3-arm design, on the standing PacBio HiFi benchmark set (full pooled, **k=7** —
already de-saturating long-read k-mer distances), same binary
(`dada2-rs 0.2.0-ed578190`). Because 0.30 turned out to churn real ASVs (below), a
second round was run at a looser 0.35. The `main` default (0.42) is the shared
baseline; its 696s dada / 757s total reproduce the standing benchmark (699s) to
within run noise.

### Performance (0.35 round, Δ% vs `main`)

| Step | main | fixed 0.35 | decoupled 0.35 |
|---|---:|---:|---:|
| learn | 24.4s | −23.6% | +1.0% |
| dada | 696.0s | −37.7% | **−38.0%** |
| **TOTAL** | 757.0s | −35.5% | **−34.9%** |

At the **tighter 0.30** the dada win is larger — decoupled dada −49.3%, total
−45.4% — so the cutoff is a direct speed/concordance dial: 0.30 buys ~10 more points
of total speedup than 0.35. `dada` dominates the pipeline (696s of 757s) and `learn`
is only ~24s, so on PacBio the learn-stage tightening is worth almost nothing at the
total level regardless — the entire win comes from `dada`. Memory is flat across all
arms (peak lives in dada's ~17.9 GB coexistence).

### Correctness — the learn-stability edge sits between 0.30 and 0.35

Post-chimera table, baseline `main` (2082 ASVs), via `dev/compare_asvs.py`:

| Arm | churn | drops ≥10 abund | novel set (abund / minH to nearest) |
|---|---:|---:|---|
| fixed 0.30 | 7 | 3 | 76/1, 12/7, 4/1177 |
| decoupled 0.30 | 5 | 3 | 76/1, 8/5 |
| fixed 0.35 | 6 | 2 | 76/1, 12/7, 4/1177 |
| decoupled 0.35 | 6 | 2 | 76/1, 12/7, 4/1177 |

Two things jump out:

1. **At 0.35, `fixed` and `decoupled` are byte-for-byte identical** (same 2082 ASVs,
   same churn, same novel set) — and this holds **pre-chimera too** (both 2815 ASVs,
   churn=8, identical novels), so it is not a chimera-filter coincidence. The only
   difference between those arms is the learn cutoff (0.35 vs 0.42), so an identical
   table means **the error model is stable across learn ∈ [0.35, 0.42]**. At 0.30,
   `fixed` ≠ `decoupled` (churn 7 vs 5), so learn=0.30 *does* perturb the model. The
   **learn-stage stability edge therefore sits between 0.30 and 0.35**:
   `--learn-kdist-cutoff` only earns its keep once you drive the dada cutoff below
   ~0.35. At 0.35 you can safely couple; below it you must decouple to protect learn.

2. **The dada-stage churn is irreducible and non-monotonic.** Every novel ASV at
   0.35 appears in the `decoupled` arm too, where learn is untouched at 0.42 — so
   they are 100% dada-stage effects, not learn perturbation. Loosening 0.30→0.35 did
   not clean them up: churn ticked *up* (5→6) and the novel set changed *identity*
   (the abund-8 at 0.30 vanished; the 12 and 4 appeared) rather than shrinking. The
   abund-76 / Hamming-1 sequence is a single mid-abundance ASV that repositions one
   base under **any** dada tightening — it is present in all four tightened arms.
   0.35 does recover one real-abundance drop (3→2 ≥10), a marginal correctness gain,
   at the ~10-point speed cost noted above.

### Does chimera removal explain the churn? Partly.

Pre-chimera churn is 8; post-chimera it is 6. So `remove-bimera` absorbs 2 of the
8 churners (one per side) — and the one it eats is a textbook bimera candidate
(abund-4, 4 bases off a neighbor of abundance 4). But the **core churners survive
chimera removal**: the abund-76/H1 flip, abund-12/H7, and abund-4/H1177 persist.
Those are genuine near-threshold partition differences, not chimera-call
disagreements — the chimera filter shaves the edge but does not explain the residual.

### What this dictates (PacBio)

- **No dada tightening is free on PacBio.** There is an irreducible ~0.3%
  (churn 5–6 of 2082) borderline-ASV floor near the OMEGA_A abundance threshold,
  non-monotonic in the cutoff. This is the "PacBio is prone to noisy borderline
  ASVs" phenomenon, and it is really an **OMEGA_A question, not a kdist question**
  (see Remaining threads) — kdist tuning cannot remove it.
- **0.30 is the speed play**: total −45%, but decoupling is *mandatory* (learn=0.30
  perturbs the model) and it costs ~2–3 near-threshold real-abundance ASVs.
- **0.35 is the conservative play**: total −35%, coupling is safe (learn stable
  ≥0.35), costs ~2. Either is defensible; the churned ASVs are all marginal
  near-threshold calls.
- **The safe cutoff is genuinely platform-specific**, and not merely a k-mer
  saturation artifact — k=7 was already de-saturating here, and 0.30 still churned.
  The kdist calibration distribution shifts with both platform and k, so the MiSeq
  0.30 was never guaranteed to transfer.

## Remaining threads

Confirm the pattern holds on binned-quality platforms (i100 / NovaSeq /
PacBio Revio), where the error-model sensitivity to the cutoff could differ again.

**Open question — platform-appropriate OMEGA_A/OMEGA_C.** On PacBio the churn
introduced by a cutoff change concentrates in the low-abundance band (≈4–76 reads),
i.e. right at the OMEGA_A abundance-p-value borderline. The abundance/error
thresholds (`OMEGA_A`, `OMEGA_C`) are fixed at `1e-40` and platform-blind, inherited
from DADA2's Illumina-calibrated defaults; HiFi's error structure is very different.
Whether they should differ for HiFi is a real question, but it is **not answerable
from concordance-vs-reference A/Bs** — those measure whether an ASV moved, not
whether it was correct. It requires a PacBio mock-community truth set (the same
resource that gates the raise-k question). Scoped and parked until that exists.