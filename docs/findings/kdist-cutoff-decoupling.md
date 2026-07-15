# KDIST cutoff decoupling

**Verdict:** the k-mer-distance (KDIST) cutoff used by `learn-errors` and the one
used by `dada` can be set independently, and they *should* be. Most of the
speedup a tighter cutoff buys lives in the **dada stage**, and that portion is
safe — it leaves the final ASV table essentially unchanged. The learn-stage
speedup is the risky part: pushing the same tight cutoff into `learn-errors`
perturbs the error model enough to churn real-abundance ASVs. `--learn-kdist-cutoff`
lets you keep the error model at its default while still tightening `dada`.

## Background

The KDIST cutoff is the k-mer-distance threshold above which a candidate unique
is screened out before Needleman–Wunsch alignment. A tighter cutoff screens more
aggressively → fewer alignments → faster. Prior work
(*kdist cutoff + error-model stability*) established that the cutoff does **not**
act by screening out error copies (the final sub-alignment re-includes them);
it acts by biasing the **partition / error model**. That predicted a clean lever:
decouple the cutoff `learn-errors` sees from the one `dada` sees. This benchmark
tests that prediction.

## Setup

Illumina MiSeq, F1000 `MiSeqSOP` full data set (384 samples), pooled mode, k=5,
24 threads on a node-exclusive run. Three arms, one run each, same binary
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

`dev/compare_asvs.py` on the final chimera-filtered table, baseline `main`
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

## What this dictates

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

The single sequence `decoupled` adds is worth pinning down, because it decides
whether the "churn=1" is a benign gain or a subtle artifact. It is a **genuinely
distinct, real-looking rare taxon** — additive promotion, exactly the mode the
mechanism decomposition predicts for a tight-inference / frozen-error-model run:

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

## Open thread

Confirm the pattern holds on binned-quality platforms (i100 / NovaSeq), where the
error-model sensitivity to the cutoff could differ.
