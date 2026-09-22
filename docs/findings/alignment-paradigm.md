# The alignment paradigm is fixed: ends-free, global, positional

**A replacement aligner has to solve the same problem, not a different one.** The
error model does not consume an alignment score or an edit distance; it consumes
a full-length positional correspondence between two reads. That requirement is
load-bearing, and it is what makes most of the alignment literature irrelevant to
this project.

The rule this yields: **a paradigm is worth evaluating if it attacks the global
ends-free problem from a genuinely different angle. It is not worth evaluating
because it is fast at a different problem.** Two candidates clear that bar —
[WFA](wfa-viability.md) and pair-HMMs — and the rest of the field does not.

## What the error model actually consumes

The constraint is visible in `compute_lambda` (`pval.rs`), and it is stronger
than "we need an alignment":

```text
λ = ∏  err_mat[transition(pos) × quality(pos)]
   pos ∈ 0..read_len
```

The product runs over **every position of the read**, not over the positions
where the two reads differ. Each position starts as a self-transition (A→A and
so on), and `Sub` overrides it only where a substitution was found. The quality
index is read at the *query* coordinate, which is recovered through
`sub.map[pos0]` — the reference-to-query position map that `al2subs` builds.

Two things follow, and they rule out whole classes of method:

- **A distance is not enough.** Edit distance, a score, a sketch similarity —
  none of them carry a `pos1` at which to read a quality. The model needs to know
  *which base* changed into *which base* at *what quality*, so score-only
  alignment modes (including WFA's default `ComputeScore`) are unusable however
  fast they are.
- **A partial alignment is not enough, and it is not merely lossy.** Every factor
  in that product is a probability below 1, so dropping positions from the domain
  makes λ *larger*. A local alignment that trims ends does not lose a little
  precision; it systematically shifts every abundance p-value toward "explainable
  as error". The bias has a direction, which is worse than noise.

## Why ends-free, and not the other two options

Amplicon reads cover the same locus, which is what makes the choice narrow:

- **Interior differences are signal.** They are the thing being classified —
  true biological variation versus a base-calling error.
- **Terminal offsets are artifacts.** Ragged ends come from trimming, primer
  removal and cycle counts, not from biology.

Pure local alignment (Smith-Waterman) would trim divergent-but-real ends, which
discards signal exactly where a variant is hardest to call. Fully global
alignment with end penalties charges for ragged ends that carry no information.
Ends-free global is the only one of the three that scores interior differences
while leaving the ends unpenalised — so it is not a tuning preference, it is the
only member of the family that matches the data.

## What indels do, and do not, do

Worth stating explicitly, because it surprises people reading `al2subs` for the
first time: **the error model is substitution-only.**

`al2subs` emits a substitution entry only when both alignment columns hold a
nucleotide. An inserted query base therefore keeps its default *match*
transition, and a deleted reference base has no query position at all. Indels
move the positional correspondence — which changes the quality each downstream
base is scored against — but they contribute no error probability of their own.

This is inherited from DADA2 deliberately, not an oversight in the port. It is
also the specific crudeness that makes a pair-HMM interesting: a generative model
would score indels and substitutions in one likelihood instead of bolting a
substitution model onto an alignment.

## The similarity assumption is load-bearing elsewhere

"These two sequences are nearly identical" is not only an assumption of the
aligner. Three other mechanisms already depend on it, and a divergence-tolerant
paradigm would invalidate all three at once:

| Mechanism | What it assumes |
|---|---|
| The [k-mer screen](kmer-size-screening.md) | dissimilar pairs can be rejected *before* alignment without loss |
| The [band](band-size-platform-defaults.md) | the optimal path hugs the diagonal, so the DP can ignore the rest of the matrix |
| The gapless shortcut (`pair_is_gapless`) | a pair with no k-mer offset shift has no indel, so the DP can be skipped |

The third one carries its own warning. That predicate used to be read off
`kdist`, which worked only while `kdist` happened to be a k-mer frequency
distance — so `kdist` was silently driving *alignment method selection*, and
swapping the screen disabled the gapless path without a word. It is now derived
from `kord` directly. See [the minimizer work](minimizer-screening.md); the
general lesson is that in this codebase the screen, the band and the aligner are
one interlocking assumption, not three independent knobs.

## What the field is actually doing (survey, 2026-06)

The conclusion that shaped the rule at the top: **global DP is not being
replaced, it is being accelerated.** WFA, KSW2's difference recurrence,
Parasail and Farrar-style striped SIMD, Edlib's bit-parallel edit distance,
GASAL2 and the GPU aligners — these are faster ways to evaluate the same
recurrence, or solutions to the distance-only problem. Neither kind changes what
`pval.rs` gets.

So "what is current in sequence alignment" mostly does not bear on the paradigm
question. It bears on the *implementation* question, which is a live and separate
one: our kernel is a vectorised banded NW that walks anti-diagonals, and [its DP
inner loop is the single largest cost in pooled denoising](compare-screen-vs-align.md).

Two candidates clear the bar:

1. **WFA** — same problem, different complexity regime: O(n·s) in edit distance
   rather than O(n·m) in read lengths, exploiting the same near-identity that
   DADA2 already relies on. Requires a traceback variant, not score-only, to feed
   `al2subs`. **Evaluated** — see [WFA as an alignment backend](wfa-viability.md).
2. **Pair-HMM** — the genuinely different angle, because it fuses alignment and
   error model into one likelihood layer instead of bolting `pval.rs` onto an NW
   alignment. Its interest is as a `learn-errors` *deliverable*, not as a faster
   aligner.

## What this dictates

- **Do not propose switching paradigms.** The productive survey space is *how to
  compute banded ends-free global alignment faster*.
- **Any candidate must emit a `Sub` with a full-length map**, with both
  nucleotides and the query quality at each substitution. That is the interface
  to check first — before benchmarking anything, since a score-only method is
  already disqualified.
- **The pair-HMM avenue is gated on cross-run transferability, not on speed.** A
  single frozen model will not travel across runs, chemistries and library preps,
  so viability depends on a hierarchical parameterisation that separates
  technology-intrinsic structure from run-specific nuisance — plus keeping
  quality-score conditioning, which is where DADA2's discriminating power comes
  from and which multiplies the parameter count. It would also deviate far enough
  from DADA2 that it may belong in a separate project rather than this one.
  Tracked in [#49](https://github.com/HPCBio/dada2-rs/issues/49).
