# Single pass vs DADA2's divisive partition

**Verdict:** a single pass over the uniques doesn't make DADA2 faster in a way
that matters, and it doesn't make it more accurate. DADA2's divisive loop stays
the default, and both single-pass options stay experimental.

- **Speed.** Every arm compares each unique with every centre, so the cost is
  about uniques × centres in all of them. A single pass drops only the
  reshuffle between splits: 7–12 s per run on the cluster, up to 27% of dada
  time on a NovaSeq run, and nothing on PacBio, where the extra centres it
  keeps cost more than it saves.
- **The error-model pass (`sp`) keeps DADA2's accuracy.** On every truth set it
  recovers exactly the alleles `dada` does. It keeps 2–49% more ASVs, all at
  low abundance (Bray-Curtis ≤ 0.03 against `dada`). One reassignment pass
  (`--reassign`) barely changes it.
- **The UNOISE rule (`--birth-rule unoise`) drifts on real data.** It finds far
  fewer ASVs than `dada` (Jaccard 0.05–0.68). On the MiSeq mocks it misses two
  minor gene copies that `dada` recovers, because their abundance falls under
  its fixed skew ratio. Against the PacBio mock's 52 alleles it recovers the
  same 43 as `dada`, with far fewer non-truth ASVs. Most of `dada`'s extra
  ASVs there are template damage.

Tested in issue [#323](https://github.com/HPCBio/dada2-rs/issues/323)
(PR [#324](https://github.com/HPCBio/dada2-rs/pull/324)), prompted by
[#303](https://github.com/HPCBio/dada2-rs/issues/303).

## The arms

All arms share the derep input, the `learn-errors` model, the global ends-free
alignment and the k-mer screen.

| Arm | What it does |
|---|---|
| `divisive` | `dada` as it is: split one new centre at a time, reshuffle members and recompute p-values after each |
| `sp` | one pass over the uniques, most abundant first. Each becomes a centre if DADA2's abundance p-value against its best existing centre is below `OMEGA_A`, and joins that centre otherwise. Nothing is moved afterwards |
| `sp_reassign` | `sp`, then one DADA2 reshuffle (`--reassign`) |
| `unoise` | the same pass with UNOISE's rule: a unique joins a centre if its skew ≤ β(*d*) = 1/2^(2*d*+1) ([Edgar 2016](https://doi.org/10.1101/081257)). *d* counts substitutions and internal gap columns. Reads are counted as UNOISE3's `-otutab` counts them |
| `vsearch` | the UNOISE3 workflow in vsearch 2.30: merge, filter at maxEE 1, `--cluster_unoise`, UCHIME3, map reads at 97% |
| `*_indel` | the DADA2 arms with the indel-aware λ (`--indel-rates`, issue [#291](https://github.com/HPCBio/dada2-rs/issues/291)) |

`unoise` chooses exactly the same centres as `vsearch --cluster_unoise
--unoise_alpha 2` when every difference is a substitution (449 of 449 on the
MiSeq SOP forward reads). It also counts reads the same way (139,495 vs 139,476
reads). The two differ in what they count as a difference: vsearch counts
mismatches only and ignores gaps (see [Indels](#indels)). The command options
are documented on the [`dada` page](../commands/dada.md).

## Data

| Run | Platform | Samples | Reads | Truth set |
|---|---|---|---|---|
| MiSeq V4 | Illumina MiSeq, paired | 362 | 2.80M | HMP mock, 25 alleles (samples `Mock`, `Mock2`) |
| PacBio 95 | PacBio HiFi, full-length 16S | 95 | 2.39M | none |
| MSA-1003 Kinnex | PacBio Sequel IIe (full qualities) and Revio (binned) | 96 each | 2.08M / 6.33M | ATCC MSA-1003, 52 alleles |
| NovaSeq 16S | Illumina NovaSeq, paired, binned qualities | 30 | 2.45M | none |
| NovaSeq ITS2 | Illumina NovaSeq, paired, binned qualities | 30 | 2.46M | none |

All runs are pooled (`dada-pooled`), with `MIN_ABUND` unset for the DADA2 arms
and `--min-abund 8` (vsearch's `--minsize` default) for `unoise` and `vsearch`.
The driver is `dev/unoise/run_arms.sh`, summarised with
`dev/unoise/summarize_arms.py`.

## Speed

The most a single pass can save is the time spent in the reshuffle, bud and
p-value refresh. Measured with `--metrics-json` before any single-pass code
existed:

| Run | Threads | Compare | Removable | Share |
|---|---|---|---|---|
| MiSeq 362, R1 | 96 | 39.3 s | 7.3 s | 15.7% |
| MiSeq 362, R2 | 96 | 27.3 s | 12.5 s | 31.2% |
| PacBio 95, k = 7 | 96 | 171.6 s | 11.5 s | 6.3% |
| PacBio 95, k = 5 | 96 | 1013.3 s | 9.7 s | 0.9% |

The removable phases are serial and take a fixed 7–12 s, so their share only
grows with threads. Measured dada time, forward plus reverse:

| Run | `divisive` | `sp` | `unoise` |
|---|---|---|---|
| MiSeq V4 | 114 s | 102 s (−10%) | 97 s (−15%) |
| PacBio 95 | 280 s | 276 s (−2%) | 334 s (+19%) |
| MSA-1003 Sequel IIe | 75 s | 79 s (+6%) | 17 s (−77%) |
| MSA-1003 Revio | 155 s | 182 s (+17%) | 83 s (−47%) |
| NovaSeq 16S | 1631 s | 1423 s (−13%) | 1341 s (−18%) |
| NovaSeq ITS2 | 297 s | 216 s (−27%) | 182 s (−39%) |

- On NovaSeq 16S every arm screens about 15e9 pairs per direction: 1.08M
  uniques against about 14k centres. `sp` saves the 221 s `divisive` spends
  reshuffling, then spends about 5% more on comparisons because it keeps more
  centres. On PacBio the extra centres cost more than the reshuffle saves.
- `unoise` is fast only where it keeps few centres. On MSA-1003 it keeps 81
  ASVs where `divisive` keeps 1194. On diverse real communities it keeps
  about as many centres as `divisive` (13,850 vs 14,256 on NovaSeq 16S
  forward reads), so the saving is small. On PacBio 95 it is slower. It has
  no greedy locks, so it aligns pairs `divisive` skips; alignment counts
  weren't recorded for that run.
- `vsearch --cluster_unoise` takes 1–25 s on the same data. It clusters only
  the uniques that reach `--minsize`, then maps every read with an indexed
  search. That index is what makes it fast, not the single pass.

## Accuracy against a truth set

| Run | Arm | Alleles recovered | ASVs |
|---|---|---|---|
| MiSeq `Mock` | `divisive`, `sp`, `sp_reassign`, `*_indel` | 22 / 25 | 42 |
| | `unoise` | 20 / 25 | 37 |
| | `vsearch` | 20 / 25 | 45 |
| MSA-1003 Sequel IIe | every arm | 43 / 52 | 61 (`vsearch`) – 1386 (`sp_indel`) |
| MSA-1003 Revio | every arm | 43 / 52 | 71 (`unoise`) – 1211 (`sp_indel`) |

`Mock2` gives the same recall. On these data, recall depends only on the birth
rule, never on the single pass:

- **The UNOISE rule misses real minor copies.** On the MiSeq mocks, `unoise`
  and `vsearch` both lose *B. vulgatus* copy 7 and *C. beijerinckii* copy 4.
  Each is one substitution from a major copy of its strain, at 8–11% of that
  copy's reads in `Mock`. A unique one difference from a centre is absorbed below a
  skew of β(1) = 1/8, so they're absorbed. DADA2's p-value is driven by the
  error model, not a fixed ratio, so it keeps them.
- **MSA-1003 doesn't separate the arms.** All recover the same 43 alleles on
  both instruments. The 9 misses are single-copy alleles where the mock's
  strains most likely differ from the reference assemblies (issue
  [#166](https://github.com/HPCBio/dada2-rs/issues/166)), so precision measured
  against this set undercounts every arm by the same amount.
- **Most of DADA2's extra ASVs on MSA-1003 are damage, not noise.** Of
  `divisive`'s 1101 near-truth ASVs on Sequel IIe, 1027 are C>T/G>A variants
  at 0.03–0.2% of a parent; on Revio, 550 of 609. That is
  [template damage](template-damage-variants.md). UNOISE's ratio absorbs it,
  which is why its ASV counts are so low there.

### Damage flagging as a secondary row

[`flag_artifacts.py --action merge --include-damage`](../scripts/flagging.md)
merges damage and slippage candidates into their parents after denoising, as
chimera removal does. It is **experimental**: the evidence that these variants
are damage comes from mocks only. On MSA-1003 Sequel IIe it reduces `divisive`
from 1194 to 120 ASVs and raises the reads on truth alleles from 82.5% to
87.3%, close to `vsearch`'s 86.8%, with recall unchanged. On the NovaSeq runs it
flags no damage at all. Report bare denoisers first and merged results
separately.

## Agreement with `divisive` on real data

ASV set identity after chimera removal: Jaccard, with Bray-Curtis on read counts
in brackets.

| Run | `sp` | `sp_reassign` | `unoise` | `vsearch` |
|---|---|---|---|---|
| MiSeq V4 | 0.91 (0.004) | 0.91 (0.003) | 0.59 (0.034) | 0.39 (0.083) |
| PacBio 95 | 0.96 (0.003) | 0.96 (0.002) | 0.65 (0.037) | 0.66 (0.022) |
| MSA-1003 Sequel IIe | 0.89 (0.006) | 0.89 (0.006) | 0.05 (0.077) | 0.05 (0.060) |
| MSA-1003 Revio | 0.66 (0.013) | 0.66 (0.013) | 0.08 (0.070) | 0.08 (0.050) |
| NovaSeq 16S | 0.82 (0.027) | 0.84 (0.025) | 0.43 (0.155) | 0.31 (0.144) |
| NovaSeq ITS2 | 0.97 (0.008) | 0.97 (0.008) | 0.68 (0.107) | 0.73 (0.053) |

- **`sp`'s extra ASVs are low-abundance near neighbours.** On NovaSeq 16S,
  `sp` has 6,415 ASVs that `divisive` lacks before chimera removal, with a
  median of 2 reads and 2.7% of reads in total. Half are one difference from
  a `divisive` ASV.
- **Reassignment can't remove them.** A reshuffle moves members between
  existing centres but never removes one, so `sp_reassign` has the same
  denoised ASVs as `sp`. Its sets are identical on PacBio and ITS2. The 12%
  shift on NovaSeq 16S comes from merging the reshuffled reads.
- **The low-abundance end is unstable whatever changes.** Adding the indel
  model to `divisive` moves NovaSeq 16S by about as much as `sp` does (Jaccard
  0.81), and the median ASV shared by `sp` and `divisive` there has 3 reads in
  total.

### `unoise` and `vsearch` disagree on paired reads

On PacBio, `unoise` and `vsearch` agree closely (Jaccard 0.97 on PacBio 95). On
paired Illumina reads they don't (0.40 on NovaSeq 16S), and `unoise` keeps
fewer reads: 67% of input on NovaSeq 16S against `vsearch`'s 77%. The
difference is the workflow, not the rule. DADA2 denoises each direction, then
merges, so a pair is lost if either read has no centre within 97%; on NovaSeq
16S that is 15% of forward and 19% of reverse reads. UNOISE3 merges first, so
each read has one chance at 97% over its full length. A merge-first `unoise`
arm would separate the two. DADA2 merges after denoising because its error
model ties error rates to quality scores, which differ between the overlap and
the single-read ends (benjjneb/dada2#2163). The UNOISE rule doesn't use quality
scores, so that objection doesn't apply to it
([#327](https://github.com/HPCBio/dada2-rs/issues/327)).

## Indels

The indel-aware λ (`--indel-rates`) is a separate axis, run alongside each
DADA2 arm because UNOISE's *d* counts gaps by definition. Neither DADA2's
default split rule (substitutions only) nor vsearch's *d* (mismatches only)
counts a difference made of gaps alone:

- **MSA-1003:** *B. pacificus* allele 1 with a T4 run read as T5 is in all 96
  samples on both instruments, at 37–39% of its parent. `unoise` and the indel
  arms keep it; `divisive` and `vsearch` merge it into its parent.
- **NovaSeq ITS2:** the indel arms keep 61–83 indel-only variants, and `unoise`
  keeps 112, holding about 2.7% of reads. `flag_artifacts.py` finds no
  artifact evidence for them: they don't track a parent across samples.
  `divisive` keeps 1 and `vsearch` 4.

## Consequences

- **Default unchanged.** `dada` keeps the divisive partition.
  `--partition single-pass`, `--reassign` and `--birth-rule unoise` remain
  experimental, for comparison only.
- **Don't retry a single pass for speed.** Its saving is bounded by the serial
  reshuffle and can turn negative when more centres are kept. Speed work
  belongs in the comparison step, which every partition pays in full: a
  cheaper screen or an index, as vsearch uses.
- **Don't adopt a fixed abundance ratio.** UNOISE's β(*d*) absorbs real minor
  copies on Illumina mocks where DADA2's error model keeps them.
- **Open:** a merge-first `unoise` arm
  ([#327](https://github.com/HPCBio/dada2-rs/issues/327)); a truth set observed
  from the MSA-1003 reads ([#166](https://github.com/HPCBio/dada2-rs/issues/166)),
  without which precision on that mock can't separate the arms.

## Caveats

- The MiSeq, PacBio 95 and MSA-1003 runs predate the `unoise` counting change
  (99d8fbc), which counts each unique toward its nearest centre at 97%. It
  leaves the centres unchanged. On NovaSeq it raised `unoise`'s reads by 4–20%
  and, through new forward/reverse pairings, added low-abundance merged ASVs,
  so the older runs' `unoise` read totals and paired-read ASV sets are low.
- `vsearch` times cover `--cluster_unoise` only, not merging or mapping.
- NovaSeq has no mock, so its numbers measure agreement, not accuracy.
