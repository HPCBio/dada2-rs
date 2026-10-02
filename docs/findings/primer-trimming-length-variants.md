# Primer trimming makes length variants

**Verdict:** when primers sit at more than one position in fixed-length reads,
removing them leaves the same molecules at **several lengths**: identical
starts, ends 1–4 nt apart. DADA2 cannot separate these at any depth, so each
organism becomes one ASV whose reported sequence is one of its length variants,
picked by a tie-break. Shortening every read to a common length after primer
removal (`cutadapt -l`) removes the variants without filtering out short
amplicons. Paired-end merging usually hides the problem, which is how it goes
unnoticed.

Found on the NovaSeq 6000 ITS2 deposit also described in
[Reading the prep before the result](reading-the-prep.md), through the
member-order experiment in [issue 157](https://github.com/HPCBio/dada2-rs/issues/157).

## What the reads show

Every raw R1 read in the sample below is 250 nt. After primer removal with
`dev/cutadapt_ITS.sh` (non-anchored `-g`, so the primer and anything before it
is removed), the five most common lengths are:

| length | reads | share | primer position in the raw read |
|---|---|---|---|
| 231 | 77,920 | 93.8% | offset 0: 250 − 19 |
| 232 | 3,195 | 3.8% | first primer base missing, or an indel: 18 removed |
| 228 | 384 | 0.5% | offset 3 |
| 227 | 677 | 0.8% | offset 4 |
| 250 | 365 | 0.4% | no primer found; kept untrimmed |

Sample `SRR39916535_1`, BioProject PRJNA1504839, gITS7 forward primer.

cutadapt removes the primer wherever it sits, so the insert always starts at
the same biological base. The read length is fixed, so a primer that starts
later leaves less room for the insert, and one that starts earlier leaves more.
The same molecule therefore ends at 227, 228, 231 or 232 nt, and each shorter
version is a **prefix** of the longer ones.

## Why DADA2 cannot split them

DADA2 aligns reads to cluster centres ends-free: gaps at either end cost
nothing ([alignment paradigm](alignment-paradigm.md)). A prefix variant differs
from its longer sibling only in such a gap, so its error probability against
it is close to 1 and its expected count is the whole cluster's. It can never be
significantly more abundant than expected, at any depth or `OMEGA_A`. The
variants collapse into one ASV, and which length names it is decided the way
any tie is (see [What R parity can and cannot
tell you](r-parity-floor-and-ceiling.md)).

In the member-order experiment on this data, pooled over 30 samples, every
large order-dependent flip was a pair of such variants:

| direction | lengths | reads (every arm) |
|---|---|---|
| R1 | 230 / 231 | 692 / 690 |
| R1 | 231 / 230 | 140 / 140 |
| R1 | 231 / 230 | 28 / 26 |
| R2 | 226 / 225 | 105 / 105 |
| R2 | 226 / 225 | 42 / 42 |
| R2 | 225 / 227 | 37 / 37 |

The reads are the same in every arm and in the same samples; only the length
of the reported sequence changes. Pooling makes this worse, not better: summing
samples gives each organism's variants similar totals, so more of them tie.

## Why it goes unnoticed

In a merged table the mate covers the ragged end, so the variants should merge
to one sequence, and `make-sequence-table` sums identical merged sequences, as
R's `makeSequenceTable` does. The variants then only appear in single-end,
per-direction or pooled per-direction output. **Expected, not yet measured on
this data**; earlier work on it used merged reads and did not see the
variants.

## The remedy

**cutadapt (the usual route).** After primer removal, shorten every read to the
shortest full-length peak:

```bash
cutadapt -l 227 -o R1.eq.fastq.gz R1.trimmed.fastq.gz
```

`-l` shortens only reads longer than the given length. Shorter reads, such as
short amplicons whose primer read-through was removed, keep their true length.
On synthetic reads with the primer at offsets −1, 0, +3 and +4, `-l 227` turns
the four lengths into one sequence and leaves a 150 nt read-through amplicon at
150 nt. `dev/cutadapt_ITS.sh` does this when `EQUALIZE_R1` / `EQUALIZE_R2` are
set, and prints the post-trim length distribution either way.

Also set `DISCARD_UNTRIMMED=1` where both primers are reliably present: reads
with no primer found are otherwise kept whole, primer included.

**The cost for paired-end use.** Equalising shortens most reads: on this data,
4 nt from R1 (231 to 227) and 7–8 nt from R2 (226–227 to 219). That is overlap
lost before merging, and on a length-variable amplicon the longest products,
already near the overlap limit, can stop merging and drop out of the table.
Since merging usually hides the variants anyway, equalise mainly for
single-end, per-direction or pooled per-direction analyses, and check merge
rates before and after if you equalise ahead of merging.

**Check what is left below the cutoff.** A peak that survives equalising (R2
had 0.4% at 219 nt after equalising to 224) is either a genuine short amplicon, which `-l` should
leave alone, or another primer-offset class, which means the cutoff is too
high. If those reads are prefixes of full-length reads from the same sample,
they are offset variants; lower the cutoff to them.

The prefix test needs a baseline: it requires an exact match, so a read with a
single sequencing error fails it even when it is the same molecule. Run the same
test on the full-length reads themselves (each cut to the shorter length,
leave-one-out). Here 62% of the 219 nt R2 reads matched, against 69% for the
224 nt reads, so about 90% of the peak were offset variants, and R2 was
equalised to 219. Lowering the cutoff costs nothing for genuine short
amplicons at that length, since `-l` leaves them as they are.

**`filter-and-trim`.** `--trunc-len` is **not** a substitute. Like R's
`truncLen`, it discards reads shorter than the cutoff, which on a
length-variable amplicon such as ITS removes whole taxa. Equalise with
`cutadapt -l` before or after it instead.

## How to check your own data

Compare read lengths before and after primer removal:

```bash
zcat R1.trimmed.fastq.gz | awk 'NR%4==2{print length}' | sort -n | uniq -c | sort -rn | head
```

Raw reads of one length and several trimmed peaks a few nt apart mean length
variants. The primer-offset diagnostic in `dev/cutadapt_ITS.sh` reports where
the primer sits, but the post-trim lengths are what show the consequence; on
this data the diagnostic's dominant offset (~94%) was first recorded as
constant.

## What this dictates

- **Check post-trim read lengths** on any amplicon whose primers can sit at
  more than one position: heterogeneity spacers, degenerate primers, or cutadapt
  error tolerance. It costs one command.
- **Equalise with `cutadapt -l`**, not `--trunc-len`, on length-variable
  amplicons.
- **Do not read an ASV's exact 3′ end as biology** on unequalised single-end or
  per-direction data. Its length may be the tie-break's choice.
- **When judging rs-vs-R or run-to-run differences**, a pair of ASVs one to a
  few nt apart in length, one a prefix of the other, at equal abundance is a
  length-variant rename, not a new organism.

## Provenance

- [Issue 157](https://github.com/HPCBio/dada2-rs/issues/157): the member-order
  experiment, and the pooled ITS2 flips above
- [Reading the prep before the result](reading-the-prep.md): the same deposit's
  primers, and the offset measurement this page corrects
- [Issue 240](https://github.com/HPCBio/dada2-rs/issues/240): identical merged
  sequences summed in the sequence table
