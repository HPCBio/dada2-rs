# Artifact flagging

`flag_artifacts.py` flags ASVs that are probably artifacts of a more abundant
ASV, after denoising and chimera removal. It is the same kind of step as
chimera removal: the denoiser judges each unique against one centre, while the
evidence that a variant is an artifact lies across the whole run. On PacBio
HiFi mocks, flagging **doubles precision with no loss of recall**
([Template damage makes C>T/G>A variants](../findings/template-damage-variants.md)).

```bash
python3 scripts/flag_artifacts.py nochim.json --out flags.tsv
python3 scripts/flag_artifacts.py nochim.json --out flags.tsv \
    --action remove --out-seqtab nochim.filtered.json
```

The input is a `make-sequence-table` or `remove-bimera-denovo` JSON (optionally
`.gz`). Python 3, standard library only. **Experimental**: the labels and
thresholds are set from PacBio full-length 16S and may change.

## What it flags

Three classes of candidate, each found without alignment:

| Class | Candidate | Evidence for an artifact |
|---|---|---|
| Homopolymer length | differs from a more abundant ASV only in run lengths, at most `--max-edits` bases (default 2) | the same edit at the same site recurs in other parents; tracking |
| One indel | one inserted or deleted base outside a run, plus any run-length changes | as above |
| Damage | one C>T or G>A substitution, at most `--damage-max-frac` (0.005) of its parent | at least `--damage-min-siblings` (3) such candidates of the same parent; tracking |

Each candidate's parent is the most abundant ASV it qualifies against.
**Tracking** means the candidate is a near-constant fraction of its parent in
every sample: Pearson dispersion near 1 and no sample with the candidate but not
its parent. A candidate that varies independently of its parent (dispersion
above 10, or any such orphan sample) is a `real-candidate` whatever else holds.
A single-sample table has no tracking evidence.

Damage is cytosine deamination in the template: a damaged molecule is copied
faithfully, so its reads carry the change at high quality and the denoiser
correctly keeps it apart. One template accumulates damage at many sites, so
the siblings are the evidence. The Illumina run tested gets no damage calls.

## Labels

| Label | Meaning |
|---|---|
| `artifact` | homopolymer or one-indel candidate: parent at least `--min-fold` (4) times more abundant, every edit explained by recurrence or the context rate, tracking does not argue against it |
| `damage` | damage candidate with enough siblings, tracking does not argue against it |
| `real-candidate` | varies independently of its parent, or is far above the rate its context shows elsewhere |
| `unresolved` | none of the above |

`--action remove` drops `artifact` and `damage` ASVs; `--action merge` adds
their reads to their parent, or to the first kept ancestor when the parent is
dropped too. Both need `--out-seqtab`. The default, `flag`, writes only the
TSV. `real-candidate` and `unresolved` ASVs are never removed.

## Output

One TSV row per candidate, most abundant first:

| Column | Content |
|---|---|
| `label` | as above |
| `reads`, `parent_reads`, `fold` | totals over all samples, and their ratio |
| `frac` | candidate / (candidate + parent), pooled over samples |
| `samples`, `dispersion`, `orphans` | tracking: samples with either, Pearson dispersion, samples with ≥2 candidate reads and no parent |
| `edits` | edited bases from the parent |
| `evidence` | per edit: change, context, recurrence (`recur=`) and context rate (`ctx=`); for damage, substitution, position, context and `siblings=` |
| `parent`, `sequence` | the two sequences |

A summary line on stderr counts each label and its reads.

## Options

| Flag | Default | Effect |
|---|---|---|
| `--max-edits` | 2 | run-length and indel bases a homopolymer or one-indel candidate may differ by |
| `--min-fold` | 4 | parent / candidate abundance needed for `artifact` |
| `--context` | 6 | bases either side that define an edit's site |
| `--max-ratio` | 4 | how far above a recurring or context fraction an edit is still explained |
| `--pair-window` | 10 | edits of opposite sign this close are a substitution, not slippage, so not a candidate |
| `--min-unexplained` | 0.1 | fraction above which an edit with no context evidence is `real-candidate` |
| `--damage-max-frac` | 0.005 | highest fraction of its parent for a damage candidate |
| `--damage-min-siblings` | 3 | damage candidates of one parent, itself included, needed for `damage` |
| `--rates`, `--min-exact` | — | an optional per-run slippage table, giving each homopolymer context's rate including variants the denoiser absorbed; maintainer tooling, so not covered here |
| `--action`, `--out-seqtab` | `flag` | as above |

## Limits

- Only one-indel and one-substitution relationships are searched. Multi-switch
  PCR recombinants between rRNA copies differ by several indels and are not
  caught (issue [#55](https://github.com/HPCBio/dada2-rs/issues/55)).
- A real allele can look like damage if it is one C>T or G>A from a much more
  abundant allele that also has damage siblings and tracks it closely. On the
  mocks no true allele was flagged.
- Tested on PacBio full-length 16S (MSA-1003 mocks, a 95-sample real run) and
  one Illumina V7–V9 run; other amplicons and platforms are untested.

Developed under issue [#291](https://github.com/HPCBio/dada2-rs/issues/291).
