# Command reference

`dada2-rs <subcommand> --help` prints a one-line summary for every flag,
grouped into categories. These pages carry the detail that does not fit on one
line: what a parameter actually does, when to change it, and what the evidence
behind its default is.

## Reading `--help`

Flags are grouped under a fixed set of headings, used consistently across every
subcommand:

| Heading | What it holds |
|---|---|
| **Input** | Input files, sample naming, Phred offset |
| **Error model** | The error model file, fitting function, and its knobs |
| **Denoising** | The DADA2 algorithm parameters (R's `setDadaOpt()` surface) |
| **Alignment** | Needleman-Wunsch scoring, band radius, aligner backend |
| **Screening** | The pre-alignment k-mer screen that decides which pairs are aligned |
| **Filtering** | Read-level trimming and quality filtering |
| **Chimera** | Bimera/trimera detection parameters |
| **Metrics** | Optional per-read metrics (`summary` only) |
| **Performance** | Threads and concurrency; never affects results |
| **Diagnostics** | Extra output for QC and debugging; never affects ASVs |
| **Experimental** | Unstable flags. Off by default, and not drop-in replacements |

Two headings are worth reading as a pair. **Denoising** parameters define the
clusters; **Screening** parameters only decide which pairs are worth aligning
in the first place. A screen flag can change runtime by orders of magnitude and
can, if set too aggressively, cost you real ASVs — but it is not part of the
algorithm's definition of an ASV.

Anything under **Performance** or **Diagnostics** is safe to change without
affecting the ASV table.

## Pages

- [`summary`](summary.md) — per-position quality metrics from a FASTQ file
- [`dada`](dada.md) — per-sample denoising
- [`learn-errors`](learn-errors.md) — fit an error model

## See also

- [Parameters: `setDadaOpt()` parity](../parameters.md) — the R-equivalency table
- [Tuning for your data](../tuning-for-your-data.md)
- [Findings](../findings/index.md) — the measurements behind the defaults
