# PacBio data

Subsampled PacBio HiFi 16S reads (raw, primered, single-end) for the concordance
guardrail — small enough that CI runs in seconds.

`SRR8557463.fastq.gz`, `SRR8557464.fastq.gz`: the first 5000 reads of two samples
from the SRA-based PacBio Sequel IIe set. Regenerate with:

```bash
for s in SRR8557463 SRR8557464; do
  gzcat data/pacbio-sqii/Raw_FASTQ/${s}.sample.fastq.gz | head -n 20000 \
    | gzip -6 > dev/concordance/data/pacbio/${s}.fastq.gz
done
```

Reads are ~1495 bp; primers are 27F (`AGRGTTYGATYMTGGCTCAG`) / 1492R
(`RGYTACCTTGTTACGACTT`). `run_pacbio.sh` denoises this to ~93 ASVs in ~20-30s.

**Depth note:** 5000 reads/sample was chosen because at 1500 the data was too
shallow to fit a stable PacBio error model — a few mid-abundance ASVs were
denoising-boundary artifacts (e.g. a Hamming-2 variant dada2-rs error-corrected
that R budded), which resolved once depth increased. At low depth neither tool is
ground truth there, so deeper toy data gives a more meaningful concordance
baseline. If you change this data, regenerate `reference/pacbio_seqtab_nochim.csv`
with `write_reference.R` on the same files.
