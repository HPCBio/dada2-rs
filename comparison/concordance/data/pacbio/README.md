# PacBio data

Subsampled PacBio HiFi 16S reads (raw, primered, single-end) for the concordance
guardrail — small enough that CI runs in seconds.

`SRR8557463.fastq.gz`, `SRR8557464.fastq.gz`: the first 1500 reads of two samples
from the SRA-based PacBio Sequel IIe set. Regenerate with:

```bash
for s in SRR8557463 SRR8557464; do
  gzcat data/pacbio-sqii/Raw_FASTQ/${s}.sample.fastq.gz | head -n 6000 \
    | gzip -6 > comparison/concordance/data/pacbio/${s}.fastq.gz
done
```

Reads are ~1495 bp; primers are 27F (`AGRGTTYGATYMTGGCTCAG`) / 1492R
(`RGYTACCTTGTTACGACTT`). `run_pacbio.sh` denoises this to ~55 ASVs in a few
seconds. If you change the above data, regenerate `reference/pacbio_seqtab_nochim.csv` with `write_reference.R` on the same files.
