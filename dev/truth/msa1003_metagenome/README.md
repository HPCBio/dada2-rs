# MSA-1003 metagenome assembly (#166)

The input for a genome-side MSA-1003 truth set: a HiFi-only metagenome
assembly of the mock's own DNA, to compare against the reference genomes and
the alleles observed in the amplicon reads. The three public Sequel II HiFi
runs (99.1 Gb) are assembled as one sample with nf-core/mag 5.5.0.

| Run | BioProject | Library | Reads | Bases | MD5 (ENA FASTQ) |
| --- | --- | --- | --- | --- | --- |
| SRR9202034 | PRJNA546278 | ATCC_MSA_1003-CCS-10kb | 2,358,257 | 19.5 Gb | `73c7bae2c5ea54d7aea9b1653a649141` |
| SRR9328980 | PRJNA546278 | ATCC_MSA_1003-CCS-10kb | 2,419,037 | 20.5 Gb | `9acfc24801f0fec21350da6e8d494c1d` |
| SRR11469785 | PRJNA597914 | PacBio HiFi, 20 Strain Staggered | 5,646,628 | 59.1 Gb | `e3086a07086368a1c160f279f5a38ed1` |

All three share BioSample `ID964_ATCC_MSA-1003`. The 10X linked-read HiSeq X
run (SRR12283286) is left out: its barcodes do not add molecules, and its
150 bp reads cannot place a variant on an rRNA copy.

## Download

```bash
for r in SRR920/004/SRR9202034 SRR932/000/SRR9328980 SRR114/085/SRR11469785; do
  wget -c "https://ftp.sra.ebi.ac.uk/vol1/fastq/$r/$(basename $r)_subreads.fastq.gz"
done
md5sum *_subreads.fastq.gz   # compare with the table
```

## Check the reads are HiFi first

ENA names the files `*_subreads.fastq.gz`, but the read counts and the
project description (Q20+, 8.3 kb mean) say CCS. Subreads would ruin a
HiFi-mode assembly, so look before launching:

```bash
seqkit stats -a -j 3 *_subreads.fastq.gz   # expect AvgQual >= 20, N50 ~ 8-11 kb
```

`longreads_min_quality: 20` in `params.yaml` drops any read below Q20 either way.

## Run

Fill in the read paths in `samplesheet.csv`, and `outdir`, `gtdb_db` (a
release matching the pipeline's GTDB-Tk, r232 by default) and `checkm2_db` in
`params.yaml`, then:

```bash
nextflow run nf-core/mag -r 5.5.0 -profile <cluster> -params-file params.yaml
```

## Afterwards

- If annotation is needed, run Bakta on the bins: mag 5.5.0 has no Bakta
  step, and Prokka is off.
- Run barrnap on each assembly (Flye and metaMDBG) and count 16S copies per
  strain against its known operon count: fewer means the assembler collapsed
  copies.
- Compare the assembled 16S copies with the amplicon-observed variants
  (*S. epidermidis* a2 group, *B. pacificus* T4/T5, *C. beijerinckii*
  a11/a12/a14) and the current 52-allele set.
- The staggered mix leaves the rarest strains at a few-fold coverage. Where
  an assembly is missing copies, extract 16S from the HiFi reads instead:
  reads of 8-10 kb span an operon plus unique flanks, which tells copies apart.
