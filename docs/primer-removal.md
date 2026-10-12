# Removing primers from paired-end reads

**Use `cutadapt` to remove primers from paired-end amplicon reads, then
dada2-rs for everything after.** dada2-rs has no paired-end, sequence-aware
primer removal, and it won't get one
([issue #113](https://github.com/HPCBio/dada2-rs/issues/113), closed as won't
fix). `cutadapt` already does the job well, and it's what we use in practice.

This matters most for libraries built with **heterogeneity spacers**: 0–N nt of
padding before the primer that varies from read to read, added to raise base
diversity on patterned flow cells. With the primer at several offsets, no fixed
trim is correct, and one biological sequence becomes several ASVs. On a NovaSeq
16S run, primers left behind spacers inflated the table 3–4× and reversed the
direction of an error-model comparison
([Reading the prep](findings/reading-the-prep.md)).

## What dada2-rs can and cannot do

| Tool | Finds the primer by sequence | Paired-end | Handles spacers |
|---|---|---|---|
| [`filter-and-trim`](commands/filter-and-trim.md) `--trim-left` | No: removes a fixed number of bases | Yes | No |
| [`remove-primers`](commands/remove-primers.md) | Yes: searches the whole read, IUPAC-aware | No: one file at a time | On single-end reads only |
| `cutadapt` | Yes | Yes: keeps pairs in sync | Yes |

`remove-primers` mirrors R's `removePrimers()`, which is designed for
full-length single-end reads such as PacBio HiFi. It discards reads with no
primer match, file by file. Run separately on R1 and R2, it would discard
different reads from each file and break the pairing that `merge-pairs` needs.

## Check for spacers first

Search the raw R1 reads for the forward primer and report where it starts. Write
the primer as a regular expression, expanding IUPAC codes; this is 341F,
`CCTACGGGNGGCWGCAG`:

```bash
zcat R1.fastq.gz | head -n 400000 | awk 'NR%4==2 {
    n++
    if (match($0, /CCTACGGG.GGC[AT]GCAG/) && RSTART <= 9) { f++; off[RSTART-1]++ }
} END {
    printf "fwd primer: %d/%d (%.1f%%)\n", f, n, 100*f/n
    for (i=0; i<10; i++) if (off[i]) printf "  offset %d: %.1f%%\n", i, 100*off[i]/f
}'
```

Repeat on R2 with the reverse primer.

- **Not found:** the primers are not the ones you expected, or were already
  removed. Check before going further; deposit metadata has named the wrong
  primers before ([Reading the prep](findings/reading-the-prep.md)).
- **One offset:** a fixed `filter-and-trim --trim-left` works, though `cutadapt`
  is still fine.
- **Several offsets:** heterogeneity spacers. Use `cutadapt`.

## Remove primers with cutadapt

```bash
cutadapt \
    -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC \
    -a GGATTAGATACCCBDGTAGTC -A CTGCWGCCNCCCGTAGG \
    -e 0.15 -n 2 --max-n 0 --minimum-length 100 \
    --discard-untrimmed \
    -j 8 -o trimmed/S1_R1.fastq.gz -p trimmed/S1_R2.fastq.gz \
    raw/S1_R1.fastq.gz raw/S1_R2.fastq.gz
```

- **`-g` / `-G`** are the forward and reverse primers. They are not anchored, so
  cutadapt removes the primer **and everything before it**, which is what clears
  a spacer of any length.
- **`-a` / `-A`** are the reverse complements of the opposite primers. They remove
  read-through at the 3′ end on amplicons shorter than the read, which is common
  on ITS. They do no harm on a fixed-length 16S amplicon.
- **`--discard-untrimmed`** drops pairs where no primer was found. Use it only
  when the check above shows both primers in nearly every read; otherwise those
  reads pass through with the primer still attached.
- **`-e 0.15`** tolerates the mismatches that degenerate primers and spacers
  produce. If many reads still fail to trim, re-check the primers rather than
  raising it.

Replace the primers above (341F / 805R) with your own. cutadapt accepts IUPAC
codes directly.

## Then check read lengths

Every raw read is the same length, so a primer at different offsets leaves
trimmed reads of slightly different lengths. Several peaks a few nt apart are
the same molecules, and per-direction or single-end tables report them as
separate ASVs:

```bash
zcat trimmed/S1_R1.fastq.gz | awk 'NR%4==2 {print length}' | sort -n | uniq -c | sort -rn | head
```

To remove them, shorten each direction to its shortest full-length peak with
`cutadapt -l N`. This doesn't affect reads that are already shorter. Merging
usually hides these variants, and shortening costs merge overlap, so mainly do
this for single-end or per-direction analyses. The trade-offs are in
[Primer trimming makes length variants](findings/primer-trimming-length-variants.md).

## Then continue with dada2-rs

Use [`filter-and-trim`](commands/filter-and-trim.md) for quality filtering only
(`--max-ee`, `--max-n`, `--trunc-q`), with `--trim-left` left at 0. Many users
skip it and filter in cutadapt instead, which has an equivalent for each of its
filters.

On amplicons whose length varies, such as ITS, do not use `--trunc-len`: it
discards every read shorter than the cutoff, which removes whole taxa. Shorten
with `cutadapt -l` instead.
