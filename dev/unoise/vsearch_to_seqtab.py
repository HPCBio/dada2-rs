#!/usr/bin/env python3
"""Convert the vsearch UNOISE3 arm's output to the files run_arms.sh scores (#323).

Inputs, from the pipeline in run_arms.sh:
  --centroids  cluster_unoise centroids with ;size= (before chimera removal)
  --zotus      uchime3_denovo non-chimeras (the ZOTUs)
  --otutab     usearch_global --otutabout table of every read against the ZOTUs
  --samples    sample names in manifest order (comma-separated)

Writes to --out-dir:
  nochim.json  sample x ZOTU counts from the OTU table, in the shape
               summarize_arms.py reads (`samples`, `sequences`, `counts`)
  nochim.fa    the ZOTUs
  seqtab.json  the centroids before chimera removal, with their pooled sizes as
               one column: vsearch has no per-sample table before chimeras
  seqtab.fa    the centroids

The JSON is tagged `vsearch-unoise3`, not as a dada2-rs command, so dada2-rs
subcommands that check the tag refuse it; run_arms.sh passes them the FASTA.
"""
import argparse
import json
import os


def read_fasta(path):
    recs, label, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if label is not None:
                    recs.append((label, "".join(seq).upper()))
                label, seq = line[1:], []
            elif line:
                seq.append(line)
    if label is not None:
        recs.append((label, "".join(seq).upper()))
    return recs


def bare(label):
    """The label without ;size= or other annotations, as the OTU table uses it."""
    return label.split(";")[0]


def size(label):
    for f in label.split(";"):
        if f.startswith("size="):
            return int(f[5:])
    return 0


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for label, s in recs:
            fh.write(f">{label}\n{s}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--centroids", required=True)
    ap.add_argument("--zotus", required=True)
    ap.add_argument("--otutab", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    samples = a.samples.split(",")

    zotus = read_fasta(a.zotus)
    seq_of = {bare(lab): s for lab, s in zotus}
    with open(a.otutab) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = header[1:]
        missing = [s for s in cols if s not in samples]
        if missing:
            raise SystemExit(f"{a.otutab}: samples not in the manifest: {missing}")
        rows = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if bare(f[0]) not in seq_of:
                raise SystemExit(f"{a.otutab}: {f[0]} is not a ZOTU in {a.zotus}")
            rows[bare(f[0])] = dict(zip(cols, (int(float(x)) for x in f[1:])))
    labels = [bare(lab) for lab, _ in zotus]
    counts = [[rows.get(z, {}).get(s, 0) for z in labels] for s in samples]
    nochim = {
        "dada2_rs_command": "vsearch-unoise3",
        "samples": samples,
        "sequence_ids": labels,
        "sequences": [seq_of[z] for z in labels],
        "counts": counts,
    }
    centroids = read_fasta(a.centroids)
    seqtab = {
        "dada2_rs_command": "vsearch-unoise3",
        "samples": ["pooled"],
        "sequence_ids": [bare(lab) for lab, _ in centroids],
        "sequences": [s for _, s in centroids],
        "counts": [[size(lab) for lab, _ in centroids]],
    }
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "nochim.json"), "w") as fh:
        json.dump(nochim, fh)
    with open(os.path.join(a.out_dir, "seqtab.json"), "w") as fh:
        json.dump(seqtab, fh)
    write_fasta(os.path.join(a.out_dir, "nochim.fa"), [(bare(l), s) for l, s in zotus])
    write_fasta(os.path.join(a.out_dir, "seqtab.fa"), [(bare(l), s) for l, s in centroids])
    print(f"vsearch arm: {len(centroids)} centroids, {len(zotus)} ZOTUs, "
          f"{sum(map(sum, counts))} reads mapped")


if __name__ == "__main__":
    main()
