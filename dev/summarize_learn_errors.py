#!/usr/bin/env python3
"""Summarize verbose output from dada2-rs `learn-errors` or R dada2::learnErrors().

Both implementations emit the same per-cluster trace lines

    New Cluster C<i>:SSS, Division (naive): Raw <r> from Bi <b>, pA=<p>
    ... , Division (prior): Raw <r> from Bi <b>, pP=<p>
    ... , No Division. Minimum pA=<p> (Raw <r> w/ <n> reads in Bi <b>).
    ALIGN: <n> aligns, <m> shrouded (<r> raw).

so this parser works on either.  Each `ALIGN:` line ends one dada run
(one sample within one self-consistency step); the script segments on
those and reports, per run, how many clusters were budded, how the buds
split between naive/prior divisions, how many shuffles ran, and the
alignment volume.  That makes it easy to diff over-budding and alignment
counts between the Rust and R implementations on identical input.

Usage:
  # dada2-rs (verbose trace goes to stderr)
  dada2-rs learn-errors --verbose ... 2> rs.log
  python3 dev/summarize_learn_errors.py rs.log

  # R dada2 (capture messages; verbose=2 for the per-cluster trace)
  Rscript -e 'dada2::learnErrors("samp.fastq.gz", verbose=2)' 2> r.log
  python3 dev/summarize_learn_errors.py r.log

  # or pipe directly
  dada2-rs learn-errors --verbose ... 2>&1 | python3 dev/summarize_learn_errors.py
"""

import argparse
import re
import sys

# Shared trace tokens (identical between dada2-rs and R dada2).
RE_NEW_CLUSTER = re.compile(r"New Cluster C(\d+):(S*)")
RE_NAIVE = re.compile(r"Division \(naive\)")
RE_PRIOR = re.compile(r"Division \(prior\)")
RE_NODIV = re.compile(r"No Division")
RE_ALIGN = re.compile(r"ALIGN:\s*(\d+)\s*aligns,\s*(\d+)\s*shrouded\s*\((\d+)\s*raw\)")

# dada2-rs learn-errors wrapper lines.
RE_RS_SAMPLE = re.compile(r"iter=(\d+)\s+sample=(\d+):\s*(\d+)\s*cluster")
RE_RS_ITER = re.compile(r"iter=(\d+):\s*max \|err_in - err_out\| = (\S+)")

# R dada2 wrapper line (learnErrors self-consistency loop).
RE_R_STEP = re.compile(r"selfConsist step (\d+)")


class Segment:
    """One dada run, delimited by a trailing `ALIGN:` line."""

    __slots__ = (
        "new_clusters",
        "naive",
        "prior",
        "nodiv",
        "shuffles",
        "aligns",
        "shrouded",
        "raw",
        "step",
        "sample",
        "reported_clusters",
    )

    def __init__(self, step):
        self.new_clusters = 0
        self.naive = 0
        self.prior = 0
        self.nodiv = 0
        self.shuffles = 0
        self.aligns = None
        self.shrouded = None
        self.raw = None
        self.step = step  # self-consistency step label, if known
        self.sample = None  # dada2-rs per-sample index, if known
        self.reported_clusters = None  # authoritative count from rs sample line

    @property
    def final_clusters(self):
        # C0 is the seed cluster; each successful division adds one.
        return self.new_clusters + 1


def parse(lines):
    segments = []
    rs_samples = []  # (iter, sample, clusters) — emitted after each ALIGN
    rs_iters = []  # (iter, err_delta)
    cur_step = None
    cur = Segment(cur_step)

    interleaved = False
    cur_last_cluster = 0  # for monotonicity check within a segment

    for line in lines:
        # Wrapper / annotation lines. NOTE: do not `continue` — a wrapper
        # marker can share a physical line with a trace token. R prints
        # "selfConsist step N" with no trailing newline, so the step's first
        # ", Division (naive): ..." lands on the same line and must still be
        # counted below.
        m = RE_R_STEP.search(line)
        if m:
            cur_step = int(m.group(1))
            cur.step = cur_step
        m = RE_RS_SAMPLE.search(line)
        if m:
            rs_samples.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
        m = RE_RS_ITER.search(line)
        if m:
            rs_iters.append((int(m.group(1)), m.group(2)))

        # Trace tokens (a physical line may carry a New Cluster + its shuffle
        # run + the division that births the *next* cluster).
        m = RE_NEW_CLUSTER.search(line)
        if m:
            cur.new_clusters += 1
            cur.shuffles += len(m.group(2))
            # Within one (single-threaded) dada run, cluster numbers increase
            # monotonically. A repeat or decrease means several samples' traces
            # are interleaved on one stream (multithreaded dada2-rs), so the
            # per-run breakdown can't be trusted — only the totals can.
            c = int(m.group(1))
            if c <= cur_last_cluster:
                interleaved = True
            cur_last_cluster = c
        cur.naive += len(RE_NAIVE.findall(line))
        cur.prior += len(RE_PRIOR.findall(line))
        cur.nodiv += len(RE_NODIV.findall(line))

        m = RE_ALIGN.search(line)
        if m:
            cur.aligns = int(m.group(1))
            cur.shrouded = int(m.group(2))
            cur.raw = int(m.group(3))
            segments.append(cur)
            cur = Segment(cur_step)
            cur_last_cluster = 0

    # If the stream ended mid-run (no trailing ALIGN) but work was recorded,
    # keep the partial segment so nothing is silently dropped.
    if cur.new_clusters or cur.naive or cur.prior or cur.nodiv:
        segments.append(cur)

    # dada2-rs emits one `iter=N sample=S: C cluster(s)` line per dada run, in
    # the same order as the ALIGN segments. The init pass (MAX_CLUST=1) adds
    # leading ALIGN segments with no sample line, so attach to the trailing
    # segments. The reported cluster count is authoritative even when the
    # per-run trace is interleaved.
    if rs_samples and len(rs_samples) <= len(segments):
        offset = len(segments) - len(rs_samples)
        for seg, (it, samp, clusters) in zip(segments[offset:], rs_samples):
            seg.step = it
            seg.sample = samp
            seg.reported_clusters = clusters

    return segments, rs_samples, rs_iters, interleaved


def fmt_int(v):
    return "-" if v is None else f"{v:,}"


def print_table(segments):
    header = [
        "run",
        "step",
        "sample",
        "clusters",
        "naive",
        "prior",
        "nodiv",
        "shuffles",
        "aligns",
        "shrouded",
        "raw",
    ]
    rows = []
    for i, seg in enumerate(segments, 1):
        sample = seg.sample
        rows.append(
            [
                str(i),
                "-" if seg.step is None else str(seg.step),
                "-" if sample is None else str(sample),
                fmt_int(seg.final_clusters),
                fmt_int(seg.naive),
                fmt_int(seg.prior),
                fmt_int(seg.nodiv),
                fmt_int(seg.shuffles),
                fmt_int(seg.aligns),
                fmt_int(seg.shrouded),
                fmt_int(seg.raw),
            ]
        )

    widths = [
        max(len(header[c]), max((len(r[c]) for r in rows), default=0))
        for c in range(len(header))
    ]
    fmt = "  ".join("{:>" + str(w) + "}" for w in widths)
    print(fmt.format(*header))
    print(fmt.format(*["-" * w for w in widths]))
    for r in rows:
        print(fmt.format(*r))


def print_sample_table(segments):
    """Authoritative per-sample cluster counts from dada2-rs `iter sample` lines.

    Reliable even when the trace is interleaved, since these come from
    discrete post-run summary lines rather than the mixed cluster trace.
    """
    rows = [s for s in segments if s.reported_clusters is not None]
    if not rows:
        return
    print()
    print("authoritative per-sample cluster counts (dada2-rs):")
    header = ["iter", "sample", "clusters"]
    body = [[str(s.step), str(s.sample), f"{s.reported_clusters:,}"] for s in rows]
    widths = [max(len(header[c]), max(len(r[c]) for r in body)) for c in range(3)]
    fmt = "  ".join("{:>" + str(w) + "}" for w in widths)
    print("  " + fmt.format(*header))
    for r in body:
        print("  " + fmt.format(*r))


def print_summary(segments, rs_iters, interleaved):
    n_runs = len(segments)
    tot_naive = sum(s.naive for s in segments)
    tot_prior = sum(s.prior for s in segments)
    tot_new = sum(s.new_clusters for s in segments)
    tot_aligns = sum(s.aligns or 0 for s in segments)
    tot_shroud = sum(s.shrouded or 0 for s in segments)

    print()
    print(f"dada runs (ALIGN segments) : {n_runs}")
    print(f"total buds (new clusters)  : {tot_new:,}  (naive {tot_naive:,}, prior {tot_prior:,})")
    print(f"total alignments           : {tot_aligns:,}  ({tot_shroud:,} shrouded)")
    if tot_aligns and n_runs:
        print(f"mean aligns / run          : {tot_aligns / n_runs:,.0f}")
    # The new-cluster/division balance check is only meaningful for a
    # non-interleaved trace; interleaving produces orphaned tokens by design.
    if not interleaved and tot_new != tot_naive + tot_prior:
        print(
            f"  note: new-cluster count ({tot_new}) != naive+prior divisions "
            f"({tot_naive + tot_prior}); trace may be truncated"
        )

    if rs_iters:
        print()
        print("self-consistency convergence (dada2-rs):")
        for it, delta in rs_iters:
            print(f"  iter {it}: max |err_in - err_out| = {delta}")


def main():
    ap = argparse.ArgumentParser(
        description="Summarize dada2-rs learn-errors / R learnErrors verbose output."
    )
    ap.add_argument(
        "input",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="verbose log file (default: stdin)",
    )
    ap.add_argument(
        "--no-table",
        action="store_true",
        help="print only the aggregate summary, not the per-run table",
    )
    args = ap.parse_args()

    segments, _rs_samples, rs_iters, interleaved = parse(args.input)

    if not segments:
        print("No dada/ALIGN trace lines found in input.", file=sys.stderr)
        print(
            "Make sure the run used --verbose (dada2-rs) or verbose=2 (R dada2).",
            file=sys.stderr,
        )
        sys.exit(1)

    if interleaved:
        print(
            "WARNING: cluster numbers are non-monotonic — this looks like a\n"
            "multithreaded run with several samples' traces interleaved on one\n"
            "stream. Per-run rows below mix samples and are NOT reliable; trust\n"
            "the totals and the authoritative per-sample counts instead.\n"
        )

    if not args.no_table:
        print_table(segments)
    print_sample_table(segments)
    print_summary(segments, rs_iters, interleaved)


if __name__ == "__main__":
    main()
