#!/usr/bin/env python3
"""Analyse a dada2-rs thread sweep from `--verbose` logs (issue #152).

Answers three questions a thread sweep is run to answer, and refuses to answer
them the wrong way:

1. **Does this pool benefit from more threads?**  `run_dada` per thread count,
   averaged over replicates, with the best count marked.  On soil ITS2 the
   answer is *no* -- 48 -> 96 threads makes it slower -- while soil 16S keeps
   improving, so this is genuinely per-dataset and not a property of the
   hardware.

2. **Where in the run does the benefit (or penalty) land?**  Marginal seconds
   per fixed *cluster* segment, from the progress lines (#150).

   Segments are aligned by **cluster index, never by wall time**.  This is the
   whole reason the script exists: at the same `t=`, a 96-thread run and a
   48-thread run are at different cluster positions and therefore in different
   phases of a workload whose serial fraction collapses over the run.  Comparing
   fixed wall-clock windows across thread counts produced a confidently reported
   conclusion with the *sign of the trend inverted* before this was caught.

3. **Why?**  The map's screen/align split.  The k-mer screen streams k-mer
   vectors (bandwidth-bound); the aligner runs DP (compute-bound).  A
   screen-dominated pool saturates at low thread counts; an align-dominated one
   keeps scaling.  Measured on four arms, the ordering is monotonic:

       16S R1   50.8% screen  -> best at 96 (still scaling there)
       16S R2   67.7% screen  -> best at 64
       ITS2 R2  83.2% screen  -> best at 48
       ITS2 R1  85.8% screen  -> best at 48

   Two cautions learned by getting them wrong.  The screen share is itself
   THREAD-DEPENDENT -- ITS2 R1 reads 80.0% at 24 threads and 85.8% at 48, because
   the bandwidth-bound half degrades faster under contention -- so the predictor
   is only comparable when pinned to one reference count (`--ref-threads`).  And
   a sweep that does not bracket the peak cannot locate it: from 48/64/96 alone,
   ITS2 looked monotonically degrading and was reported as "already past the knee
   below 48".  Adding 24 showed 48 is the peak, beating 24 by 12-17%.

   It predicts *within* a pool as well as across pools -- 16S's two reads differ
   by 17 points of screen share and have different knees -- so this is a property
   of the workload's arithmetic intensity, not of the dataset's name.  The split
   is printed by every verbose run, so it is usable up front.  The thread numbers
   are specific to this machine (EPYC 7713); the ordering is what travels.

Also reports the serial block, which is invariant to thread count -- its *share*
rises as the map gets faster, so map optimisation walks toward the Amdahl wall
rather than away from it.  That is the number to quote when deciding between
bandwidth work and serial work.

Usage:
  dev/analyze_thread_sweep.py tmp/issue-152/full-pooling-novaseq-ITS
  dev/analyze_thread_sweep.py <sweep-dir> --segments 6 --read R1

Expected layout (as produced by the sweep job scripts):
  <sweep-dir>/rep<N>/threads<T>/dada/dada.<READ>.log
"""

import argparse
import collections
import glob
import os
import re
import statistics
import sys

RE_PROGRESS = re.compile(r"progress t=(\d+)s cluster (\d+) ")
RE_RUNDADA = re.compile(r"run_dada=([\d.]+)s")
RE_PHASES = re.compile(
    r"compare=([\d.]+)s \(map=([\d.]+)s parallel, store=([\d.]+)s serial\)\s+"
    r"shuffle=([\d.]+)s\s+bud=([\d.]+)s\s+p_update=([\d.]+)s"
)
RE_MAPEFF = re.compile(r"map parallel efficiency: (\d+)% \(busy=(\d+)s")
RE_SCREEN = re.compile(r"kmer screen\s+([\d.]+)s \(\s*([\d.]+)%\)")
RE_ALIGN = re.compile(r"align total\s+([\d.]+)s \(\s*([\d.]+)%\).*?\(([\d.]+)% passed")


class Run:
    """One log: one (rep, threads, read)."""

    def __init__(self, path):
        self.path = path
        self.progress = []          # (wall_s, clusters)
        self.run_dada = None
        self.phases = None          # compare, map, store, shuffle, bud, pupdate
        self.map_eff = self.busy = None
        self.screen_pct = self.align_pct = self.pass_pct = None
        self._parse()

    def _parse(self):
        with open(self.path, errors="ignore") as fh:
            for line in fh:
                m = RE_PROGRESS.search(line)
                if m:
                    self.progress.append((float(m.group(1)), float(m.group(2))))
                    continue
                m = RE_RUNDADA.search(line)
                if m:
                    self.run_dada = float(m.group(1))
                    continue
                m = RE_PHASES.search(line)
                if m:
                    self.phases = tuple(float(x) for x in m.groups())
                    continue
                m = RE_MAPEFF.search(line)
                if m:
                    self.map_eff, self.busy = int(m.group(1)), float(m.group(2))
                    continue
                m = RE_SCREEN.search(line)
                if m:
                    self.screen_pct = float(m.group(2))
                    continue
                m = RE_ALIGN.search(line)
                if m:
                    self.align_pct, self.pass_pct = float(m.group(2)), float(m.group(3))

    @property
    def serial(self):
        """Serial block: store + shuffle + bud + p_update (everything but the map)."""
        if not self.phases:
            return None
        _, _, store, shuffle, bud, pupd = self.phases
        return store + shuffle + bud + pupd

    def time_to(self, clusters):
        """Wall seconds to reach `clusters`, linearly interpolated between points."""
        prev = (0.0, 0.0)
        for wall, cl in self.progress:
            if cl >= clusters:
                span = cl - prev[1]
                f = (clusters - prev[1]) / span if span else 0.0
                return prev[0] + f * (wall - prev[0])
            prev = (wall, cl)
        return None


def load(sweep_dir):
    runs = collections.defaultdict(list)   # (read, threads) -> [Run]
    pat = os.path.join(sweep_dir, "rep*", "threads*", "dada", "dada.*.log")
    for path in sorted(glob.glob(pat)):
        parts = path.split(os.sep)
        threads = int([p for p in parts if p.startswith("threads")][0][7:])
        read = os.path.basename(path).split(".")[1]
        runs[(read, threads)].append(Run(path))
    if not runs:
        sys.exit(f"no logs matched {pat}")
    return runs


def mean(xs):
    xs = [x for x in xs if x is not None]
    return statistics.mean(xs) if xs else None


def spread(xs):
    """Max-min as a percentage of the mean -- the replicate noise floor.

    Printed next to every averaged figure so a difference smaller than the
    replicate spread is visibly not a difference.
    """
    xs = [x for x in xs if x is not None]
    if len(xs) < 2:
        return None
    m = statistics.mean(xs)
    return (max(xs) - min(xs)) / m * 100 if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir")
    ap.add_argument("--segments", type=int, default=4,
                    help="number of equal cluster segments (default 4)")
    ap.add_argument("--read", help="restrict to one read (R1/R2)")
    ap.add_argument("--ref-threads", type=int, default=48,
                    help="thread count at which to report the map's screen/align "
                         "split (default 48). The split is itself thread-dependent "
                         "-- ITS2 R1 reads 80.0%% screen at 24 threads and 85.8%% at "
                         "48, because the bandwidth-bound half degrades faster under "
                         "contention -- so the predictor is only comparable across "
                         "arms when pinned to one reference count.")
    args = ap.parse_args()

    runs = load(args.sweep_dir)
    reads = sorted({r for r, _ in runs} if not args.read else {args.read})
    threads = sorted({t for _, t in runs})

    for read in reads:
        present = [t for t in threads if (read, t) in runs]
        if not present:
            continue
        print(f"\n{'=' * 78}\n{os.path.basename(args.sweep_dir.rstrip('/'))}  {read}\n{'=' * 78}")

        # --- 1. does it scale? -------------------------------------------
        print("\n-- run_dada by thread count (mean of reps; spread = max-min of reps)")
        print(f"{'threads':>8} {'run_dada':>10} {'spread':>8} {'vs best':>9} {'reps':>5}")
        vals = {t: mean([r.run_dada for r in runs[(read, t)]]) for t in present}
        best = min((v, t) for t, v in vals.items() if v is not None)[1]
        for t in present:
            sp = spread([r.run_dada for r in runs[(read, t)]])
            d = (vals[t] - vals[best]) / vals[best] * 100 if vals[t] else 0
            mark = "  <- best" if t == best else ""
            spc = f"{sp:>7.1f}%" if sp is not None else "      --"
            print(f"{t:>8} {vals[t]:>9.1f}s {spc} {d:>+8.1f}% {len(runs[(read, t)]):>5}{mark}")
        if all(spread([r.run_dada for r in runs[(read, t)]]) is None for t in present):
            print("  NOTE: single replicate per arm -- differences below ~2% are not"
                  " distinguishable from run-to-run variation.")

        # --- 2. serial block ---------------------------------------------
        print("\n-- serial block (store+shuffle+bud+p_update): invariant in seconds,"
              " rising in share")
        print(f"{'threads':>8} {'serial':>10} {'share':>8} {'map':>10} {'busy':>10} {'map eff':>8}")
        for t in present:
            rs = runs[(read, t)]
            s, rd = mean([r.serial for r in rs]), vals[t]
            mp = mean([r.phases[1] for r in rs if r.phases])
            bz, me = mean([r.busy for r in rs]), mean([r.map_eff for r in rs])
            if None in (s, rd, mp):
                continue
            print(f"{t:>8} {s:>9.1f}s {s / rd * 100:>7.1f}% {mp:>9.1f}s "
                  f"{bz:>9.0f} {me:>7.0f}%")
        print("  Amdahl cap on further map work = run_dada / serial.")

        # --- 3. where does it land? --------------------------------------
        maxcl = min(max(cl for _, cl in r.progress)
                    for t in present for r in runs[(read, t)] if r.progress)
        edges = [round(maxcl * i / args.segments) for i in range(args.segments + 1)]
        print(f"\n-- marginal seconds per cluster segment (aligned by CLUSTER, not by t=)")
        head = "".join(f"{t:>9}" for t in present)
        print(f"{'segment':>15}{head}   {present[-1]} vs {present[0]}")
        for a, b in zip(edges, edges[1:]):
            row, ok = {}, True
            for t in present:
                xs = []
                for r in runs[(read, t)]:
                    ta = r.time_to(a) if a else 0.0
                    tb = r.time_to(b)
                    if ta is None or tb is None:
                        ok = False
                    else:
                        xs.append(tb - ta)
                row[t] = mean(xs)
                ok = ok and row[t] is not None
            if ok:
                d = (row[present[-1]] - row[present[0]]) / row[present[0]] * 100
                print(f"{a:>7}-{b:<7}" + "".join(f"{row[t]:>9.1f}" for t in present)
                      + f"   {d:>+7.1f}%")

        # --- 4. why? -----------------------------------------------------
        ref = args.ref_threads if (read, args.ref_threads) in runs else present[0]
        r0 = runs[(read, ref)][0]
        if r0.screen_pct is not None:
            print(f"\n-- map composition at {ref} threads (predicts the knee)")
            if ref != args.ref_threads:
                print(f"   WARNING: {args.ref_threads} threads not in this sweep; "
                      f"using {ref}. The split is thread-dependent, so this is NOT "
                      "comparable to the calibration below.")
            print(f"   k-mer screen {r0.screen_pct:>5.1f}% of busy   (streams k-mer"
                  " vectors: bandwidth-bound)")
            print(f"   alignment    {r0.align_pct:>5.1f}% of busy   (DP kernel:"
                  " compute-bound)")
            print(f"   screen pass  {r0.pass_pct:>5.2f}%")
            # Calibration, not a rule: four arms measured on one machine
            # (EPYC 7713, 2 NUMA domains, 8 CCDs). The ORDERING has held on all
            # four; the thread numbers are hardware-specific and will not travel.
            print("   observed on this machine -- best thread count vs screen share")
            print("   (all measured at 48 threads; 2 reps per arm):")
            print("     16S R1   50.8% screen  -> 96 (still scaling at 96)")
            print("     16S R2   67.7% screen  -> 64")
            print("     ITS2 R2  83.2% screen  -> 48")
            print("     ITS2 R1  85.8% screen  -> 48")
            near = min(
                [(50.8, "16S R1"), (67.7, "16S R2"), (83.2, "ITS2 R2"), (85.8, "ITS2 R1")],
                key=lambda x: abs(x[0] - r0.screen_pct),
            )
            print(f"   => this arm at {r0.screen_pct:.1f}% sits nearest {near[1]}"
                  f" ({near[0]}%)")


if __name__ == "__main__":
    main()
