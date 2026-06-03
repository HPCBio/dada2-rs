#!/usr/bin/env python3
"""bench_pooled.py — head-to-head pooled-denoising benchmark: R DADA2 vs dada2-rs.

Focuses on POOLED denoising (R `pool=TRUE` / dada2-rs `dada-pooled`), the
worst case for runtime and memory, on real data. Runs the full pipeline for
each stack, timing every step and capturing per-process peak RSS, then prints a
symmetric side-by-side table.

Two platforms:
  illumina : paired-end. filter -> learn(F,R) -> dada-pooled(F,R) ->
             merge-pairs -> make-sequence-table -> remove-bimera-denovo
  pacbio   : single-end long reads. filter -> learn(pacbio errfun) ->
             dada-pooled(band 32, homo-gap -1, k=7) -> make-sequence-table ->
             remove-bimera-denovo

PER-STEP timing AND peak RSS for BOTH stacks. Every step runs as its own
process and is wrapped with os.wait4(): ru_maxrss is the kernel's peak
resident-set high-water mark for that child. No /usr/bin/time, no /proc
polling — portable to the cluster and macOS. (Linux reports ru_maxrss in kB;
macOS in bytes; we normalize to kB.)

The R side is run step-by-step (each its own `Rscript bench_step.R` process,
state passed via .rds files) precisely so its per-step RSS is comparable to the
Rust side. Caveat: each R step reloads R + dada2 (~150-200 MB baseline), so
small R steps floor at that baseline — the honest cost of invoking R per step.

The dada2-rs binary must be named explicitly with --dada2rs (no auto-discovery):
the build target materially affects the numbers. Use target/release for
reproducible, cross-node-comparable results, or target/release-native for
best-case single-machine performance.

Usage:
  # dada2-rs only, Illumina:
  python3 bench_pooled.py illumina /path/to/miseq_raw \\
      --dada2rs target/release/dada2-rs

  # both stacks, PacBio HiFi, 8 threads:
  python3 bench_pooled.py pacbio /path/to/hifi_raw --run-r --threads 8 \\
      --dada2rs target/release/dada2-rs

Run `python3 bench_pooled.py --help` for all options.
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent   # comparison/benchmark/ -> repo root

R_STEPS = {
    "illumina": ["filter", "learn_fwd", "learn_rev", "dada_fwd", "dada_rev",
                 "merge", "make_table", "remove_bimera"],
    "pacbio": ["filter", "learn", "dada", "make_table", "remove_bimera"],
}


def find_binary(explicit):
    # No auto-discovery: the build target (release vs release-native) materially
    # affects the numbers, so it must be chosen explicitly rather than guessed.
    if not explicit:
        sys.exit("--dada2rs PATH is required (no auto-discovery). Choose the build\n"
                 "explicitly, e.g. target/release/dada2-rs (reproducible) or\n"
                 "target/release-native/dada2-rs (best-case, machine-specific).")
    if not Path(explicit).is_file():
        sys.exit(f"--dada2rs: not a file: {explicit}")
    return explicit


def run_step(name, cmd, logf, results, append_log=False):
    """Run cmd as one process; record (name, wall_s, maxrss_kb, rc). Returns rc."""
    print(f"  ==> {name}: {' '.join(str(c) for c in cmd)}", flush=True)
    start = time.time()
    with open(logf, "ab" if append_log else "wb") as lf:
        proc = subprocess.Popen([str(c) for c in cmd],
                                stdout=lf if append_log else subprocess.DEVNULL,
                                stderr=subprocess.STDOUT if append_log else lf)
        _pid, status, rusage = os.wait4(proc.pid, 0)
    wall = time.time() - start
    rc = os.waitstatus_to_exitcode(status)
    maxrss = rusage.ru_maxrss
    maxrss_kb = maxrss / 1024 if sys.platform == "darwin" else maxrss
    results.append({"step": name, "wall_s": wall, "maxrss_kb": maxrss_kb, "rc": rc})
    if rc != 0:
        print(f"      FAILED (rc={rc}); see {logf}", file=sys.stderr)
    return rc


# --------------------------------------------------------------------------
# dada2-rs pipelines
# --------------------------------------------------------------------------
def rust_illumina(args, bin_, outdir, results):
    filt = outdir / "filtered"
    filt.mkdir(parents=True, exist_ok=True)
    fwd = sorted(f for f in glob.glob(str(Path(args.input) / f"*{args.fwd_pattern}*"))
                 if re.search(r"\.f(ast)?q(\.gz)?$", f))
    if not fwd:
        sys.exit(f"no forward reads matching *{args.fwd_pattern}* in {args.input}")
    filtFs, filtRs = [], []
    tl = args.trunc_len.split(",")
    ee = args.max_ee.split(",")
    for f in fwd:
        r = f.replace(args.fwd_pattern, args.rev_pattern)
        if not os.path.exists(r):
            sys.exit(f"missing reverse mate for {f}")
        name = os.path.basename(f).split(args.fwd_pattern)[0]
        ff = filt / f"{name}_F_filt.fastq.gz"
        fr = filt / f"{name}_R_filt.fastq.gz"
        filtFs.append(ff); filtRs.append(fr)
        run_step(f"filter:{name}", [bin_, "filter-and-trim",
                 "--fwd", f, "--filt", ff, "--rev", r, "--filt-rev", fr,
                 "--trunc-len", *tl, "--max-n", str(args.max_n),
                 "--max-ee", *ee, "--trunc-q", str(args.trunc_q),
                 "--compress"], outdir / f"filter_{name}.log", results)

    errF = outdir / "errors_fwd.json"; errR = outdir / "errors_rev.json"
    run_step("learn_fwd", [bin_, "learn-errors", *map(str, filtFs),
             "--nbases", str(int(args.nbases)), "--errfun", "loess",
             "--threads", str(args.threads), "-o", errF],
             outdir / "learn_fwd.log", results)
    run_step("learn_rev", [bin_, "learn-errors", *map(str, filtRs),
             "--nbases", str(int(args.nbases)), "--errfun", "loess",
             "--threads", str(args.threads), "-o", errR],
             outdir / "learn_rev.log", results)

    ddF = outdir / "dada_fwd"; ddR = outdir / "dada_rev"
    ddF.mkdir(exist_ok=True); ddR.mkdir(exist_ok=True)
    run_step("dada_fwd", [bin_, "dada-pooled", *map(str, filtFs),
             "--error-model", errF, "--output-dir", ddF,
             "--threads", str(args.threads)], outdir / "dada_fwd.log", results)
    run_step("dada_rev", [bin_, "dada-pooled", *map(str, filtRs),
             "--error-model", errR, "--output-dir", ddR,
             "--threads", str(args.threads)], outdir / "dada_rev.log", results)

    merged = outdir / "merged.json"
    run_step("merge", [bin_, "merge-pairs",
             "--fwd-dada", *sorted(map(str, ddF.glob("*.json"))),
             "--rev-dada", *sorted(map(str, ddR.glob("*.json"))),
             "--fwd-fastq", *sorted(map(str, filtFs)),
             "--rev-fastq", *sorted(map(str, filtRs)),
             "-o", merged], outdir / "merge.log", results)

    seqtab = outdir / "seqtab.json"
    run_step("make_table", [bin_, "make-sequence-table", merged, "-o", seqtab],
             outdir / "make_table.log", results)
    run_step("remove_bimera", [bin_, "remove-bimera-denovo", seqtab,
             "--method", "consensus", "-o", outdir / "seqtab_nochim.json"],
             outdir / "bimera.log", results)


def rust_pacbio(args, bin_, outdir, results):
    filt = outdir / "filtered"
    filt.mkdir(parents=True, exist_ok=True)
    reads = sorted(f for f in glob.glob(str(Path(args.input) / "*"))
                   if re.search(r"\.f(ast)?q(\.gz)?$", f))
    if not reads:
        sys.exit(f"no reads in {args.input}")
    filts = []
    for f in reads:
        name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", os.path.basename(f))
        ff = filt / f"{name}_filt.fastq.gz"
        filts.append(ff)
        run_step(f"filter:{name}", [bin_, "filter-and-trim",
                 "--fwd", f, "--filt", ff,
                 "--min-len", str(int(args.min_len)), "--max-len", str(int(args.max_len)),
                 "--max-n", str(args.max_n), "--max-ee", str(args.max_ee.split(",")[0]),
                 "--trunc-q", str(args.trunc_q), "--compress"],
                 outdir / f"filter_{name}.log", results)

    err = outdir / "errors_pacbio.json"
    run_step("learn", [bin_, "learn-errors", *map(str, filts),
             "--nbases", str(int(args.nbases)), "--errfun", "pacbio",
             "--band", str(args.band), "--threads", str(args.threads), "-o", err],
             outdir / "learn.log", results)

    dd = outdir / "dada"; dd.mkdir(exist_ok=True)
    run_step("dada", [bin_, "dada-pooled", *map(str, filts),
             "--error-model", err, "--output-dir", dd,
             "--band", str(args.band), "--homo-gap-p", str(args.homo_gap),
             "--kmer-size", str(args.kmer_size),
             "--threads", str(args.threads)], outdir / "dada.log", results)

    seqtab = outdir / "seqtab.json"
    run_step("make_table", [bin_, "make-sequence-table",
             *sorted(map(str, dd.glob("*.json"))), "-o", seqtab],
             outdir / "make_table.log", results)
    run_step("remove_bimera", [bin_, "remove-bimera-denovo", seqtab,
             "--method", "consensus", "-o", outdir / "seqtab_nochim.json"],
             outdir / "bimera.log", results)


# --------------------------------------------------------------------------
# R DADA2 pipeline — one Rscript process per step (symmetric per-step RSS)
# --------------------------------------------------------------------------
def r_common_args(args, statedir):
    a = [f"platform={args.platform}", f"statedir={statedir}",
         f"threads={args.threads}", f"nbases={args.nbases}", f"input={args.input}",
         f"max_ee={args.max_ee}", f"trunc_q={args.trunc_q}", f"max_n={args.max_n}"]
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}", f"homo_gap={args.homo_gap}"]
    return a


def run_r_split(args, outdir, results):
    rscript = args.rscript or "Rscript"
    statedir = outdir / "Rstate"; statedir.mkdir(parents=True, exist_ok=True)
    log = outdir / "r_pipeline.log"
    if log.exists():
        log.unlink()
    n_asv = None
    for step in R_STEPS[args.platform]:
        cmd = [rscript, str(HERE / "bench_step.R"), f"step={step}", *r_common_args(args, statedir)]
        rc = run_step(step, cmd, log, results, append_log=True)
        if rc != 0:
            break
    text = log.read_text(errors="replace")
    m = re.search(r"^BENCH_RESULT\tn_asv\t(\d+)", text, re.M)
    if m:
        n_asv = int(m.group(1))
    return n_asv


def run_r_single(args, outdir):
    """Run the whole R pipeline as ONE process (realistic single session).

    Gives a fair end-to-end wall time and a single overall peak RSS (one
    accumulating process, unlike the split mode's per-step reload). Per-step
    wall comes from the R script's system.time() BENCH_STEP lines; per-step RSS
    is not available in this mode (one process spans all steps).
    """
    rscript = args.rscript or "Rscript"
    rdir = outdir / "Rsingle"; rdir.mkdir(parents=True, exist_ok=True)
    a = [f"platform={args.platform}", f"input={args.input}", f"outdir={rdir}",
         f"threads={args.threads}", f"nbases={args.nbases}",
         f"max_ee={args.max_ee}", f"trunc_q={args.trunc_q}", f"max_n={args.max_n}"]
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}", f"homo_gap={args.homo_gap}"]
    log = outdir / "r_single.log"
    if log.exists():
        log.unlink()
    tmp = []
    run_step("R-single (whole pipeline)",
             [rscript, str(HERE / "run_dada2_pooled.R"), *a], log, tmp, append_log=True)
    text = log.read_text(errors="replace")
    steps = [{"step": m.group(1), "wall_s": float(m.group(2))}
             for m in re.finditer(r"^BENCH_STEP\t(\S+)\t([\d.]+)", text, re.M)]
    m = re.search(r"^BENCH_RESULT\tn_asv\t(\d+)", text, re.M)
    return {"total_w": tmp[0]["wall_s"], "peak": tmp[0]["maxrss_kb"],
            "steps": steps, "n_asv": int(m.group(1)) if m else None}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def fmt_rss(kb):
    return f"{kb/1024:.0f} MB" if kb else "—"


def collapse(results):
    """Collapse per-sample filter:* rows into a single 'filter' row; keep order."""
    filt = [r for r in results if r["step"].startswith("filter")]
    rows = []
    if filt:
        rows.append({"step": "filter",
                     "wall_s": sum(r["wall_s"] for r in filt),
                     "maxrss_kb": max(r["maxrss_kb"] for r in filt)})
    rows += [r for r in results if not r["step"].startswith("filter")]
    return rows


def print_stack(label, results, cf):
    rows = collapse(results)
    total_w = sum(r["wall_s"] for r in rows)
    peak = max((r["maxrss_kb"] for r in rows), default=0)
    print(f"\n{label} (per step):")
    print(f"  {'step':<18}{'wall_s':>10}{'peak_rss':>12}")
    for r in rows:
        print(f"  {r['step']:<18}{r['wall_s']:>10.2f}{fmt_rss(r['maxrss_kb']):>12}")
        cf.write(f"{label},{r['step']},{r['wall_s']:.2f},{r['maxrss_kb']:.0f}\n")
    print(f"  {'TOTAL':<18}{total_w:>10.2f}{fmt_rss(peak):>12}")
    cf.write(f"{label},TOTAL,{total_w:.2f},{peak:.0f}\n")
    dada_w = sum(r["wall_s"] for r in rows if r["step"].startswith("dada"))
    return {"total_w": total_w, "peak": peak, "dada_w": dada_w}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("platform", choices=["illumina", "pacbio"])
    p.add_argument("input", help="directory of raw FASTQ files")
    p.add_argument("--outdir", default="bench_pooled_out")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--nbases", type=float, default=1e8)
    p.add_argument("--dada2rs", help="path to dada2-rs binary (REQUIRED; e.g. "
                   "target/release/dada2-rs or target/release-native/dada2-rs)")
    p.add_argument("--rscript", help="path to Rscript (default: Rscript on PATH)")
    p.add_argument("--run-r", action="store_true", help="also run the R DADA2 pipeline")
    p.add_argument("--r-mode", choices=["split", "single", "both"], default="both",
                   help="how to run R: 'split' = one process per step (per-step RSS, "
                        "but each step reloads dada2 so wall is inflated); 'single' = "
                        "one process for the whole pipeline (fair end-to-end wall + one "
                        "overall RSS, per-step wall from system.time, no per-step RSS); "
                        "'both' (default) runs each and reports side by side")
    p.add_argument("--no-run-rust", action="store_true", help="skip the dada2-rs pipeline")
    # Illumina filter params
    p.add_argument("--fwd-pattern", default="_R1")
    p.add_argument("--rev-pattern", default="_R2")
    p.add_argument("--trunc-len", default="240,160")
    # PacBio params
    p.add_argument("--min-len", type=float, default=1000)
    p.add_argument("--max-len", type=float, default=1600)
    p.add_argument("--band", type=int, default=32)
    p.add_argument("--homo-gap", type=int, default=-1)
    p.add_argument("--kmer-size", type=int, default=7)
    # shared filter params
    p.add_argument("--max-ee", default="2,2", help="single value or F,R")
    p.add_argument("--trunc-q", type=int, default=None,
                   help="default: 2 for illumina, 0 for pacbio")
    p.add_argument("--max-n", type=int, default=0)
    args = p.parse_args()

    if args.trunc_q is None:
        args.trunc_q = 2 if args.platform == "illumina" else 0

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rust_results, r_split_results, r_split_nasv, r_single = [], [], None, None
    do_split = args.r_mode in ("split", "both")
    do_single = args.r_mode in ("single", "both")

    if not args.no_run_rust:
        bin_ = find_binary(args.dada2rs)
        print(f"=== dada2-rs ({args.platform}) — {bin_} ===", flush=True)
        rust_out = outdir / "rust"; rust_out.mkdir(exist_ok=True)
        (rust_illumina if args.platform == "illumina" else rust_pacbio)(args, bin_, rust_out, rust_results)

    if args.run_r:
        r_out = outdir / "R"; r_out.mkdir(exist_ok=True)
        if do_split:
            print(f"\n=== R DADA2 ({args.platform}) — split (one process per step) ===", flush=True)
            r_split_nasv = run_r_split(args, r_out, r_split_results)
        if do_single:
            print(f"\n=== R DADA2 ({args.platform}) — single process (whole pipeline) ===", flush=True)
            r_single = run_r_single(args, r_out)

    # ---- report ----
    print("\n" + "=" * 56)
    print(f"POOLED BENCHMARK SUMMARY — {args.platform}, {args.threads} thread(s)")
    print("=" * 56)
    csv_path = outdir / "summary.csv"
    with open(csv_path, "w") as cf:
        cf.write("stack,step,wall_s,maxrss_kb\n")
        rs = print_stack("dada2-rs", rust_results, cf) if rust_results else None
        rr_split = print_stack("R-split", r_split_results, cf) if r_split_results else None

        rr_single = None
        if r_single:
            rr_single = {"total_w": r_single["total_w"], "peak": r_single["peak"]}
            print("\nR-single (whole pipeline in one process):")
            print(f"  {'step':<18}{'wall_s':>10}{'peak_rss':>12}")
            for s in r_single["steps"]:
                print(f"  {s['step']:<18}{s['wall_s']:>10.2f}{'—':>12}")
                cf.write(f"R-single,{s['step']},{s['wall_s']:.2f},\n")
            print(f"  {'TOTAL':<18}{r_single['total_w']:>10.2f}{fmt_rss(r_single['peak']):>12}")
            cf.write(f"R-single,TOTAL,{r_single['total_w']:.2f},{r_single['peak']:.0f}\n")
            print("  (per-step RSS unavailable: one accumulating process spans all steps)")

        nasv = r_split_nasv if r_split_nasv is not None else (r_single or {}).get("n_asv")
        if nasv is not None:
            print(f"\n  R final ASVs (post-chimera): {nasv}")

        # HEADLINE: prefer the single-process R run for a fair end-to-end wall +
        # overall RSS; fall back to split if single wasn't run.
        rr = rr_single or rr_split
        if rs and rr:
            print("\nHEADLINE (R vs dada2-rs):")
            if rr_single:
                print("  [R = single-process run; fair end-to-end wall]")
            else:
                print("  [R = split run; wall inflated by per-step dada2 reloads]")
            # Prefer single-process per-step dada wall (no reload inflation);
            # fall back to the split run's dada steps.
            r_dada = None
            if r_single:
                r_dada = sum(s["wall_s"] for s in r_single["steps"] if s["step"].startswith("dada"))
            elif rr_split:
                r_dada = rr_split["dada_w"]
            if rs["dada_w"] > 0 and r_dada:
                print(f"  pooled denoise : R {r_dada:.1f}s vs dada2-rs "
                      f"{rs['dada_w']:.1f}s  →  {r_dada/rs['dada_w']:.1f}× faster")
            print(f"  end-to-end     : R {rr['total_w']:.1f}s vs dada2-rs "
                  f"{rs['total_w']:.1f}s  →  {rr['total_w']/rs['total_w']:.1f}× faster")
            if rs["peak"] and rr["peak"]:
                print(f"  peak RSS       : R {fmt_rss(rr['peak'])} vs dada2-rs "
                      f"{fmt_rss(rs['peak'])}  →  {rr['peak']/rs['peak']:.1f}× less")
        if rr_split and not rr_single:
            print("\n  (note: each R split step reloads R + dada2 (~150-200 MB baseline);"
                  " small R steps floor there. Use --r-mode both for a fair wall.)")
    print(f"\nWrote {csv_path}")


if __name__ == "__main__":
    main()
