#!/usr/bin/env python3
"""bench_pooled.py — head-to-head denoising benchmark: R DADA2 vs dada2-rs.

Runs the full pipeline for each stack, timing every step and capturing
per-process peak RSS, then prints a symmetric side-by-side table.

Two denoising modes via --pool:
  true  (default): POOLED — R `pool=TRUE` / dada2-rs `dada-pooled`. One giant
                   inference over all samples; the worst case for runtime/memory.
  false          : PER-SAMPLE — R `pool=FALSE` (multithread across samples) /
                   dada2-rs per-sample `dada` fanned across workers. The regime
                   where independent small inferences favor dada2-rs most.
  pseudo         : PSEUDO-POOLING — R `dada(pool="pseudo")` vs dada2-rs
                   `dada-pseudo` (one call each). dada-pseudo runs samples
                   serially while R parallelizes across them, so pseudo
                   wall-time favors R at high thread counts; dada2-rs's edge
                   here is per-sample memory. Thresholds: --pseudo-prevalence
                   / --pseudo-min-abundance (R PSEUDO_PREVALENCE/ABUNDANCE).

Two platforms:
  illumina : paired-end. filter -> learn(F,R) -> dada-pooled(F,R) ->
             merge-pairs -> make-sequence-table -> remove-bimera-denovo
  pacbio   : single-end long reads. remove-primers(+orient+filter) ->
             learn(pacbio errfun, k=7) -> dada-pooled(band 32, homo-gap -1, k=7)
             -> make-sequence-table -> remove-bimera-denovo. (R does primer
             removal and filtering as two steps: removePrimers -> filterAndTrim.)
             Input is RAW, primered reads; pass --primer-fwd/--primer-rev.

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
import concurrent.futures
import glob
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent   # R helper scripts live alongside this file

R_STEPS = {
    "illumina": ["filter", "learn_fwd", "learn_rev", "dada_fwd", "dada_rev",
                 "merge", "make_table", "remove_bimera"],
    # R does primer removal and filtering as two functions (removePrimers ->
    # filterAndTrim); the dada2-rs side consolidates both into remove-primers.
    "pacbio": ["remove_primers", "filter", "learn", "dada", "make_table", "remove_bimera"],
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


def run_phase_concurrent(name, jobs, results, max_workers):
    """Run per-sample subprocesses CONCURRENTLY (up to max_workers at once),
    mirroring R's filterAndTrim/removePrimers multithread=N (one core per sample).

    Records ONE row for the phase: wall_s = wall-clock of the whole batch (not
    the sum of per-sample times), maxrss_kb = max peak RSS of any single child.
    jobs = list of (cmd, logf). Returns the worst rc."""
    print(f"  ==> {name}: {len(jobs)} samples, up to {max_workers} concurrent", flush=True)

    def one(cmd, logf):
        with open(logf, "wb") as lf:
            proc = subprocess.Popen([str(c) for c in cmd],
                                    stdout=subprocess.DEVNULL, stderr=lf)
            _pid, status, rusage = os.wait4(proc.pid, 0)
        rc = os.waitstatus_to_exitcode(status)
        rss = rusage.ru_maxrss / 1024 if sys.platform == "darwin" else rusage.ru_maxrss
        return rc, rss

    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, max_workers)) as ex:
        res = list(ex.map(lambda j: one(*j), jobs))
    wall = time.time() - start
    rc = max((r[0] for r in res), default=0)
    peak = max((r[1] for r in res), default=0)
    results.append({"step": name, "wall_s": wall, "maxrss_kb": peak, "rc": rc})
    if rc != 0:
        print(f"      FAILED ({name}); check logs", file=sys.stderr)
    return rc


# --------------------------------------------------------------------------
# dada2-rs pipelines
# --------------------------------------------------------------------------
def rust_dada_step(args, bin_, step, filts, names, errmodel, ddir, outdir, results,
                   extra=()):
    """Denoise. pool=true -> one dada-pooled call; pool=false -> per-sample dada
    fanned across workers (mirrors R dada(pool=FALSE, multithread=N)). Either way
    ddir ends up with one JSON per sample and the phase records a single row."""
    if args.pool == "true":
        run_step(step, [bin_, "dada-pooled", *map(str, filts),
                 "--error-model", errmodel, "--output-dir", ddir, *extra,
                 "--threads", str(args.threads)], outdir / f"{step}.log", results)
    elif args.pool == "pseudo":
        cmd = [bin_, "dada-pseudo", *map(str, filts),
               "--error-model", errmodel, "--output-dir", ddir, *extra,
               "--pseudo-prevalence", str(args.pseudo_prevalence),
               "--threads", str(args.threads)]
        if args.pseudo_min_abundance is not None:
            cmd += ["--pseudo-min-abundance", str(args.pseudo_min_abundance)]
        run_step(step, cmd, outdir / f"{step}.log", results)
    else:
        jobs = []
        for fp, name in zip(filts, names):
            jobs.append(([bin_, "dada", str(fp), "--error-model", str(errmodel),
                          *extra, "-o", ddir / f"{name}.json"],
                         outdir / f"{step}_{name}.log"))
        run_phase_concurrent(step, jobs, results, args.threads)


def rust_illumina(args, bin_, outdir, results):
    filt = outdir / "filtered"
    filt.mkdir(parents=True, exist_ok=True)
    fwd = sorted(f for f in glob.glob(str(Path(args.input) / f"*{args.fwd_pattern}*"))
                 if re.search(r"\.f(ast)?q(\.gz)?$", f))
    if not fwd:
        sys.exit(f"no forward reads matching *{args.fwd_pattern}* in {args.input}")
    filtFs, filtRs, names, jobs = [], [], [], []
    tl = args.trunc_len.split(",")
    ee = args.max_ee.split(",")
    for f in fwd:
        r = f.replace(args.fwd_pattern, args.rev_pattern)
        if not os.path.exists(r):
            sys.exit(f"missing reverse mate for {f}")
        name = os.path.basename(f).split(args.fwd_pattern)[0]
        ff = filt / f"{name}_F_filt.fastq.gz"
        fr = filt / f"{name}_R_filt.fastq.gz"
        filtFs.append(ff); filtRs.append(fr); names.append(name)
        # filter-and-trim is single-core per sample (--threads only affects bgzf),
        # so we fan samples across workers to match R's filterAndTrim(multithread).
        jobs.append(([bin_, "filter-and-trim",
                      "--fwd", f, "--filt", ff, "--rev", r, "--filt-rev", fr,
                      "--trunc-len", *tl, "--max-n", str(args.max_n),
                      "--max-ee", *ee, "--trunc-q", str(args.trunc_q),
                      "--compress"], outdir / f"filter_{name}.log"))
    run_phase_concurrent("filter", jobs, results, args.threads)

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
    rust_dada_step(args, bin_, "dada_fwd", filtFs, names, errF, ddF, outdir, results)
    rust_dada_step(args, bin_, "dada_rev", filtRs, names, errR, ddR, outdir, results)

    merged = outdir / "merged.json"
    run_step("merge", [bin_, "merge-pairs",
             "--fwd-dada", *sorted(map(str, ddF.glob("*.json"))),
             "--rev-dada", *sorted(map(str, ddR.glob("*.json"))),
             "--fwd-fastq", *sorted(map(str, filtFs)),
             "--rev-fastq", *sorted(map(str, filtRs)),
             "--threads", str(args.threads),
             "-o", merged], outdir / "merge.log", results)

    seqtab = outdir / "seqtab.json"
    run_step("make_table", [bin_, "make-sequence-table", merged, "-o", seqtab],
             outdir / "make_table.log", results)
    run_step("remove_bimera", [bin_, "remove-bimera-denovo", seqtab,
             "--method", "consensus", "--threads", str(args.threads),
             "-o", outdir / "seqtab_nochim.json"],
             outdir / "bimera.log", results)


def rust_pacbio(args, bin_, outdir, results):
    filt = outdir / "filtered"
    filt.mkdir(parents=True, exist_ok=True)
    reads = sorted(f for f in glob.glob(str(Path(args.input) / "*"))
                   if re.search(r"\.f(ast)?q(\.gz)?$", f))
    if not reads:
        sys.exit(f"no reads in {args.input}")
    # Consolidated remove-primers: trims primers, orients, AND applies the same
    # length/quality filters as filter-and-trim, in one pass. This mirrors R's
    # removePrimers() + filterAndTrim() (two functions) as a single dada2-rs step.
    filts, names, jobs = [], [], []
    for f in reads:
        name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", os.path.basename(f))
        ff = filt / f"{name}_filt.fastq.gz"
        filts.append(ff); names.append(name)
        # one sample per worker, mirroring R's removePrimers/filterAndTrim threading
        jobs.append(([bin_, "remove-primers", f,
                      "--fout", ff, "--primer-fwd", args.primer_fwd,
                      "--primer-rev", args.primer_rev, "--max-mismatch", str(args.max_mismatch),
                      "--trim-fwd", "--trim-rev", "--orient",
                      "--min-len", str(int(args.min_len)), "--max-len", str(int(args.max_len)),
                      "--max-n", str(args.max_n), "--max-ee", str(args.max_ee.split(",")[0]),
                      "--trunc-q", str(args.trunc_q), "--compress",
                      "-o", outdir / f"primers_{name}.json"],
                     outdir / f"primers_{name}.log"))
    run_phase_concurrent("remove_primers", jobs, results, args.threads)

    err = outdir / "errors_pacbio.json"
    run_step("learn", [bin_, "learn-errors", *map(str, filts),
             "--nbases", str(int(args.nbases)), "--errfun", "pacbio",
             "--band", str(args.band), "--kmer-size", str(args.kmer_size),
             "--threads", str(args.threads), "-o", err],
             outdir / "learn.log", results)

    dd = outdir / "dada"; dd.mkdir(exist_ok=True)
    rust_dada_step(args, bin_, "dada", filts, names, err, dd, outdir, results,
                   extra=["--band", str(args.band), "--homo-gap-p", str(args.homo_gap),
                          "--kmer-size", str(args.kmer_size)])

    seqtab = outdir / "seqtab.json"
    run_step("make_table", [bin_, "make-sequence-table",
             *sorted(map(str, dd.glob("*.json"))), "-o", seqtab],
             outdir / "make_table.log", results)
    run_step("remove_bimera", [bin_, "remove-bimera-denovo", seqtab,
             "--method", "consensus", "--threads", str(args.threads),
             "-o", outdir / "seqtab_nochim.json"],
             outdir / "bimera.log", results)


# --------------------------------------------------------------------------
# R DADA2 pipeline — one Rscript process per step (symmetric per-step RSS)
# --------------------------------------------------------------------------
def r_common_args(args, statedir):
    # Illumina maxEE is paired (e.g. "2,2"); PacBio is a single value.
    mee = args.max_ee if args.platform == "illumina" else args.max_ee.split(",")[0]
    a = [f"platform={args.platform}", f"statedir={statedir}",
         f"threads={args.threads}", f"nbases={args.nbases}", f"input={args.input}",
         f"max_ee={mee}", f"trunc_q={args.trunc_q}", f"max_n={args.max_n}",
         f"pool={args.pool}", f"pseudo_prevalence={args.pseudo_prevalence}"]
    if args.pseudo_min_abundance is not None:
        a.append(f"pseudo_min_abundance={args.pseudo_min_abundance}")
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}", f"homo_gap={args.homo_gap}",
              f"primer_fwd={args.primer_fwd}", f"primer_rev={args.primer_rev}",
              f"max_mismatch={args.max_mismatch}"]
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
    m = re.search(r"BENCH_RESULT\tn_asv\t(\d+)", text, re.M)
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
    mee = args.max_ee if args.platform == "illumina" else args.max_ee.split(",")[0]
    a = [f"platform={args.platform}", f"input={args.input}", f"outdir={rdir}",
         f"threads={args.threads}", f"nbases={args.nbases}",
         f"max_ee={mee}", f"trunc_q={args.trunc_q}", f"max_n={args.max_n}",
         f"pool={args.pool}", f"pseudo_prevalence={args.pseudo_prevalence}"]
    if args.pseudo_min_abundance is not None:
        a.append(f"pseudo_min_abundance={args.pseudo_min_abundance}")
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}", f"homo_gap={args.homo_gap}",
              f"primer_fwd={args.primer_fwd}", f"primer_rev={args.primer_rev}",
              f"max_mismatch={args.max_mismatch}"]
    log = outdir / "r_single.log"
    if log.exists():
        log.unlink()
    tmp = []
    run_step("R-single (whole pipeline)",
             [rscript, str(HERE / "run_dada2_pooled.R"), *a], log, tmp, append_log=True)
    text = log.read_text(errors="replace")
    steps = [{"step": m.group(1), "wall_s": float(m.group(2))}
             for m in re.finditer(r"BENCH_STEP\t(\S+)\t([\d.]+)", text, re.M)]
    m = re.search(r"BENCH_RESULT\tn_asv\t(\d+)", text, re.M)
    return {"total_w": tmp[0]["wall_s"], "peak": tmp[0]["maxrss_kb"],
            "steps": steps, "n_asv": int(m.group(1)) if m else None}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def fmt_rss(kb):
    return f"{kb/1024:.0f} MB" if kb else "—"


def collapse(results):
    """Collapse per-sample 'name:sample' rows into one 'name' row (sum wall, max
    RSS), preserving first-seen order. Non-namespaced steps pass through."""
    grouped, out = {}, []
    for r in results:
        if ":" in r["step"]:
            key = r["step"].split(":", 1)[0]
            if key not in grouped:
                grouped[key] = {"step": key, "wall_s": 0.0, "maxrss_kb": 0.0}
                out.append(grouped[key])
            grouped[key]["wall_s"] += r["wall_s"]
            grouped[key]["maxrss_kb"] = max(grouped[key]["maxrss_kb"], r["maxrss_kb"])
        else:
            out.append(r)
    return out


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
    p.add_argument("--pool", choices=["true", "false", "pseudo"], default="true",
                   help="denoising mode: 'true' = pooled (R pool=TRUE / dada-pooled, "
                        "the worst case); 'false' = per-sample (R pool=FALSE / "
                        "per-sample dada run concurrently); 'pseudo' = pseudo-pooling "
                        "(R dada(pool=\"pseudo\") vs dada2-rs dada-pseudo, one call each). "
                        "Default true. Note: dada-pseudo runs samples serially while R "
                        "parallelizes across them, so pseudo wall-time favors R at high "
                        "thread counts; dada2-rs's win there is per-sample memory.")
    p.add_argument("--pseudo-prevalence", type=int, default=2,
                   help="pseudo: min samples a sequence must appear in to seed round-2 "
                        "priors (R PSEUDO_PREVALENCE). Default 2")
    p.add_argument("--pseudo-min-abundance", type=int, default=None,
                   help="pseudo: min total abundance to seed a prior (R PSEUDO_ABUNDANCE). "
                        "Default off (R's Inf)")
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
    # PacBio primer removal (mirrors removePrimers); defaults are 27F / 1492R
    p.add_argument("--primer-fwd", default="AGRGTTYGATYMTGGCTCAG",
                   help="forward primer, 5'->3' (PacBio; default 27F)")
    p.add_argument("--primer-rev", default="RGYTACCTTGTTACGACTT",
                   help="reverse primer, 5'->3' catalog direction (PacBio; default 1492R)")
    p.add_argument("--max-mismatch", type=int, default=2,
                   help="max mismatches when matching each primer")
    # shared filter params
    p.add_argument("--max-ee", default=None,
                   help="max expected errors. Default: '2,2' for illumina "
                        "(per-direction F,R), '2' for pacbio (single-end)")
    p.add_argument("--trunc-q", type=int, default=None,
                   help="default: 2 for illumina, 0 for pacbio")
    p.add_argument("--max-n", type=int, default=0)
    args = p.parse_args()

    if args.trunc_q is None:
        args.trunc_q = 2 if args.platform == "illumina" else 0
    if args.max_ee is None:
        args.max_ee = "2,2" if args.platform == "illumina" else "2"

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
    mode = {"true": "pooled", "false": "per-sample", "pseudo": "pseudo"}[args.pool]
    print(f"BENCHMARK SUMMARY — {args.platform}, {mode} denoise, {args.threads} thread(s)")
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
                dlabel = {"true": "pooled denoise", "false": "per-sample dada",
                          "pseudo": "pseudo denoise"}[args.pool]
                print(f"  {dlabel:<14} : R {r_dada:.1f}s vs dada2-rs "
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
