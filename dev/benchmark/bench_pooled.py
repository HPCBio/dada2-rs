#!/usr/bin/env python3
"""bench_pooled.py — head-to-head denoising benchmark: R DADA2 vs dada2-rs.

Runs the full pipeline for each stack, timing every step and capturing
per-process peak RSS, then prints a symmetric side-by-side table.

Two denoising modes via --pool:
  true  (default): POOLED — R `pool=TRUE` / dada2-rs `dada-pooled`. One giant
                   inference over all samples; the worst case for runtime/memory.
  false          : PER-SAMPLE — R `pool=FALSE` (multithread across samples) /
                   dada2-rs multi-input `dada` (one process, --sample-jobs across
                   samples). The regime where independent inferences favor
                   dada2-rs most.
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
             learn(pacbio errfun, k=7) -> dada-pooled(band 32, k=7; homo-gap
             falls back to --gap-p unless --homo-gap is set)
             -> make-sequence-table -> remove-bimera-denovo. (R does primer
             removal and filtering as two steps: removePrimers -> filterAndTrim.)
             Input is RAW, primered reads; pass --primer-fwd/--primer-rev.

PER-STEP timing, peak RSS, AND effective-core usage for BOTH stacks. Every step
runs as its own process and is wrapped with os.wait4(): ru_maxrss is the
kernel's peak resident-set high-water mark, and ru_utime+ru_stime is CPU time.
The `cores` column = CPU/wall (effective cores used; ideal ≈ threads for an
in-process step), which quantifies thread under-utilization. No /usr/bin/time,
no /proc polling — portable to the cluster and macOS. (Linux reports ru_maxrss
in kB; macOS in bytes; we normalize to kB.)

For a focused scaling study, --thread-sweep N,N,... prepares inputs once
(filter+learn) then runs ONLY the denoise step at each thread count, reporting
wall / cores / speedup / parallel efficiency (dada2-rs only; skips R).
--sample-jobs-sweep N,N,... does the same at FIXED --threads, varying
--sample-jobs (samples-in-flight; pseudo or false) to find the best wall_s.

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
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent   # R helper scripts live alongside this file

# Set from --verbose in main(). When true, dada2-rs steps get --verbose so their
# progress is captured in that step's OWN per-step log (each run_step writes its
# own logf; we do not merge into one big log).
VERBOSE = False

# dada2-rs subcommands that accept --verbose (make-sequence-table does not). The
# benchmark auto-injects --verbose into these when --verbose is given, keyed on
# the subcommand token (cmd[1]); R steps never match, so they're left alone.
RS_VERBOSE_SUBCMDS = {
    "filter-and-trim", "remove-primers", "learn-errors", "dada", "dada-pooled",
    "dada-pseudo", "merge-pairs", "remove-bimera-denovo",
}


def maybe_verbose(cmd):
    """Append --verbose to a dada2-rs subcommand (keyed on cmd[1]) when the global
    --verbose is set, so its progress is written to that step's own log. Leaves R
    commands and make-sequence-table untouched."""
    if VERBOSE and len(cmd) > 1 and str(cmd[1]) in RS_VERBOSE_SUBCMDS:
        return [*cmd, "--verbose"]
    return cmd


# Set from --align-backend in main(). The dada2-rs alignment backend ("nw" or
# "wfa2"); injected into every alignment-using subcommand so a whole-pipeline
# A/B is one flag. R DADA2 has no such switch (always NW), so a wfa2 run is an
# honest "R-NW vs dada2-rs-WFA" comparison.
ALIGN_BACKEND = "nw"

# Set from --wfa-max-edits in main() (or per-iteration by the cap sweep). The WFA
# edit-budget cap (issue #51); None = leave the binary's own default. Only ever
# injected alongside the wfa2 backend.
WFA_MAX_EDITS = None

# dada2-rs subcommands that actually align (and accept --align-backend /
# --wfa-max-edits). Note the chimera step (remove-bimera-denovo) aligns too, on
# its own path, so it's here.
RS_ALIGN_SUBCMDS = {
    "learn-errors", "dada", "dada-pooled", "dada-pseudo", "remove-bimera-denovo",
}


def maybe_align_backend(cmd):
    """Append --align-backend (and --wfa-max-edits, if set) to a dada2-rs
    alignment subcommand (keyed on cmd[1]) when a non-default backend is
    selected. Only injected for wfa2 so default (nw) runs are byte-identical to
    before and work with any binary; R commands and non-aligning steps (filter,
    merge, make-table) are left untouched. The cap rides with wfa2 only — it is
    meaningless for nw, and the binary ignores it there anyway."""
    if ALIGN_BACKEND != "nw" and len(cmd) > 1 and str(cmd[1]) in RS_ALIGN_SUBCMDS:
        cmd = [*cmd, "--align-backend", ALIGN_BACKEND]
        if WFA_MAX_EDITS is not None:
            cmd = [*cmd, "--wfa-max-edits", str(WFA_MAX_EDITS)]
    return cmd

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
    cmd = maybe_align_backend(maybe_verbose(cmd))
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
    cpu_s = rusage.ru_utime + rusage.ru_stime
    results.append({"step": name, "wall_s": wall, "maxrss_kb": maxrss_kb,
                    "cpu_s": cpu_s, "rc": rc})
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
        cmd = maybe_align_backend(maybe_verbose(cmd))
        with open(logf, "wb") as lf:
            proc = subprocess.Popen([str(c) for c in cmd],
                                    stdout=subprocess.DEVNULL, stderr=lf)
            _pid, status, rusage = os.wait4(proc.pid, 0)
        rc = os.waitstatus_to_exitcode(status)
        rss = rusage.ru_maxrss / 1024 if sys.platform == "darwin" else rusage.ru_maxrss
        return rc, rss, rusage.ru_utime + rusage.ru_stime

    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, max_workers)) as ex:
        res = list(ex.map(lambda j: one(*j), jobs))
    wall = time.time() - start
    rc = max((r[0] for r in res), default=0)
    peak = max((r[1] for r in res), default=0)
    # CPU across all concurrent children / phase wall = effective cores used.
    cpu_s = sum(r[2] for r in res)
    results.append({"step": name, "wall_s": wall, "maxrss_kb": peak,
                    "cpu_s": cpu_s, "rc": rc})
    if rc != 0:
        print(f"      FAILED ({name}); check logs", file=sys.stderr)
    return rc


# --------------------------------------------------------------------------
# dada2-rs pipelines
# --------------------------------------------------------------------------
def rust_dada_step(args, bin_, step, filts, names, errmodel, ddir, outdir, results,
                   extra=(), threads=None, sample_jobs="default"):
    """Denoise as a single dada2-rs call, the direct analog of the R function:
    pool=true -> dada-pooled (R pool=TRUE); pool=false -> multi-input dada (R
    pool=FALSE); pool=pseudo -> dada-pseudo (R pool="pseudo"). All three are one
    process using in-process across-sample concurrency (--sample-jobs), so the
    comparison is function-vs-function, not harness-vs-harness. ddir ends up with
    one JSON per sample and the phase records a single row. `threads` overrides
    args.threads (thread-sweep); `sample_jobs` overrides --sample-jobs
    (jobs-sweep); "default" means use args.sample_jobs (None -> let the command
    pick round(threads/4))."""
    t = args.threads if threads is None else threads
    sj = args.sample_jobs if sample_jobs == "default" else sample_jobs
    if args.pool == "true":
        run_step(step, [bin_, "dada-pooled", *map(str, filts),
                 "--error-model", errmodel, "--output-dir", ddir, *extra,
                 "--threads", str(t)], outdir / f"{step}.log", results)
    elif args.pool == "pseudo":
        cmd = [bin_, "dada-pseudo", *map(str, filts),
               "--error-model", errmodel, "--output-dir", ddir, *extra,
               "--pseudo-prevalence", str(args.pseudo_prevalence),
               "--threads", str(t)]
        if args.pseudo_min_abundance is not None:
            cmd += ["--pseudo-min-abundance", str(args.pseudo_min_abundance)]
        if sj is not None:
            cmd += ["--sample-jobs", str(sj)]
        # dada-pseudo streams by default (faster + lighter); --cache-samples
        # opts back into the all-in-memory mode for comparison.
        if args.cache_samples:
            cmd += ["--cache-samples"]
        run_step(step, cmd, outdir / f"{step}.log", results)
    else:
        # Multi-input dada: one process, per-sample (R pool=FALSE), fanning
        # samples in-process via --sample-jobs (mirrors how pseudo is run).
        cmd = [bin_, "dada", *map(str, filts),
               "--error-model", errmodel, "--output-dir", ddir, *extra,
               "--threads", str(t)]
        if sj is not None:
            cmd += ["--sample-jobs", str(sj)]
        run_step(step, cmd, outdir / f"{step}.log", results)


def prepare_illumina(args, bin_, outdir, results):
    """Filter + learn errors. Returns the inputs the denoise step needs."""
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
    return {"filtFs": filtFs, "filtRs": filtRs, "names": names,
            "errF": errF, "errR": errR}


def rust_illumina(args, bin_, outdir, results):
    prep = prepare_illumina(args, bin_, outdir, results)
    filtFs, filtRs, names = prep["filtFs"], prep["filtRs"], prep["names"]

    ddF = outdir / "dada_fwd"; ddR = outdir / "dada_rev"
    ddF.mkdir(exist_ok=True); ddR.mkdir(exist_ok=True)
    rust_dada_step(args, bin_, "dada_fwd", filtFs, names, prep["errF"], ddF, outdir, results)
    rust_dada_step(args, bin_, "dada_rev", filtRs, names, prep["errR"], ddR, outdir, results)

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


def pacbio_dada_extra(args):
    extra = ["--band", str(args.band), "--kmer-size", str(args.kmer_size)]
    # Only pass --homo-gap-p when explicitly set; otherwise dada2-rs falls back to
    # --gap-p (R's HOMOPOLYMER_GAP_PENALTY=NULL -> GAP_PENALTY semantics). HiFi
    # should leave it at default (HOMOPOLYMER_GAP_PENALTY=-1 is DADA2's 454 rec,
    # not PacBio); maintainer confirms:
    # https://github.com/benjjneb/dada2/issues/1663#issuecomment-1359905397
    if args.homo_gap is not None:
        extra += ["--homo-gap-p", str(args.homo_gap)]
    return extra


def prepare_pacbio(args, bin_, outdir, results):
    """Remove-primers (+orient+filter) + learn errors. Returns denoise inputs."""
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
    # Learn with the SAME alignment params used for denoising (band, homo-gap-p,
    # kmer) so the error model matches the dada step — if --homo-gap is set but
    # only passed to dada, the model would be learned at a different homo-gap-p
    # and dada2-rs would warn. When --homo-gap is unset, both learn and dada
    # fall back to --gap-p, staying consistent.
    learn_cmd = [bin_, "learn-errors", *map(str, filts),
                 "--nbases", str(int(args.nbases)), "--errfun", "pacbio",
                 "--band", str(args.band), "--kmer-size", str(args.kmer_size),
                 "--threads", str(args.threads), "-o", err]
    if args.homo_gap is not None:
        learn_cmd += ["--homo-gap-p", str(args.homo_gap)]
    run_step("learn", learn_cmd, outdir / "learn.log", results)
    return {"filts": filts, "names": names, "err": err}


def rust_pacbio(args, bin_, outdir, results):
    prep = prepare_pacbio(args, bin_, outdir, results)
    dd = outdir / "dada"; dd.mkdir(exist_ok=True)
    rust_dada_step(args, bin_, "dada", prep["filts"], prep["names"], prep["err"],
                   dd, outdir, results, extra=pacbio_dada_extra(args))

    seqtab = outdir / "seqtab.json"
    run_step("make_table", [bin_, "make-sequence-table",
             *sorted(map(str, dd.glob("*.json"))), "-o", seqtab],
             outdir / "make_table.log", results)
    run_step("remove_bimera", [bin_, "remove-bimera-denovo", seqtab,
             "--method", "consensus", "--threads", str(args.threads),
             "-o", outdir / "seqtab_nochim.json"],
             outdir / "bimera.log", results)


# --------------------------------------------------------------------------
# Thread-scaling sweep (dada2-rs only)
# --------------------------------------------------------------------------
def _denoise_measure(args, bin_, prep, tdir, threads, sample_jobs="default"):
    """Run only the denoise step(s) into tdir at the given threads/sample_jobs.
    Returns (wall_s, cpu_s, peak_rss_kb): wall/cpu summed across fwd+rev (illumina)
    or the single call; peak_rss_kb is the max single-step peak RSS."""
    tdir.mkdir(exist_ok=True)
    res = []
    if args.platform == "illumina":
        ddF, ddR = tdir / "dada_fwd", tdir / "dada_rev"
        ddF.mkdir(exist_ok=True); ddR.mkdir(exist_ok=True)
        rust_dada_step(args, bin_, "dada_fwd", prep["filtFs"], prep["names"],
                       prep["errF"], ddF, tdir, res, threads=threads, sample_jobs=sample_jobs)
        rust_dada_step(args, bin_, "dada_rev", prep["filtRs"], prep["names"],
                       prep["errR"], ddR, tdir, res, threads=threads, sample_jobs=sample_jobs)
    else:
        dd = tdir / "dada"
        dd.mkdir(exist_ok=True)
        rust_dada_step(args, bin_, "dada", prep["filts"], prep["names"], prep["err"],
                       dd, tdir, res, threads=threads, sample_jobs=sample_jobs,
                       extra=pacbio_dada_extra(args))
    return (sum(r["wall_s"] for r in res),
            sum(r.get("cpu_s", 0.0) for r in res),
            max((r["maxrss_kb"] for r in res), default=0))


def _prepare_for_sweep(args, bin_, outdir):
    prep_dir = outdir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== prepare (filter+learn) once at {args.threads} threads ===", flush=True)
    return (prepare_illumina if args.platform == "illumina" else prepare_pacbio)(
        args, bin_, prep_dir, [])


def thread_sweep(args, bin_, outdir):
    """Prepare inputs ONCE, then run only the DENOISE step at each thread count
    to characterize in-process scaling. dada2-rs only. `cores` = CPU/wall."""
    sweep = [int(t) for t in args.thread_sweep.split(",")]
    prep = _prepare_for_sweep(args, bin_, outdir)
    rows = []
    for t in sweep:
        print(f"\n=== denoise @ {t} thread(s) ===", flush=True)
        # Pin sample_jobs=1 so this isolates single-sample thread scaling
        # (pseudo only; ignored for pooled/per-sample). Use --sample-jobs-sweep
        # to explore samples-in-flight at fixed threads.
        wall, cpu, _rss = _denoise_measure(args, bin_, prep, outdir / f"t{t}", t,
                                           sample_jobs=1)
        rows.append({"threads": t, "wall_s": wall, "cpu_s": cpu,
                     "cores": cpu / wall if wall > 0 else 0.0})

    base_t, base_w = rows[0]["threads"], rows[0]["wall_s"]
    print("\n" + "=" * 60)
    print(f"THREAD SWEEP — {args.platform}, {args.pool} denoise (speedup/eff vs {base_t}t)")
    print("=" * 60)
    print(f"  {'threads':>8}{'wall_s':>10}{'cores':>8}{'speedup':>9}{'efficiency':>12}")
    csv_path = outdir / "sweep.csv"
    with open(csv_path, "w") as cf:
        cf.write("threads,wall_s,cpu_s,cores,speedup,efficiency\n")
        for r in rows:
            speedup = base_w / r["wall_s"] if r["wall_s"] > 0 else 0.0
            ideal = r["threads"] / base_t
            eff = speedup / ideal if ideal > 0 else 0.0
            print(f"  {r['threads']:>8}{r['wall_s']:>10.2f}{r['cores']:>8.1f}"
                  f"{speedup:>8.2f}×{eff*100:>11.0f}%")
            cf.write(f"{r['threads']},{r['wall_s']:.2f},{r['cpu_s']:.2f},"
                     f"{r['cores']:.2f},{speedup:.3f},{eff:.3f}\n")
    print("\n  cores = CPU/wall (ideal ≈ threads for in-process pooled/pseudo).")
    print(f"  Wrote {csv_path}")


def sample_jobs_sweep(args, bin_, outdir):
    """Prepare inputs ONCE, then run the DENOISE step at FIXED --threads for each
    --sample-jobs value, to find the best samples-in-flight. Read the minimum
    wall_s (the optimum) and watch cpu_s climb back up when the sub-pools get too
    small (per-map overhead returns). pseudo or false (both have --sample-jobs);
    dada2-rs only."""
    if args.pool not in ("pseudo", "false"):
        sys.exit("--sample-jobs-sweep only applies to --pool pseudo or false "
                 "(dada-pooled has no --sample-jobs)")
    sweep = [int(j) for j in args.sample_jobs_sweep.split(",")]
    prep = _prepare_for_sweep(args, bin_, outdir)
    rows = []
    for j in sweep:
        tps = max(1, args.threads // j)
        print(f"\n=== denoise @ {j} sample-job(s) × ~{tps} thread(s), "
              f"{args.threads} total ===", flush=True)
        wall, cpu, rss = _denoise_measure(args, bin_, prep, outdir / f"j{j}",
                                          args.threads, sample_jobs=j)
        rows.append({"jobs": j, "wall_s": wall, "cpu_s": cpu, "maxrss_kb": rss,
                     "cores": cpu / wall if wall > 0 else 0.0})

    base_w = rows[0]["wall_s"]
    best = min(rows, key=lambda r: r["wall_s"])
    print("\n" + "=" * 68)
    print(f"SAMPLE-JOBS SWEEP — {args.platform} {args.pool} denoise, {args.threads} threads")
    print("=" * 68)
    print(f"  {'jobs':>6}{'~t/job':>8}{'wall_s':>10}{'cpu_s':>10}{'cores':>8}"
          f"{'peak_rss':>11}{'speedup':>9}")
    csv_path = outdir / "sweep_jobs.csv"
    with open(csv_path, "w") as cf:
        cf.write("sample_jobs,threads_per_job,wall_s,cpu_s,cores,maxrss_kb,speedup\n")
        for r in rows:
            tps = max(1, args.threads // r["jobs"])
            speedup = base_w / r["wall_s"] if r["wall_s"] > 0 else 0.0
            mark = "  <- fastest" if r is best else ""
            print(f"  {r['jobs']:>6}{tps:>8}{r['wall_s']:>10.2f}{r['cpu_s']:>10.1f}"
                  f"{r['cores']:>8.1f}{fmt_rss(r['maxrss_kb']):>11}{speedup:>8.2f}×{mark}")
            cf.write(f"{r['jobs']},{tps},{r['wall_s']:.2f},{r['cpu_s']:.2f},"
                     f"{r['cores']:.2f},{r['maxrss_kb']:.0f},{speedup:.3f}\n")
    print(f"\n  best wall: {best['jobs']} job(s) ({max(1, args.threads // best['jobs'])} "
          f"threads/job) at {best['wall_s']:.1f}s, {fmt_rss(best['maxrss_kb'])} peak")
    print("  (peak_rss is the max single denoise process; more jobs = more "
          "concurrent working sets = higher peak)")
    print(f"  Wrote {csv_path}")


def load_table(nochim_json):
    """Load a make-sequence-table / remove-bimera-denovo JSON into a dict
    {sequence -> {sample -> count}} (per-sample, order-independent). Returns None
    if the file is missing/unreadable."""
    try:
        with open(nochim_json) as fh:
            d = json.load(fh)
        samples, seqs, counts = d["samples"], d["sequences"], d["counts"]
    except (OSError, ValueError, KeyError):
        return None
    # counts[i][j] = sample i, sequence j
    table = {}
    for j, seq in enumerate(seqs):
        table[seq] = {samples[i]: counts[i][j] for i in range(len(samples))}
    return table


def table_equiv(a, b):
    """Compare two {seq -> {sample -> count}} tables. Returns
    (jaccard, exact_frac, identical) where jaccard is over the ASV sequence sets,
    exact_frac is the fraction of shared ASVs whose full per-sample count vector
    matches, and identical is True iff the tables are byte-equal (same ASV set
    AND every count). None-safe: returns (None, None, None) if either is None."""
    if a is None or b is None:
        return (None, None, None)
    sa, sb = set(a), set(b)
    union = sa | sb
    shared = sa & sb
    jaccard = 1.0 if not union else len(shared) / len(union)
    exact = sum(1 for s in shared if a[s] == b[s])
    exact_frac = 1.0 if not shared else exact / len(shared)
    identical = (sa == sb) and all(a[s] == b[s] for s in shared)
    return (jaccard, exact_frac, identical)


def _run_full_pipeline(args, bin_, outdir):
    """Run the whole dada2-rs pipeline once under the current ALIGN_BACKEND /
    WFA_MAX_EDITS globals; return (wall_s, cpu_s, peak_kb, table)."""
    run_pipeline = rust_illumina if args.platform == "illumina" else rust_pacbio
    res = []
    run_pipeline(args, bin_, outdir, res)
    rows_c = collapse(res)
    wall = sum(r["wall_s"] for r in rows_c)
    cpu = sum(r.get("cpu_s", 0.0) for r in rows_c)
    peak = max((r["maxrss_kb"] for r in rows_c), default=0)
    return wall, cpu, peak, load_table(outdir / "seqtab_nochim.json")


def wfa_max_edits_sweep(args, bin_, outdir):
    """Sweep the WFA edit-budget cap (--wfa-max-edits, issue #51) over the FULL
    dada2-rs pipeline, plus an NW reference arm, reporting end-to-end wall /
    cores / peak RSS and TABLE-LEVEL equivalence (ASV sequence set + per-sample
    counts) against both the NW reference (backend equivalence) and the unbounded
    WFA arm (cap correctness-neutrality).

    Unlike the thread / sample-jobs sweeps this re-runs the whole pipeline per
    arm (not just denoise): the cap also governs learn-errors, where the WFA
    O(n·s) cost concentrates on divergent pairs, so a denoise-only sweep would
    miss most of the effect.

    Two questions in one run:
      Q1 (cap)     — WFA(cap) vs WFA(unbounded): should be IDENTICAL (the cap
                     only routes over-budget pairs to NW).
      Q2 (backend) — WFA vs NW: should be ASV-equivalent (jaccard≈1, counts≈)."""
    global ALIGN_BACKEND, WFA_MAX_EDITS
    # 0 = unbounded; accept it as a sweep point alongside positive caps. Sort so
    # the unbounded arm (if requested) is the canonical WFA reference.
    sweep = [int(e) for e in args.wfa_max_edits_sweep.split(",")]

    # --- NW reference arm ---
    print("\n=== full pipeline @ NW (reference) ===", flush=True)
    ALIGN_BACKEND, WFA_MAX_EDITS = "nw", None
    nw_dir = outdir / "nw"
    nw_dir.mkdir(parents=True, exist_ok=True)
    nw_wall, nw_cpu, nw_peak, nw_tab = _run_full_pipeline(args, bin_, nw_dir)
    rows = [{"arm": "nw", "edits": None, "wall_s": nw_wall, "cpu_s": nw_cpu,
             "maxrss_kb": nw_peak, "cores": nw_cpu / nw_wall if nw_wall > 0 else 0.0,
             "table": nw_tab}]

    # --- WFA arms (one per cap) ---
    ALIGN_BACKEND = "wfa2"
    for e in sweep:
        WFA_MAX_EDITS = e
        label = "unbounded" if e == 0 else f"{e} edits"
        print(f"\n=== full pipeline @ wfa2 --wfa-max-edits {e} ({label}) ===", flush=True)
        edir = outdir / (f"e{e}")
        edir.mkdir(parents=True, exist_ok=True)
        wall, cpu, peak, tab = _run_full_pipeline(args, bin_, edir)
        rows.append({"arm": "wfa2", "edits": e, "wall_s": wall, "cpu_s": cpu,
                     "maxrss_kb": peak, "cores": cpu / wall if wall > 0 else 0.0,
                     "table": tab})

    # Canonical WFA reference for the cap-neutrality column: the unbounded arm if
    # present, else the first WFA arm.
    wfa_rows = [r for r in rows if r["arm"] == "wfa2"]
    wfa_ref = next((r for r in wfa_rows if r["edits"] == 0), wfa_rows[0] if wfa_rows else None)

    def nasv(r):
        return "—" if r["table"] is None else str(len(r["table"]))

    print("\n" + "=" * 92)
    print(f"WFA-MAX-EDITS SWEEP — {args.platform} {args.pool}, {args.threads} thread(s), "
          "full pipeline")
    print("=" * 92)
    print(f"  {'arm':>10}{'wall_s':>10}{'cores':>7}{'peak_rss':>11}{'n_asv':>7}"
          f"{'jac/cnt vs NW':>15}{'≡ unbnd WFA':>13}")
    csv_path = outdir / "sweep_wfa_edits.csv"
    with open(csv_path, "w") as cf:
        cf.write("arm,wfa_max_edits,wall_s,cpu_s,cores,maxrss_kb,n_asv,"
                 "jaccard_vs_nw,count_frac_vs_nw,identical_to_unbounded_wfa\n")
        for r in rows:
            arm = "nw" if r["arm"] == "nw" else ("wfa unbnd" if r["edits"] == 0
                                                 else f"wfa {r['edits']}")
            jac, cnt, _ = table_equiv(r["table"], nw_tab)
            vsnw = "—" if jac is None else f"{jac:.3f}/{cnt:.3f}"
            if r["arm"] == "nw" or wfa_ref is None:
                equ = "—"
            else:
                _, _, ident = table_equiv(r["table"], wfa_ref["table"])
                equ = "—" if ident is None else ("yes" if ident else "NO")
            print(f"  {arm:>10}{r['wall_s']:>10.2f}{r['cores']:>7.1f}"
                  f"{fmt_rss(r['maxrss_kb']):>11}{nasv(r):>7}{vsnw:>15}{equ:>13}")
            cf.write(f"{r['arm']},{'' if r['edits'] is None else r['edits']},"
                     f"{r['wall_s']:.2f},{r['cpu_s']:.2f},{r['cores']:.2f},"
                     f"{r['maxrss_kb']:.0f},{nasv(r)},"
                     f"{'' if jac is None else f'{jac:.4f}'},"
                     f"{'' if cnt is None else f'{cnt:.4f}'},"
                     f"{'' if (r['arm'] == 'nw' or wfa_ref is None) else ident}\n")
    print("\n  jac/cnt vs NW : ASV-set Jaccard / per-sample exact-count fraction vs the")
    print("                  NW reference (backend equivalence; expect ≈1.0/≈1.0).")
    print("  ≡ unbnd WFA   : table byte-identical to the unbounded WFA arm — the cap's")
    print("                  correctness-neutrality. 'yes' across all caps = the default")
    print("                  has margin (no real error-copy was truncated). A 'NO' marks")
    print("                  the cap value where an ASV first changed.")
    print(f"  Wrote {csv_path}")


# --------------------------------------------------------------------------
# R DADA2 pipeline — one Rscript process per step (symmetric per-step RSS)
# --------------------------------------------------------------------------
def r_common_args(args, statedir):
    # Illumina maxEE is paired (e.g. "2,2"); PacBio is a single value.
    mee = args.max_ee if args.platform == "illumina" else args.max_ee.split(",")[0]
    a = [f"platform={args.platform}", f"statedir={statedir}",
         f"threads={args.threads}", f"nbases={args.nbases}", f"input={args.input}",
         f"max_ee={mee}", f"trunc_q={args.trunc_q}", f"max_n={args.max_n}",
         f"pool={args.pool}", f"pseudo_prevalence={args.pseudo_prevalence}",
         f"r_derep_mode={args.r_derep_mode}"]
    if args.pseudo_min_abundance is not None:
        a.append(f"pseudo_min_abundance={args.pseudo_min_abundance}")
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}",
              f"primer_fwd={args.primer_fwd}", f"primer_rev={args.primer_rev}",
              f"max_mismatch={args.max_mismatch}"]
        # omit homo_gap when unset -> R's getn(...,NULL) -> HOMOPOLYMER_GAP_PENALTY
        # falls back to GAP_PENALTY.
        if args.homo_gap is not None:
            a.append(f"homo_gap={args.homo_gap}")
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
         f"pool={args.pool}", f"pseudo_prevalence={args.pseudo_prevalence}",
         f"r_derep_mode={args.r_derep_mode}"]
    if args.pseudo_min_abundance is not None:
        a.append(f"pseudo_min_abundance={args.pseudo_min_abundance}")
    if args.platform == "illumina":
        a += [f"fwd_pattern={args.fwd_pattern}", f"rev_pattern={args.rev_pattern}",
              f"trunc_len={args.trunc_len}"]
    else:
        a += [f"min_len={args.min_len}", f"max_len={args.max_len}",
              f"band={args.band}",
              f"primer_fwd={args.primer_fwd}", f"primer_rev={args.primer_rev}",
              f"max_mismatch={args.max_mismatch}"]
        if args.homo_gap is not None:
            a.append(f"homo_gap={args.homo_gap}")
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
            "total_cpu": tmp[0].get("cpu_s", 0.0),
            "steps": steps, "n_asv": int(m.group(1)) if m else None}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def fmt_rss(kb):
    return f"{kb/1024:.0f} MB" if kb else "—"


def cores(r):
    """Effective cores used = CPU time / wall time (1.0 ≈ one core saturated)."""
    return (r.get("cpu_s", 0.0) / r["wall_s"]) if r.get("wall_s", 0) > 0 else 0.0


def collapse(results):
    """Collapse per-sample 'name:sample' rows into one 'name' row (sum wall, sum
    CPU, max RSS), preserving first-seen order. Non-namespaced steps pass through."""
    grouped, out = {}, []
    for r in results:
        if ":" in r["step"]:
            key = r["step"].split(":", 1)[0]
            if key not in grouped:
                grouped[key] = {"step": key, "wall_s": 0.0, "maxrss_kb": 0.0, "cpu_s": 0.0}
                out.append(grouped[key])
            grouped[key]["wall_s"] += r["wall_s"]
            grouped[key]["cpu_s"] += r.get("cpu_s", 0.0)
            grouped[key]["maxrss_kb"] = max(grouped[key]["maxrss_kb"], r["maxrss_kb"])
        else:
            out.append(r)
    return out


def print_stack(label, results, cf):
    rows = collapse(results)
    total_w = sum(r["wall_s"] for r in rows)
    total_cpu = sum(r.get("cpu_s", 0.0) for r in rows)
    peak = max((r["maxrss_kb"] for r in rows), default=0)
    print(f"\n{label} (per step):")
    print(f"  {'step':<18}{'wall_s':>10}{'cores':>8}{'peak_rss':>12}")
    for r in rows:
        print(f"  {r['step']:<18}{r['wall_s']:>10.2f}{cores(r):>8.1f}{fmt_rss(r['maxrss_kb']):>12}")
        cf.write(f"{label},{r['step']},{r['wall_s']:.2f},{r.get('cpu_s', 0.0):.2f},"
                 f"{cores(r):.2f},{r['maxrss_kb']:.0f}\n")
    tot_cores = total_cpu / total_w if total_w > 0 else 0.0
    print(f"  {'TOTAL':<18}{total_w:>10.2f}{tot_cores:>8.1f}{fmt_rss(peak):>12}")
    cf.write(f"{label},TOTAL,{total_w:.2f},{total_cpu:.2f},{tot_cores:.2f},{peak:.0f}\n")
    dada_w = sum(r["wall_s"] for r in rows if r["step"].startswith("dada"))
    dada_cpu = sum(r.get("cpu_s", 0.0) for r in rows if r["step"].startswith("dada"))
    return {"total_w": total_w, "peak": peak, "dada_w": dada_w, "dada_cpu": dada_cpu}


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
                        "multi-input dada, one process with --sample-jobs); 'pseudo' "
                        "= pseudo-pooling "
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
    p.add_argument("--sample-jobs", type=int, default=None,
                   help="pseudo/false: samples denoised concurrently (passed to "
                        "dada-pseudo / multi-input dada --sample-jobs). Default: let "
                        "the command decide (round(threads/4))")
    p.add_argument("--cache-samples", action="store_true",
                   help="pseudo: pass --cache-samples to dada-pseudo to hold all samples'"
                        " uniques in memory across both rounds (the old behavior). "
                        "dada-pseudo now STREAMS by default (re-reads per round; caps peak"
                        " at --sample-jobs samples in flight), which benchmarked faster AND"
                        " lighter — use this flag to measure the cached mode for comparison.")
    p.add_argument("--thread-sweep", default=None, metavar="N,N,...",
                   help="dada2-rs-only scaling study: prepare inputs once (filter+learn), "
                        "then run only the denoise step at each comma-separated thread "
                        "count, reporting wall/cores/speedup/efficiency. Skips R and "
                        "downstream steps. e.g. --thread-sweep 1,2,4,8,16,24")
    p.add_argument("--sample-jobs-sweep", default=None, metavar="N,N,...",
                   help="scaling study at FIXED --threads (--pool pseudo or false): prepare "
                        "once, then run the denoise step at each --sample-jobs value to find "
                        "the best samples-in-flight (min wall_s; watch cpu_s climb when "
                        "sub-pools get too small). e.g. --sample-jobs-sweep 1,2,3,4,6,8")
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
    p.add_argument("--r-derep-mode", choices=["filenames", "objects"], default="filenames",
                   help="how R's dada() receives its input. 'filenames' (default): pass "
                        "file paths so dada() dereplicates on the fly, one sample resident "
                        "at a time (streamed; ~ dada2-rs streaming). 'objects': derepFastq() "
                        "all samples up front and pass the list so dada() holds ALL derep "
                        "objects resident (~ dada2-rs --cache-samples). Use 'objects' for a "
                        "like-for-like preloaded-RSS comparison against --cache-samples; the "
                        "derepFastq cost is counted inside the dada step.")
    p.add_argument("--align-backend", choices=["nw", "wfa2"], default="nw",
                   help="dada2-rs pairwise-alignment backend, threaded into every "
                        "alignment-using step (learn-errors, dada*, remove-bimera-denovo). "
                        "'nw' (default) = Needleman-Wunsch; 'wfa2' = experimental WFA. "
                        "R DADA2 always uses NW, so --align-backend wfa2 --run-r compares "
                        "R-NW vs dada2-rs-WFA. Default nw leaves command lines unchanged.")
    p.add_argument("--wfa-max-edits", type=int, default=None, metavar="N",
                   help="dada2-rs WFA edit-budget cap (issue #51), passed to every "
                        "alignment step when --align-backend wfa2. 0 = unbounded. "
                        "Default: leave the binary's own default (50). Only meaningful "
                        "with wfa2.")
    p.add_argument("--wfa-max-edits-sweep", default=None, metavar="N,N,...",
                   help="run the FULL pipeline once as an NW reference plus once per "
                        "WFA cap value, reporting end-to-end wall / cores / peak RSS "
                        "and TABLE-LEVEL equivalence (ASV set + per-sample counts) vs "
                        "both NW (backend equivalence) and the unbounded WFA arm (cap "
                        "correctness-neutrality). Include 0 for the unbounded WFA "
                        "reference. Watch that '≡ unbnd WFA' stays 'yes' across caps "
                        "while wall/peak drop. e.g. --wfa-max-edits-sweep 0,30,50,80")
    p.add_argument("--no-run-rust", action="store_true", help="skip the dada2-rs pipeline")
    p.add_argument("--verbose", action="store_true",
                   help="pass --verbose to each dada2-rs step (filter/remove-primers/"
                        "learn-errors/dada*/merge-pairs/remove-bimera-denovo) so its "
                        "progress is captured in that step's own per-step log "
                        "(e.g. learn_fwd.log, dada.log). Does not merge logs and does "
                        "not affect R steps; timing impact is negligible.")
    # Illumina filter params
    p.add_argument("--fwd-pattern", default="_R1")
    p.add_argument("--rev-pattern", default="_R2")
    p.add_argument("--trunc-len", default="240,160")
    # PacBio params
    p.add_argument("--min-len", type=float, default=1000)
    p.add_argument("--max-len", type=float, default=1600)
    p.add_argument("--band", type=int, default=32)
    p.add_argument("--homo-gap", type=int, default=None,
                   help="PacBio HOMOPOLYMER_GAP_PENALTY. Default: unset, so both "
                        "stacks fall back to the gap penalty (dada2-rs --gap-p; R "
                        "NULL -> GAP_PENALTY = -8). Pass e.g. -1 only to opt into the "
                        "legacy 454/CCS homopolymer-tolerant value (HiFi reads don't "
                        "need it).")
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

    global VERBOSE, ALIGN_BACKEND, WFA_MAX_EDITS
    VERBOSE = args.verbose
    ALIGN_BACKEND = args.align_backend
    WFA_MAX_EDITS = args.wfa_max_edits

    if args.trunc_q is None:
        args.trunc_q = 2 if args.platform == "illumina" else 0
    if args.max_ee is None:
        args.max_ee = "2,2" if args.platform == "illumina" else "2"

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if args.thread_sweep:
        bin_ = find_binary(args.dada2rs)
        print(f"=== thread sweep: dada2-rs ({args.platform}, {args.pool}) — {bin_} ===",
              flush=True)
        thread_sweep(args, bin_, outdir)
        return

    if args.sample_jobs_sweep:
        bin_ = find_binary(args.dada2rs)
        print(f"=== sample-jobs sweep: dada2-rs ({args.platform}, {args.pool}) — {bin_} ===",
              flush=True)
        sample_jobs_sweep(args, bin_, outdir)
        return

    if args.wfa_max_edits_sweep:
        bin_ = find_binary(args.dada2rs)
        print(f"=== wfa-max-edits sweep: dada2-rs ({args.platform}, {args.pool}) — {bin_} ===",
              flush=True)
        wfa_max_edits_sweep(args, bin_, outdir)
        return

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
    if args.pool == "pseudo":
        mode += " cached" if args.cache_samples else " streaming"
    rmode = f", R derep={args.r_derep_mode}" if args.run_r else ""
    bemode = f", dada2-rs align={args.align_backend}" if args.align_backend != "nw" else ""
    print(f"BENCHMARK SUMMARY — {args.platform}, {mode} denoise, "
          f"{args.threads} thread(s){bemode}{rmode}")
    print("=" * 56)
    print(f"  cores = CPU/wall (effective cores; ideal ≈ {args.threads} for an "
          "in-process step, ≈ min(#samples, threads) for fanned steps)")
    csv_path = outdir / "summary.csv"
    with open(csv_path, "w") as cf:
        cf.write("stack,step,wall_s,cpu_s,cores,maxrss_kb\n")
        rs = print_stack("dada2-rs", rust_results, cf) if rust_results else None
        rr_split = print_stack("R-split", r_split_results, cf) if r_split_results else None

        rr_single = None
        if r_single:
            rr_single = {"total_w": r_single["total_w"], "peak": r_single["peak"]}
            tot_cores = (r_single["total_cpu"] / r_single["total_w"]
                         if r_single["total_w"] > 0 else 0.0)
            print("\nR-single (whole pipeline in one process):")
            print(f"  {'step':<18}{'wall_s':>10}{'cores':>8}{'peak_rss':>12}")
            for s in r_single["steps"]:
                # per-step wall from system.time(); cores/RSS only available process-wide
                print(f"  {s['step']:<18}{s['wall_s']:>10.2f}{'—':>8}{'—':>12}")
                cf.write(f"R-single,{s['step']},{s['wall_s']:.2f},,,\n")
            print(f"  {'TOTAL':<18}{r_single['total_w']:>10.2f}{tot_cores:>8.1f}"
                  f"{fmt_rss(r_single['peak']):>12}")
            cf.write(f"R-single,TOTAL,{r_single['total_w']:.2f},{r_single['total_cpu']:.2f},"
                     f"{tot_cores:.2f},{r_single['peak']:.0f}\n")
            print("  (per-step cores/RSS unavailable: one accumulating process spans all steps)")

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
    print("  Compare runs (peak RSS / wall / cores) with compare_bench.py:")
    print(f"    compare_bench.py --baseline a={csv_path.parent} --compare b=OTHER_OUTDIR")


if __name__ == "__main__":
    main()
