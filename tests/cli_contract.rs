//! The CLI's contract with callers: argument guards, input validation, flag
//! defaults recorded in output, and the stderr error format.

mod common;

use std::path::Path;
use std::process::Command;

use common::{BIN, err_model, fixture, param_i64, run, run_expect_err, scratch};

#[test]
fn dada_input_output_guards() {
    let dir = scratch("guards");
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    // No error model needed: these must fail during argument validation,
    // before any denoising. We point --error-model at a path that exists so
    // the guard (not a missing-file error) is what trips.
    let err = err_model();

    // >1 input with -o is rejected.
    let e = run_expect_err(&[
        "dada",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        dir.join("x.json").to_str().unwrap(),
    ]);
    assert!(e.contains("--output"), "unexpected error: {e}");

    // Single input with --output-dir is rejected.
    let e = run_expect_err(&[
        "dada",
        s1.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--output-dir",
        dir.join("d").to_str().unwrap(),
    ]);
    assert!(e.contains("--output-dir"), "unexpected error: {e}");
}

/// Non-ACGT input must fail with a usable error, not a panic (issue #101).
///
/// Before the fix, an `N` reached `compute_lambda` and aborted a rayon worker
/// with a message naming neither the sample nor the sequence -- after a
/// potentially long dereplication, and repeated once per worker. Non-ACGT input
/// is a user error with a known remedy (DADA2's workflow requires maxN=0), so it
/// belongs with the other input guards.
#[test]
fn dada_rejects_non_acgt_with_a_clear_error() {
    let dir = scratch("non_acgt");
    let err = err_model();
    let read = |id: &str, seq: &str| format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len()));

    // One clean unique plus one carrying an N, both long enough to clear the
    // k-mer-size guard so the ACGT check is what trips.
    let clean = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let withn = "ACGTACGTACGTACNTACGTACGTACGTACGT";
    let mut fq = String::new();
    for i in 0..4 {
        fq.push_str(&read(&format!("clean{i}"), clean));
    }
    for i in 0..3 {
        fq.push_str(&read(&format!("ambig{i}"), withn));
    }
    let fq_path = dir.join("withN.fastq");
    std::fs::write(&fq_path, fq).unwrap();

    let e = run_expect_err(&[
        "dada",
        fq_path.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--threads",
        "1",
        "-o",
        dir.join("out.json").to_str().unwrap(),
    ]);

    // Must be the validation error, and must be actionable: name the base, and
    // point at the fix rather than just refusing.
    assert!(e.contains("non-ACGT"), "unexpected error: {e}");
    assert!(e.contains('N'), "error should name the offending base: {e}");
    assert!(
        e.contains("max-n") || e.contains("filter-and-trim"),
        "error should say how to fix it: {e}"
    );
    // And it must be an error, not a panic escaping a worker thread.
    assert!(
        !e.contains("panicked"),
        "should not panic on non-ACGT input: {e}"
    );
}

/// A quality score above the error model's column count must extend the model and
/// warn, not panic with an index-out-of-bounds (issue #102).
///
/// Mirrors R DADA2 (`dada.R:302-312`), which repeats the last column rather than
/// failing -- a model learned on one run, or on a `--nbases` subsample that missed
/// the top quality bin, would otherwise abort on data R handles. Unlike R the
/// warning is not verbose-gated: extrapolated rates move low-abundance calls, so a
/// run that extrapolates says so.
#[test]
fn dada_extends_error_model_for_out_of_range_quality() {
    let dir = scratch("q_extend");
    let read = |id: &str, seq: &str, q: char| {
        format!("@{id}\n{seq}\n+\n{}\n", q.to_string().repeat(seq.len()))
    };
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";

    // Learn a model from reads capped at Phred 30 ('?'), so err_ncol covers Q0..Q30.
    // Every read carries the same Q, so exactly ONE quality column is populated
    // -- and `interpolate`, the default surface since #205, cannot fit that:
    // its kd-tree has no room to place a vertex away from the lone data point,
    // so the tricube zeroes it and the fit fails (loudly, suggesting
    // `--errfun binned-qual`). `direct` needs only one point. The subject here
    // is quality EXTENSION, not the surface, so pin the setup.
    let learn_fq = dir.join("learn.fastq");
    let mut fq = String::new();
    for i in 0..40 {
        fq.push_str(&read(&format!("l{i}"), seq, '?'));
    }
    std::fs::write(&learn_fq, fq).unwrap();
    let err = dir.join("err_lowq.json");
    run(&[
        "learn-errors",
        learn_fq.to_str().unwrap(),
        "--errfun",
        "loess",
        "--loess-surface",
        "direct",
        "--threads",
        "1",
        "-o",
        err.to_str().unwrap(),
    ]);

    // Denoise reads at Phred 40 ('I') -- above what the model covers.
    let hi_fq = dir.join("hi.fastq");
    let mut fq = String::new();
    for i in 0..6 {
        fq.push_str(&read(&format!("h{i}"), seq, 'I'));
    }
    std::fs::write(&hi_fq, fq).unwrap();
    let out = dir.join("hi.json");

    let res = Command::new(BIN)
        .args([
            "dada",
            hi_fq.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&res.stderr);

    assert!(
        res.status.success(),
        "out-of-range quality should extend the model, not fail: {stderr}"
    );
    assert!(
        !stderr.contains("index out of bounds"),
        "must not panic on out-of-range quality: {stderr}"
    );
    // The extrapolation must be announced without needing --verbose.
    assert!(
        stderr.contains("warning") && stderr.contains("Extending"),
        "extension should warn unconditionally: {stderr}"
    );
    assert!(
        stderr.contains("extrapolated"),
        "warning should say the rates are extrapolated: {stderr}"
    );
    assert!(out.exists(), "no output written");
}

/// CLI failures use a stable, parseable one-line format on stderr:
///
/// ```text
/// dada2-rs: error[<ErrorKind>]: <message>
/// ```
///
/// Pinned because consumers (pipelines, nf-core wrappers) may match on it, and
/// because the default `main() -> io::Result` behaviour it replaced printed the
/// `Debug` of io::Error -- `Error: Custom { kind: Other, error: "..." }` -- which
/// buried the message in a struct dump.
#[test]
fn cli_errors_use_the_documented_format() {
    let dir = scratch("err_format");
    let err = err_model();
    let missing = dir.join("does_not_exist.json");

    let e = run_expect_err(&[
        "dada",
        missing.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        dir.join("out.json").to_str().unwrap(),
    ]);
    let line = e.lines().next().unwrap_or_default();
    assert!(
        line.starts_with("dada2-rs: error[") && line.contains("]: "),
        "error line does not match `dada2-rs: error[Kind]: message`: {line:?}"
    );
    // The kind must be the ErrorKind token, not prose, so it can be matched on.
    assert!(
        line.starts_with("dada2-rs: error[NotFound]: "),
        "expected NotFound for a missing input: {line:?}"
    );
    // The path belongs in the message (WithPath), not just the kind.
    assert!(
        line.contains("does_not_exist.json"),
        "message should name the offending path: {line:?}"
    );
    // No struct dump from the std runtime's Debug formatting.
    assert!(!e.contains("Custom {"), "raw io::Error Debug leaked: {e}");
}

/// `--homo-gap-p`, when unset, must default to `--gap-p` (R's
/// `HOMOPOLYMER_GAP_PENALTY = NULL` semantics); both are recorded in the output
/// `params` block. With neither flag the defaults are -8/-8 (unchanged).
#[test]
fn dada_homo_gap_defaults_to_gap_penalty() {
    let dir = scratch("gap_penalty");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");

    let run_dada = |out: &Path, extra: &[&str]| {
        let mut args = vec![
            "dada",
            s1.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        run(&args);
    };

    // Default: both -8.
    let def = dir.join("def.json");
    run_dada(&def, &[]);
    assert_eq!(param_i64(&def, "gap_p"), -8);
    assert_eq!(param_i64(&def, "homo_gap_p"), -8);

    // --gap-p set, --homo-gap-p unset: homo falls back to gap.
    let g = dir.join("g.json");
    run_dada(&g, &["--gap-p", "-4"]);
    assert_eq!(param_i64(&g, "gap_p"), -4);
    assert_eq!(param_i64(&g, "homo_gap_p"), -4);

    // Both set: independent.
    let gh = dir.join("gh.json");
    run_dada(&gh, &["--gap-p", "-4", "--homo-gap-p", "-1"]);
    assert_eq!(param_i64(&gh, "gap_p"), -4);
    assert_eq!(param_i64(&gh, "homo_gap_p"), -1);

    // Positive penalties are normalized to negative (R dada.R:223-227): a
    // positive --gap-p flips sign and homo falls back to the normalized value;
    // a positive --homo-gap-p flips independently.
    let pos = dir.join("pos.json");
    run_dada(&pos, &["--gap-p", "8"]);
    assert_eq!(param_i64(&pos, "gap_p"), -8);
    assert_eq!(param_i64(&pos, "homo_gap_p"), -8);

    let posh = dir.join("posh.json");
    run_dada(&posh, &["--gap-p", "-4", "--homo-gap-p", "1"]);
    assert_eq!(param_i64(&posh, "gap_p"), -4);
    assert_eq!(param_i64(&posh, "homo_gap_p"), -1);
}
