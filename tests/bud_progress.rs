//! The bud-round progress record is one atomic line (issue #172).
//!
//! DADA2's divisive progress output was ported from R as a *partial* line built
//! from several `eprint!` calls. That is safe in R, where `dada()` denoises one
//! sample per C++ call. dada2-rs runs up to `--sample-jobs` samples
//! concurrently in one process, where one sample's dangling fragment and
//! another's `eprintln!` land on the same line — and a filter anchored on
//! `^\[dada\]` then discards the collided line *and the message glued to it*.
//!
//! These tests pin the two properties that fix depends on:
//!
//! * a record is emitted as a whole line, never left dangling, so nothing can
//!   be glued to it;
//! * the text stays byte-identical to the R-derived format when one `run_dada`
//!   is in flight, and is tagged only when samples run concurrently, which is
//!   the only case where the line would otherwise be unattributable.

use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_dada2-rs");

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("d2rs_budprog_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn learn_errors(dir: &Path) -> PathBuf {
    let errs = dir.join("errs.json");
    let out = Command::new(BIN)
        .args(["learn-errors", "--threads", "2", "-o"])
        .arg(&errs)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("learn-errors");
    assert!(out.status.success(), "learn-errors failed");
    errs
}

/// Returns stderr from a verbose denoise.
fn denoise_stderr(dir: &Path, errs: &Path, jobs: &str, pooled: bool) -> String {
    let out_dir = dir.join(format!("out{jobs}{}", u8::from(pooled)));
    let mut cmd = Command::new(BIN);
    if pooled {
        cmd.args(["dada-pooled", "--threads", "4"]);
    } else {
        cmd.args(["dada", "--threads", "4", "--sample-jobs", jobs]);
    }
    let out = cmd
        .arg("--error-model")
        .arg(errs)
        .arg("--output-dir")
        .arg(&out_dir)
        .arg("--verbose")
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("denoise");
    assert!(
        out.status.success(),
        "denoise failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// Every progress record must be a complete line: it opens with `New Cluster`
/// and closes with the division that ended the round. A dangling record is
/// exactly what let another sample's message be glued on and filtered away.
#[test]
fn every_record_is_a_complete_line() {
    let dir = tmpdir("complete");
    let errs = learn_errors(&dir);

    for (jobs, pooled) in [("1", false), ("2", false), ("1", true)] {
        let err = denoise_stderr(&dir, &errs, jobs, pooled);
        let records: Vec<&str> = err.lines().filter(|l| l.contains("New Cluster")).collect();
        assert!(
            !records.is_empty(),
            "no progress records at jobs={jobs} pooled={pooled}"
        );
        for r in &records {
            // A record ends with the bud outcome that closed its round. The
            // last round ends with `No Division`, every other with a Division.
            assert!(
                r.contains("Division"),
                "record was left dangling (no bud outcome) at jobs={jobs}: {r}"
            );
            // Nothing else may share the line -- neither another message...
            assert!(
                !r.contains("[dada]"),
                "another message was glued onto a record at jobs={jobs}: {r}"
            );
            // ...nor another record. Two on one line means a flush failed to
            // terminate, which is the dangling-partial-line bug itself.
            assert_eq!(
                r.matches("New Cluster").count(),
                1,
                "two records shared a line at jobs={jobs}: {r}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// Serial and pooled runs keep the R-derived text exactly: no tag, because a
/// single `run_dada` has nothing to disambiguate.
#[test]
fn untagged_when_one_run_is_in_flight() {
    let dir = tmpdir("untagged");
    let errs = learn_errors(&dir);

    for (jobs, pooled) in [("1", false), ("1", true)] {
        let err = denoise_stderr(&dir, &errs, jobs, pooled);
        for line in err.lines().filter(|l| l.contains("New Cluster")) {
            assert!(
                line.starts_with("New Cluster"),
                "serial/pooled record should be untagged, got: {line}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// With samples in flight together, each record carries the sample that
/// produced it — otherwise `New Cluster C5:` is unattributable.
#[test]
fn tagged_when_samples_run_concurrently() {
    let dir = tmpdir("tagged");
    let errs = learn_errors(&dir);
    let err = denoise_stderr(&dir, &errs, "2", false);

    let records: Vec<&str> = err.lines().filter(|l| l.contains("New Cluster")).collect();
    assert!(!records.is_empty(), "no progress records");
    for r in &records {
        assert!(
            r.starts_with("[sam1F]") || r.starts_with("[sam2F]"),
            "concurrent record should carry its sample, got: {r}"
        );
    }
    // Both samples should be represented, or the tag is not actually per-sample.
    assert!(records.iter().any(|r| r.starts_with("[sam1F]")));
    assert!(records.iter().any(|r| r.starts_with("[sam2F]")));

    let _ = std::fs::remove_dir_all(&dir);
}
