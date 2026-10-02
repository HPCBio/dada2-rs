//! `DADA2RS_MEMBER_ORDER` (issue #157): the experimental member-order arms.
//!
//! These run the debug binary, so `b_bud_incremental`'s cross-check against a
//! full scan is live on every bud round. A reorder that left the positional bud
//! cache stale would panic there.

mod common;

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;

use common::{BIN, fixture, scratch};

/// Run the binary with `DADA2RS_MEMBER_ORDER` set to `arm` (unset if `None`).
fn run(arm: Option<&str>, args: &[&str]) -> Output {
    let mut cmd = Command::new(BIN);
    cmd.args(args).env_remove("DADA2RS_MEMBER_ORDER");
    if let Some(a) = arm {
        cmd.env("DADA2RS_MEMBER_ORDER", a);
    }
    cmd.output().unwrap()
}

/// Error model learned once, with the gate unset.
fn err_model() -> PathBuf {
    static ERR: OnceLock<PathBuf> = OnceLock::new();
    ERR.get_or_init(|| {
        let err = scratch("member_order_err").join("err.json");
        let out = run(
            None,
            &[
                "learn-errors",
                fixture("sam1F.fastq.gz").to_str().unwrap(),
                "--threads",
                "1",
                "-o",
                err.to_str().unwrap(),
            ],
        );
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        err
    })
    .clone()
}

fn dada(arm: Option<&str>, out: &Path) -> Output {
    run(
        arm,
        &[
            "dada",
            fixture("sam1F.fastq.gz").to_str().unwrap(),
            "--error-model",
            err_model().to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ],
    )
}

#[test]
fn every_arm_runs_with_the_bud_cache_cross_check_live() {
    let dir = scratch("member_order_arms");
    for arm in ["sorted", "shuffle:1", "shuffle:2"] {
        let out = dada(Some(arm), &dir.join(format!("{arm}.json")));
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(out.status.success(), "{arm}: {stderr}");
        assert!(
            stderr.contains("CHANGES RESULTS"),
            "{arm}: the result-changing warning must print"
        );
    }
}

#[test]
fn explicit_insertion_is_byte_identical_to_unset() {
    let dir = scratch("member_order_default");
    let (unset, insertion) = (dir.join("unset.json"), dir.join("insertion.json"));
    assert!(dada(None, &unset).status.success());
    let out = dada(Some("insertion"), &insertion);
    assert!(out.status.success());
    assert!(!String::from_utf8_lossy(&out.stderr).contains("CHANGES RESULTS"));
    assert_eq!(
        std::fs::read(&unset).unwrap(),
        std::fs::read(&insertion).unwrap(),
        "insertion must be the untouched baseline"
    );
}

#[test]
fn an_unrecognised_arm_fails_rather_than_running_the_baseline() {
    let dir = scratch("member_order_bad");
    let out = dada(Some("shufle:1"), &dir.join("x.json"));
    assert!(!out.status.success(), "a mistyped arm must not run");
    assert!(String::from_utf8_lossy(&out.stderr).contains("not recognised"));
}
