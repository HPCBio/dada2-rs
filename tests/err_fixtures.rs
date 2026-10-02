//! The committed error models must still be what `learn-errors` produces.
//!
//! Every other integration test loads `tests/fixtures/errs*.json` instead of
//! learning a model: a debug-build learn on the fixtures takes ~30 s on one
//! thread, and the suite used to run one per test (#103). This is the one
//! place that still learns, so a change to `learn-errors` output fails here,
//! by name, rather than silently leaving the other tests on a stale model.
//!
//! After a deliberate change to the error model, regenerate with
//! `DADA2RS_BLESS=1 cargo test --test err_fixtures` and commit the result.

mod common;

use std::path::Path;
use std::process::Command;

use serde_json::Value;

use common::{BIN, fixture, scratch};

/// Floats may differ by this relative amount. Debug and release builds already
/// disagree in the last ulp (~1.5e-16) on two cells of `errs.json`, and libm
/// can differ across platforms; a real change to the model is far larger.
const REL_TOL: f64 = 1e-12;

/// Paths at which `a` and `b` differ, ignoring the version tag and float noise
/// within [`REL_TOL`].
fn differences(a: &Value, b: &Value, path: &str, out: &mut Vec<String>) {
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            for k in x.keys().chain(y.keys().filter(|k| !x.contains_key(*k))) {
                if k != "dada2_rs_version" {
                    let (u, v) = (
                        x.get(k).unwrap_or(&Value::Null),
                        y.get(k).unwrap_or(&Value::Null),
                    );
                    differences(u, v, &format!("{path}.{k}"), out);
                }
            }
        }
        (Value::Array(x), Value::Array(y)) if x.len() == y.len() => {
            for (i, (u, v)) in x.iter().zip(y).enumerate() {
                differences(u, v, &format!("{path}[{i}]"), out);
            }
        }
        (Value::Number(x), Value::Number(y)) if x.is_f64() || y.is_f64() => {
            let (x, y) = (x.as_f64().unwrap(), y.as_f64().unwrap());
            if (x - y).abs() > REL_TOL * x.abs().max(y.abs()) {
                out.push(format!("{path}: {x} vs {y}"));
            }
        }
        _ if a != b => out.push(format!("{path}: {a} vs {b}")),
        _ => {}
    }
}

fn load(path: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(path).unwrap())
        .unwrap_or_else(|e| panic!("parse {}: {e}", path.display()))
}

/// Learn a model from `inputs` and compare it with (or, when blessing, write
/// it to) the committed fixture `name`.
fn check(name: &str, inputs: &[&str]) {
    let fresh = scratch(&format!("err_fixture_{name}")).join(name);
    let mut cmd = Command::new(BIN);
    // Output is identical at any thread count; 4 keeps the debug build quick.
    cmd.args(["learn-errors", "--threads", "4", "-o"])
        .arg(&fresh);
    for i in inputs {
        cmd.arg(fixture(i));
    }
    // learn-errors denoises internally, so a tuning gate left set in the
    // caller's shell would read as drift.
    for (k, _) in std::env::vars().filter(|(k, _)| k.starts_with("DADA2RS_")) {
        cmd.env_remove(k);
    }
    let out = cmd.output().expect("learn-errors");
    assert!(
        out.status.success(),
        "learn-errors failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let committed = fixture(name);
    if std::env::var_os("DADA2RS_BLESS").is_some() {
        std::fs::copy(&fresh, &committed).unwrap();
        return;
    }
    let mut diffs = Vec::new();
    differences(&load(&committed), &load(&fresh), "", &mut diffs);
    assert!(
        diffs.is_empty(),
        "{name} no longer matches a fresh learn-errors run ({} values differ, first: {}); \
         if the error model changed deliberately, regenerate with \
         `DADA2RS_BLESS=1 cargo test --test err_fixtures` and commit it",
        diffs.len(),
        diffs[..diffs.len().min(3)].join("; ")
    );
}

#[test]
fn committed_error_model_matches_a_fresh_learn() {
    check("errs.json", &["sam1F.fastq.gz", "sam2F.fastq.gz"]);
}

#[test]
fn committed_sam1f_error_model_matches_a_fresh_learn() {
    check("errs_sam1F.json", &["sam1F.fastq.gz"]);
}
