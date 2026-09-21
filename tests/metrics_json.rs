//! `--metrics-json` contract tests (issue #162).
//!
//! Three things have to hold, and none of them are obvious from reading the
//! schema:
//!
//! 1. **Measuring never changes results.** The whole point of separating
//!    instrumentation from printing is that you can leave `--metrics-json` on.
//!    If it perturbed inference it would be worse than useless, because the
//!    perturbation would show up as a finding.
//! 2. **Absent is not zero.** A field the run did not measure must be omitted,
//!    never serialized as `0.0`. A consumer that averages a missing split into
//!    its numbers gets a wrong answer with no error.
//! 3. **The optimisation projections are actually populated.** These are the
//!    #132 / #136 / #139 counters. They only fire on a workload with enough
//!    clusters and bud rounds, so a single-sample smoke test says nothing about
//!    them — which is exactly how a schema field ends up permanently zero
//!    without anyone noticing. The pooled fixtures do exercise them.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

const BIN: &str = env!("CARGO_BIN_EXE_dada2-rs");

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("d2rs_metrics_{tag}_{}", std::process::id()));
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
    assert!(
        out.status.success(),
        "learn-errors failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    errs
}

/// Run `dada-pooled` over both fixtures, returning the parsed metrics document.
fn run_pooled(dir: &Path, errs: &Path, out_sub: &str, extra: &[&str]) -> Value {
    let out = dir.join(out_sub);
    let metrics = dir.join(format!("{out_sub}.metrics.json"));
    let mut cmd = Command::new(BIN);
    cmd.args(["dada-pooled", "--threads", "2", "--error-model"])
        .arg(errs)
        .arg("--output-dir")
        .arg(&out)
        .arg("--metrics-json")
        .arg(&metrics)
        .args(extra)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"));
    let res = cmd.output().expect("dada-pooled");
    assert!(
        res.status.success(),
        "dada-pooled failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    serde_json::from_slice(&std::fs::read(&metrics).expect("metrics file")).expect("valid JSON")
}

/// Concatenated per-sample outputs, with the version tag (git hash) stripped so
/// two builds compare equal.
fn outputs_digest(dir: &Path) -> String {
    let mut names: Vec<_> = std::fs::read_dir(dir)
        .expect("read output dir")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".json") || n.ends_with(".json.gz"))
        .collect();
    names.sort();
    assert!(!names.is_empty(), "no outputs produced");
    let mut all = String::new();
    for n in &names {
        let raw = dada2_rs::misc::read_all_maybe_gz(&dir.join(n)).expect("read output");
        let text = String::from_utf8_lossy(&raw).into_owned();
        let stripped = match text.find("\"dada2_rs_version\"") {
            Some(i) => {
                let end = text[i..].find(',').map(|j| i + j).unwrap_or(text.len());
                format!("{}{}", &text[..i], &text[end..])
            }
            None => text,
        };
        all.push_str(&stripped);
    }
    all
}

/// Measuring must not perturb inference, at either level.
#[test]
fn measurement_does_not_change_results() {
    let dir = tmpdir("neutral");
    let errs = learn_errors(&dir);

    // Baseline with no metrics flags at all.
    let plain = dir.join("plain");
    let res = Command::new(BIN)
        .args(["dada-pooled", "--threads", "2", "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(&plain)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("dada-pooled");
    assert!(res.status.success());

    run_pooled(&dir, &errs, "cheap", &[]);
    run_pooled(&dir, &errs, "full", &["--metrics-attribution"]);

    let base = outputs_digest(&plain);
    assert_eq!(
        base,
        outputs_digest(&dir.join("cheap")),
        "--metrics-json changed the ASV output"
    );
    assert_eq!(
        base,
        outputs_digest(&dir.join("full")),
        "--metrics-attribution changed the ASV output"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The cheap level must omit what it did not measure, and keep what it did.
#[test]
fn cheap_level_omits_attribution_but_keeps_the_free_counters() {
    let dir = tmpdir("levels");
    let errs = learn_errors(&dir);

    let cheap = run_pooled(&dir, &errs, "cheap", &[]);
    assert_eq!(cheap["measure_level"], "phases");
    let r = &cheap["runs"][0];

    // Absent, because the per-comparison timers never ran.
    assert!(
        r["compare"].get("split").is_none(),
        "cheap level leaked a compare split"
    );
    assert!(r["compare"].get("busy").is_none());

    // Present, because these cost nothing: the store-loop fold counts them
    // regardless, and the minimizer index probe depends on `screened`.
    assert!(r["compare"]["screened"].as_u64().unwrap() > 0);
    assert!(r["compare"]["aligned"].as_u64().unwrap() > 0);
    assert!(r["compare"]["attribution"]["map"].as_f64().unwrap() > 0.0);
    assert!(r["footprint"]["screen_vector_bytes"].as_u64().unwrap() > 0);

    let full = run_pooled(&dir, &errs, "full", &["--metrics-attribution"]);
    assert_eq!(full["measure_level"], "attribution");
    let split = &full["runs"][0]["compare"]["split"];
    assert!(!split.is_null(), "attribution level dropped the split");
    assert!(split["screen"].as_f64().unwrap() >= 0.0);
    assert!(split["dp_kernel"].as_f64().unwrap() > 0.0);
    assert!(
        full["runs"][0]["compare"]["map_parallel_efficiency"]
            .as_f64()
            .unwrap()
            > 0.0
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The #132 / #136 / #139 projections, plus bud and p-update churn. These are
/// the fields most likely to sit silently at zero: they need a workload with
/// real bud rounds, which a single-sample smoke test does not provide.
#[test]
fn optimisation_projections_are_populated() {
    let dir = tmpdir("projections");
    let errs = learn_errors(&dir);
    let doc = run_pooled(&dir, &errs, "proj", &["--metrics-attribution"]);
    let r = &doc["runs"][0];

    // A pooled run over both fixtures reaches every one of these paths.
    for (group, field) in [
        ("reconcile", "affected"),
        ("reconcile", "rescan_comps"),
        ("move_pruning", "passes"),
        ("move_pruning", "dirty"),
        ("carry_87", "first_reconcile_calls"),
        ("bud", "calls"),
        ("bud", "raws_scanned"),
        ("p_update", "rounds"),
        ("p_update", "repriced"),
    ] {
        let v = r[group][field]
            .as_u64()
            .unwrap_or_else(|| panic!("{group}.{field} missing or not an integer: {}", r[group]));
        assert!(
            v > 0,
            "{group}.{field} is zero -- is it still being collected?"
        );
    }

    // Alignment-work totals, which the sweep reads.
    assert!(r["run"]["nalign"].as_u64().unwrap() > 0);
    assert!(r["run"]["nshroud"].as_u64().unwrap() > 0);
    assert_eq!(r["run"]["multithread"], true);

    // Pooled denoises the merged table once, so exactly one run entry.
    assert_eq!(doc["runs"].as_array().unwrap().len(), 1);
    assert_eq!(r["sample"], "__pooled__");
    // All four pipeline stages are timed for pooled.
    for stage in ["derep", "merge", "dada", "output"] {
        assert!(
            doc["pipeline"][stage].as_f64().unwrap() >= 0.0,
            "pipeline.{stage} missing"
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

/// The index decision is measured, not predicted, and must be recorded so a run
/// that declines the index is not mistaken for one that never probed.
#[test]
fn minimizer_index_decision_is_recorded() {
    let dir = tmpdir("index");
    let errs = learn_errors(&dir);
    let doc = run_pooled(
        &dir,
        &errs,
        "mz",
        &["--screen-backend", "minimizer", "--kdist-cutoff", "0.63"],
    );
    let idx = &doc["runs"][0]["index"];
    assert!(!idx.is_null(), "minimizer run recorded no index decision");
    assert!(idx["probed_clusters"].as_u64().unwrap() > 0);
    assert!(idx["use_index"].is_boolean());
    // Not forced, so the probe decided on its own.
    assert!(idx["forced"].is_null());
    assert_eq!(doc["runs"][0]["run"]["screen_backend"], "minimizer");

    let _ = std::fs::remove_dir_all(&dir);
}

/// The timers must never run with nowhere to write.
#[test]
fn attribution_requires_metrics_json() {
    let out = Command::new(BIN)
        .args(["dada-pooled", "--metrics-attribution", "--error-model", "x"])
        .arg(fixture("sam1F.fastq.gz"))
        .output()
        .expect("spawn");
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("--metrics-json"),
        "expected a clear requirement error, got: {err}"
    );
}

/// `--verbose` carries the run's shape and results; the attribution tables live
/// in `--metrics-json` and nowhere else (#162 item 2).
///
/// The point of the split is that stderr answers "how big is this run and is it
/// configured sanely" while the JSON answers "where did the time go". A
/// multi-line `ns/comp` table creeping back into stderr would undo that, and
/// would also re-create the collision surface #172 closed.
#[test]
fn verbose_carries_run_shape_not_attribution_tables() {
    let dir = tmpdir("quiet");
    let errs = learn_errors(&dir);

    let out = Command::new(BIN)
        .args(["dada-pooled", "--threads", "2", "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(dir.join("out"))
        .arg("--verbose")
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("dada-pooled");
    assert!(out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);

    // Moved to the JSON.
    for gone in [
        "compare attribution (of",
        "compare split (of",
        "map parallel efficiency",
        "shuffle phases (",
        "shuffle scan split",
        "shuffle redundancy",
        "bud redundancy",
        "p-update churn",
        "move pruning (#132)",
        "reconcile incremental (#136)",
        "#87 carry (#139)",
    ] {
        assert!(
            !err.contains(gone),
            "`{gone}` is back in --verbose; it belongs in --metrics-json"
        );
    }

    // Kept: the run's shape and its results.
    for kept in [
        "alignment backend",
        "cpu allocation",
        "tuning gates",
        "resident Raw footprint",
        "phase times",
        "ALIGN:",
    ] {
        assert!(err.contains(kept), "--verbose lost `{kept}`");
    }

    // And it says where the detail went, rather than just dropping it.
    assert!(
        err.contains("--metrics-json"),
        "--verbose should point at --metrics-json for the attribution"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
