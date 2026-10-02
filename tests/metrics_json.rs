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

mod common;

use std::path::Path;
use std::process::Command;

use serde_json::Value;

use common::{BIN, err_model, fixture, json_outputs, normalized, scratch};

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
    let names = json_outputs(dir);
    assert!(!names.is_empty(), "no outputs produced");
    names.iter().map(|n| normalized(&dir.join(n))).collect()
}

/// Measuring must not perturb inference, at either level.
#[test]
fn measurement_does_not_change_results() {
    let dir = scratch("neutral");
    let errs = err_model();

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
    let dir = scratch("levels");
    let errs = err_model();

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
    let dir = scratch("projections");
    let errs = err_model();
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
    let dir = scratch("index");
    let errs = err_model();
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

    // Union equality with the index is a unit test in `minimizers`; here the
    // union is read from the index, so comparing them would prove nothing.
    let occ = &doc["runs"][0]["screen_occupancy"]["minimizer"];
    assert!(occ["density"].as_f64().unwrap() > 0.0);
    assert!(occ["mean_sharing"].as_f64().unwrap() > 1.0);
    assert!(doc["runs"][0]["screen_occupancy"].get("kmer").is_none());

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
    let dir = scratch("quiet");
    let errs = err_model();

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

/// Pull the number that follows `after` on the first line containing `line`.
fn prose_num(stderr: &str, line: &str, after: &str) -> f64 {
    let l = stderr
        .lines()
        .find(|l| l.contains(line))
        .unwrap_or_else(|| panic!("no `{line}` line in --verbose"));
    let at = l
        .find(after)
        .unwrap_or_else(|| panic!("no `{after}` in `{l}`"));
    let num: String = l[at + after.len()..]
        .trim_start()
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '.')
        .collect();
    num.parse()
        .unwrap_or_else(|_| panic!("no number after `{after}` in `{l}`"))
}

/// Prose rounds to `digits` decimals; the JSON must agree to within that.
fn assert_rounds_to(prose: f64, json: f64, digits: i32) {
    let half = 0.5 * 10f64.powi(-digits) + 1e-9;
    assert!((prose - json).abs() <= half, "prose {prose} vs json {json}");
}

/// The prose-only measurements of #178 have JSON homes carrying the same
/// numbers. Checked field by field: a topic-level check is what nearly lost
/// `shuf_comps_scanned` in #162.
#[test]
fn prose_only_measurements_have_json_homes() {
    let dir = scratch("homes");
    let errs = err_model();
    let metrics = dir.join("homes.metrics.json");
    let res = Command::new(BIN)
        .args([
            "dada-pooled",
            "--threads",
            "2",
            "--verbose",
            "--error-model",
        ])
        .arg(&errs)
        .arg("--output-dir")
        .arg(dir.join("out"))
        .arg("--metrics-json")
        .arg(&metrics)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("dada-pooled");
    assert!(res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    let doc: Value = serde_json::from_slice(&std::fs::read(&metrics).unwrap()).unwrap();
    let f = |v: &Value| v.as_f64().unwrap_or_else(|| panic!("missing field: {v}"));

    let occ = &doc["runs"][0]["screen_occupancy"]["kmer"];
    for (line, after, field, digits) in [
        ("kmer8 fill", "mean", "mean_distinct_kmers", 0),
        ("kmer8 fill", "/", "mean_positional", 0),
        ("kmer8 fill", "(", "positional_pct", 1),
        ("kmer8 fill", "max),", "dense_pct", 1),
        ("kmer8 pooled diversity", ":", "union", 0),
        ("kmer8 pooled diversity", "(", "union_pct", 1),
        ("kmer8 pooled diversity", "sharing", "mean_sharing", 0),
    ] {
        assert_rounds_to(prose_num(&err, line, after), f(&occ[field]), digits);
    }
    assert!(
        doc["runs"][0]["screen_occupancy"]
            .get("minimizer")
            .is_none()
    );

    // Prose truncates kB / 1024 to an integer MB.
    let rss = &doc["pipeline"]["peak_rss_mb"];
    for (line, field) in [
        ("peak RSS after derep+merge", "after_derep_merge"),
        ("peak RSS after merge", "after_merge"),
        ("peak RSS after dada", "after_dada"),
    ] {
        assert_eq!(
            prose_num(&err, line, ":"),
            f(&rss[field]).floor(),
            "{field}"
        );
    }

    let d = &doc["pipeline"]["derep_detail"];
    assert_eq!(d["input_kind"], "fastq");
    assert_eq!(d["samples"], 2);
    assert!(d["bytes"].as_u64().unwrap() > 0);
    // FASTQ is one pass: there is no read / parse split to report.
    assert!(d.get("read").is_none() && d.get("parse").is_none());
    let s = &d["per_sample"];
    assert_rounds_to(
        prose_num(&err, "per-sample", "median"),
        f(&s["median"]) * 1e3,
        0,
    );
    assert!(f(&s["min"]) <= f(&s["median"]) && f(&s["median"]) <= f(&s["max"]));

    // Per-sample derep counts moved here from the unnamed `[derep]` lines.
    assert!(
        !err.contains("[derep]"),
        "pooled still prints unnamed [derep] lines"
    );
    let inputs = doc["pipeline"]["inputs"]
        .as_array()
        .expect("pipeline.inputs");
    let names: Vec<&str> = inputs
        .iter()
        .map(|i| i["sample"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["sam1F", "sam2F"]);
    let reads: u64 = inputs.iter().map(|i| i["reads"].as_u64().unwrap()).sum();
    assert_eq!(reads as f64, prose_num(&err, "merged unique(s),", ","));
    for i in inputs {
        assert!(i["uniques"].as_u64().unwrap() > 0);
        assert!(i["uniques"].as_u64() <= i["reads"].as_u64());
        assert!(f(&i["load_seconds"]) > 0.0);
    }

    let _ = std::fs::remove_dir_all(&dir);
}

/// JSON derep inputs carry the read vs parse split #133 turned on.
#[test]
fn derep_detail_splits_read_and_parse_for_json_inputs() {
    let dir = scratch("derep_json");
    let errs = err_model();
    let mut inputs = Vec::new();
    for s in ["sam1F", "sam2F"] {
        let out = dir.join(format!("{s}.derep.json"));
        let res = Command::new(BIN)
            .arg("derep")
            .arg("-o")
            .arg(&out)
            .arg(fixture(&format!("{s}.fastq.gz")))
            .output()
            .expect("derep");
        assert!(res.status.success());
        inputs.push(out);
    }
    let metrics = dir.join("m.json");
    let res = Command::new(BIN)
        .args(["dada-pooled", "--threads", "2", "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(dir.join("out"))
        .arg("--metrics-json")
        .arg(&metrics)
        .args(&inputs)
        .output()
        .expect("dada-pooled");
    assert!(res.status.success());
    let doc: Value = serde_json::from_slice(&std::fs::read(&metrics).unwrap()).unwrap();
    let d = &doc["pipeline"]["derep_detail"];
    assert_eq!(d["input_kind"], "json");
    assert!(d["read"].as_f64().unwrap() >= 0.0);
    assert!(d["parse"].as_f64().unwrap() > 0.0);
    assert!(d["mb_per_s"].as_f64().unwrap() > 0.0);

    // JSON inputs never printed a `[derep]` line; their counts exist only here.
    let inputs = doc["pipeline"]["inputs"]
        .as_array()
        .expect("pipeline.inputs");
    assert_eq!(inputs.len(), 2);
    assert!(inputs.iter().all(|i| i["reads"].as_u64().unwrap() > 0));

    let _ = std::fs::remove_dir_all(&dir);
}
