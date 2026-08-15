//! Equivalence tests for the dirty-cluster move pass (issue #132).
//!
//! The shuffle's move pass visits only the clusters holding raws whose best
//! cluster changed in the preceding reconcile, instead of every cluster every
//! iteration. That is a pure work reduction: the partition it produces must be
//! **identical** to the full scan's, because the two differ only in which
//! clusters are examined, never in what happens to a raw once examined.
//!
//! This is the invariant #124's pruning arm nearly broke — a byte-identical
//! partition depends on ascending-`ci` order, strict `>`, and a lowest-`ci`
//! tie-break, and a subtle violation there survived both benchmarking and ASV
//! concordance before exact-equality testing caught it. So these tests compare
//! the full per-sample output, not ASV counts.
//!
//! `DADA2RS_SHUFFLE_NO_PRUNE=1` forces the full scan, so both arms run from one
//! binary — which also removes the failure mode where an A/B is built from the
//! wrong checkout and silently measures the same code twice.

use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_dada2-rs");

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixture(name: &str) -> PathBuf {
    manifest_dir().join("tests/fixtures").join(name)
}

/// Strip the version tag, which embeds the git hash and so differs between any
/// two builds without meaning the results differ.
fn normalized(path: &Path) -> String {
    let raw = dada2_rs::misc::read_all_maybe_gz(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let text = String::from_utf8_lossy(&raw).into_owned();
    let re_start = text.find("\"dada2_rs_version\"");
    match re_start {
        Some(i) => {
            let rest = &text[i..];
            let end = rest.find(',').map(|j| i + j).unwrap_or(text.len());
            format!("{}{}", &text[..i], &text[end..])
        }
        None => text,
    }
}

/// Run `dada-pooled` on the two committed fixtures, returning the output dir.
fn run_pooled(out: &Path, prune: bool, threads: &str) {
    let errs = out.join("errs.json");
    let learn = Command::new(BIN)
        .args(["learn-errors", "--threads", threads, "-o"])
        .arg(&errs)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"))
        .output()
        .expect("learn-errors");
    assert!(
        learn.status.success(),
        "learn-errors failed: {}",
        String::from_utf8_lossy(&learn.stderr)
    );

    let mut cmd = Command::new(BIN);
    cmd.args(["dada-pooled", "--threads", threads, "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(out)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"));
    if !prune {
        cmd.env("DADA2RS_SHUFFLE_NO_PRUNE", "1");
    }
    let run = cmd.output().expect("dada-pooled");
    assert!(
        run.status.success(),
        "dada-pooled failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
}

fn assert_same_outputs(a: &Path, b: &Path, label: &str) {
    let mut names: Vec<_> = std::fs::read_dir(a)
        .expect("read output dir")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".json") || n.ends_with(".json.gz"))
        .collect();
    names.sort();
    assert!(!names.is_empty(), "{label}: no outputs produced");
    for n in &names {
        if n == "errs.json" {
            continue;
        }
        assert_eq!(
            normalized(&a.join(n)),
            normalized(&b.join(n)),
            "{label}: {n} differs between pruned and unpruned move pass"
        );
    }
}

/// The pruned and unpruned move passes must produce identical output.
#[test]
fn dirty_cluster_prune_matches_full_scan() {
    let tmp = std::env::temp_dir().join(format!("d2rs_prune_{}", std::process::id()));
    let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
    std::fs::create_dir_all(&pruned).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&pruned, true, "1");
    run_pooled(&full, false, "1");
    assert_same_outputs(&pruned, &full, "single-threaded");

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
/// stored in, which is upstream of the move pass — so this exercises the prune
/// against a different (still deterministic) input ordering.
#[test]
fn dirty_cluster_prune_matches_full_scan_threaded() {
    let tmp = std::env::temp_dir().join(format!("d2rs_prune_mt_{}", std::process::id()));
    let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
    std::fs::create_dir_all(&pruned).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&pruned, true, "4");
    run_pooled(&full, false, "4");
    assert_same_outputs(&pruned, &full, "4 threads");

    let _ = std::fs::remove_dir_all(&tmp);
}
