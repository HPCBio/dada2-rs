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

mod common;

use std::path::Path;
use std::process::Command;

use common::{BIN, assert_same_outputs, fixture, learn_errors, scratch};

/// Run `dada-pooled` on the two committed fixtures into `out`.
fn run_pooled(out: &Path, prune: bool, threads: &str) {
    let errs = learn_errors(out, threads);
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

/// The pruned and unpruned move passes must produce identical output.
#[test]
fn dirty_cluster_prune_matches_full_scan() {
    let tmp = scratch("prune");
    let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
    std::fs::create_dir_all(&pruned).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&pruned, true, "1");
    run_pooled(&full, false, "1");
    assert_same_outputs(&pruned, &full, "pruned vs full move pass, single-threaded");

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
/// stored in, which is upstream of the move pass — so this exercises the prune
/// against a different (still deterministic) input ordering.
#[test]
fn dirty_cluster_prune_matches_full_scan_threaded() {
    let tmp = scratch("prune_mt");
    let (pruned, full) = (tmp.join("pruned"), tmp.join("full"));
    std::fs::create_dir_all(&pruned).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&pruned, true, "4");
    run_pooled(&full, false, "4");
    assert_same_outputs(&pruned, &full, "pruned vs full move pass, 4 threads");

    let _ = std::fs::remove_dir_all(&tmp);
}
