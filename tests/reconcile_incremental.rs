//! Equivalence tests for the incremental reconcile (issue #136).
//!
//! The reconcile updates each affected raw's best cluster incrementally --
//! testing only the *changed* candidates against a carried incumbent -- instead
//! of rescanning every affected raw's whole candidate list. It reaches the same
//! answer by a different route, so the equality is not structural the way the
//! move-pass prune's was: it depends on the replacement rule reproducing
//! `best_from_cands` exactly, including the lowest-`ci` tie-break on exact
//! float equality.
//!
//! That tie-break is exactly where #124's pruning arm hid a latent bug, which
//! survived both benchmarking and ASV concordance. So these tests compare full
//! per-sample output, not ASV counts -- and the stronger check is the
//! `debug_assertions` invariant in `b_shuffle_converge`, which compares against
//! a real full rescan every iteration.
//!
//! `DADA2RS_RECONCILE_FULL=1` forces the old full-rescan path, so both arms run
//! from one binary -- which also removes the failure mode where an A/B is built
//! from the wrong checkout and silently measures the same code twice.

mod common;

use std::path::Path;
use std::process::Command;

use common::{BIN, assert_same_outputs, fixture, learn_errors, scratch};

/// Run `dada-pooled` on the two committed fixtures into `out`.
fn run_pooled(out: &Path, incremental: bool, threads: &str) {
    let errs = learn_errors(out, threads);
    let mut cmd = Command::new(BIN);
    cmd.args(["dada-pooled", "--threads", threads, "--error-model"])
        .arg(&errs)
        .arg("--output-dir")
        .arg(out)
        .arg(fixture("sam1F.fastq.gz"))
        .arg(fixture("sam2F.fastq.gz"));
    if !incremental {
        cmd.env("DADA2RS_RECONCILE_FULL", "1");
    }
    let run = cmd.output().expect("dada-pooled");
    assert!(
        run.status.success(),
        "dada-pooled failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
}

/// The incremental and full reconcilees must produce identical output.
#[test]
fn incremental_reconcile_matches_full_rescan() {
    let tmp = scratch("recon");
    let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
    std::fs::create_dir_all(&incremental).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&incremental, true, "1");
    run_pooled(&full, false, "1");
    assert_same_outputs(
        &incremental,
        &full,
        "incremental vs full reconcile, single-threaded",
    );

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
/// stored in, which is upstream of the reconcile -- so this exercises the
/// incremental rule against a different (still deterministic) candidate order,
/// where a tie-break error is more likely to surface.
#[test]
fn incremental_reconcile_matches_full_rescan_threaded() {
    let tmp = scratch("recon_mt");
    let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
    std::fs::create_dir_all(&incremental).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&incremental, true, "4");
    run_pooled(&full, false, "4");
    assert_same_outputs(
        &incremental,
        &full,
        "incremental vs full reconcile, 4 threads",
    );

    let _ = std::fs::remove_dir_all(&tmp);
}
