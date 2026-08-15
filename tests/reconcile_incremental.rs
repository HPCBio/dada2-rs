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
fn run_pooled(out: &Path, incremental: bool, threads: &str) {
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
            "{label}: {n} differs between incremental and full reconcile"
        );
    }
}

/// The incremental and full reconcilees must produce identical output.
#[test]
fn incremental_reconcile_matches_full_rescan() {
    let tmp = std::env::temp_dir().join(format!("d2rs_recon_{}", std::process::id()));
    let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
    std::fs::create_dir_all(&incremental).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&incremental, true, "1");
    run_pooled(&full, false, "1");
    assert_same_outputs(&incremental, &full, "single-threaded");

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
/// stored in, which is upstream of the reconcile -- so this exercises the
/// incremental rule against a different (still deterministic) candidate order,
/// where a tie-break error is more likely to surface.
#[test]
fn incremental_reconcile_matches_full_rescan_threaded() {
    let tmp = std::env::temp_dir().join(format!("d2rs_recon_mt_{}", std::process::id()));
    let (incremental, full) = (tmp.join("incremental"), tmp.join("full"));
    std::fs::create_dir_all(&incremental).unwrap();
    std::fs::create_dir_all(&full).unwrap();

    run_pooled(&incremental, true, "4");
    run_pooled(&full, false, "4");
    assert_same_outputs(&incremental, &full, "4 threads");

    let _ = std::fs::remove_dir_all(&tmp);
}
