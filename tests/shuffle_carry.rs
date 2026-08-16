//! Equivalence tests for carrying `compmax` across bud rounds (issue #139).
//!
//! The shuffle normally rebuilds its `compmax`/`emax` map from scratch on every
//! `b_shuffle_converge` call — one full cluster-major pass over every comp in
//! the pool, per bud. `DADA2RS_SHUFFLE_CARRY=1` keeps the map alive across buds
//! instead, letting the first reconcile of the next call repair the parts the
//! bud invalidated (the parent cluster's fallen reads, the new cluster's
//! never-folded comps).
//!
//! That is a pure work reduction *if and only if* the carried map ends up
//! exactly where a rebuild would have put it: same argmax, same lowest-`ci`
//! tie-break. The partition must therefore be **byte-identical** between arms.
//!
//! Why this needs more than a fixture A/B: #136's four seeded mutations all
//! survived end-to-end equivalence on these same fixtures, because a small pool
//! never reaches the states that separate the routes. So the real guard is
//! `DADA2RS_RECONCILE_VERIFY=1`, which asserts after *every* reconcile that
//! `compmax`/`emax` equal a full rescan — turned on here, and the only check
//! that can see a carry-specific divergence the moment it happens rather than
//! if it survives to the output. The carry makes that assertion strictly
//! stronger, because the state it validates now spans bud rounds.

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
    match text.find("\"dada2_rs_version\"") {
        Some(i) => {
            let rest = &text[i..];
            let end = rest.find(',').map(|j| i + j).unwrap_or(text.len());
            format!("{}{}", &text[..i], &text[end..])
        }
        None => text,
    }
}

/// Run `dada-pooled` on the two committed fixtures, with the carry on or off.
fn run_pooled(out: &Path, carry: bool, threads: &str, verify: bool) {
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
    if carry {
        cmd.env("DADA2RS_SHUFFLE_CARRY", "1");
    }
    if verify {
        cmd.env("DADA2RS_RECONCILE_VERIFY", "1");
    }
    let run = cmd.output().expect("dada-pooled");
    assert!(
        run.status.success(),
        "dada-pooled failed (carry={carry}, verify={verify}): {}",
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
            "{label}: {n} differs between carried and rebuilt shuffle"
        );
    }
}

/// The carried and rebuilt maps must produce identical output.
#[test]
fn carried_compmax_matches_rebuild() {
    let tmp = std::env::temp_dir().join(format!("d2rs_carry_{}", std::process::id()));
    let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
    std::fs::create_dir_all(&carried).unwrap();
    std::fs::create_dir_all(&rebuilt).unwrap();

    run_pooled(&carried, true, "1", false);
    run_pooled(&rebuilt, false, "1", false);
    assert_same_outputs(&carried, &rebuilt, "single-threaded");

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Same, multi-threaded. `b_compare`'s parallel map changes the order comps are
/// stored in, which is upstream of the map — so this exercises the carry against
/// a different (still deterministic) input ordering.
#[test]
fn carried_compmax_matches_rebuild_threaded() {
    let tmp = std::env::temp_dir().join(format!("d2rs_carry_mt_{}", std::process::id()));
    let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
    std::fs::create_dir_all(&carried).unwrap();
    std::fs::create_dir_all(&rebuilt).unwrap();

    run_pooled(&carried, true, "4", false);
    run_pooled(&rebuilt, false, "4", false);
    assert_same_outputs(&carried, &rebuilt, "4 threads");

    let _ = std::fs::remove_dir_all(&tmp);
}

/// The carry under the full-rescan invariant.
///
/// This is the test that can actually fail for a carry-specific reason. The two
/// above compare *outputs*, which a divergence only reaches if it survives to
/// the partition; this asserts the carried `compmax`/`emax` equal a full rescan
/// after every single reconcile, including the first one of each bud round —
/// the one where the carry's freshly-relocated work lands.
#[test]
fn carried_compmax_survives_full_rescan_verify() {
    let tmp = std::env::temp_dir().join(format!("d2rs_carry_vfy_{}", std::process::id()));
    let out = tmp.join("carried");
    std::fs::create_dir_all(&out).unwrap();

    run_pooled(&out, true, "1", true);

    let _ = std::fs::remove_dir_all(&tmp);
}

/// The carry must also compose with the #132 move-pass prune disabled.
///
/// The two interact: the carry decides what the reconcile has to repair, and the
/// prune decides which clusters the *next* move pass visits based on what that
/// reconcile marked dirty. A carry bug that under-marks would be masked by the
/// unpruned scan, so both combinations need to agree with the baseline.
#[test]
fn carried_compmax_matches_rebuild_unpruned() {
    let tmp = std::env::temp_dir().join(format!("d2rs_carry_np_{}", std::process::id()));
    let (carried, rebuilt) = (tmp.join("carried"), tmp.join("rebuilt"));
    std::fs::create_dir_all(&carried).unwrap();
    std::fs::create_dir_all(&rebuilt).unwrap();

    for (dir, carry) in [(&carried, true), (&rebuilt, false)] {
        let errs = dir.join("errs.json");
        let learn = Command::new(BIN)
            .args(["learn-errors", "--threads", "1", "-o"])
            .arg(&errs)
            .arg(fixture("sam1F.fastq.gz"))
            .arg(fixture("sam2F.fastq.gz"))
            .output()
            .expect("learn-errors");
        assert!(learn.status.success());

        let mut cmd = Command::new(BIN);
        cmd.args(["dada-pooled", "--threads", "1", "--error-model"])
            .arg(&errs)
            .arg("--output-dir")
            .arg(dir)
            .arg(fixture("sam1F.fastq.gz"))
            .arg(fixture("sam2F.fastq.gz"))
            .env("DADA2RS_SHUFFLE_NO_PRUNE", "1");
        if carry {
            cmd.env("DADA2RS_SHUFFLE_CARRY", "1");
        }
        let run = cmd.output().expect("dada-pooled");
        assert!(
            run.status.success(),
            "dada-pooled failed: {}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    assert_same_outputs(&carried, &rebuilt, "unpruned move pass");

    let _ = std::fs::remove_dir_all(&tmp);
}
