//! Scaffolding shared by the integration-test binaries (issue #103).
//!
//! Each file under `tests/` compiles to its own binary and pulls this in with
//! `mod common;`, so a helper one binary leaves unused is not dead code.
#![allow(dead_code)]

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

pub const BIN: &str = env!("CARGO_BIN_EXE_dada2-rs");

pub fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn fixture(name: &str) -> PathBuf {
    // Fixtures live under tests/ (tracked); the repo's /data dir is gitignored
    // and so is absent on CI.
    manifest_dir().join("tests/fixtures").join(name)
}

/// Per-test scratch dir under the system tmp area; cleaned and recreated.
/// The pid keeps concurrently running test binaries apart, so `tag` need only
/// be unique within one binary.
pub fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dada2rs_{}_{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run the binary; panic with stderr on a non-zero exit.
pub fn run(args: &[&str]) {
    let out = Command::new(BIN)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn {BIN}: {e}"));
    assert!(
        out.status.success(),
        "command failed: dada2-rs {}\n--- stderr ---\n{}",
        args.join(" "),
        String::from_utf8_lossy(&out.stderr),
    );
}

/// Run the binary expecting failure; return stderr.
pub fn run_expect_err(args: &[&str]) -> String {
    let out = Command::new(BIN).args(args).output().unwrap();
    assert!(
        !out.status.success(),
        "expected failure but command succeeded: dada2-rs {}",
        args.join(" "),
    );
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// The error model learned from both forward fixtures, committed so tests do not
/// each spend ~30 s (debug build) learning it; `err_fixtures.rs` checks it is
/// still current.
pub fn err_model() -> PathBuf {
    fixture("errs.json")
}

/// An integer field from a dada output JSON's `params` block.
pub fn param_i64(path: &Path, key: &str) -> i64 {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    v["params"][key]
        .as_i64()
        .unwrap_or_else(|| panic!("no integer params.{key} in {}", path.display()))
}

/// Sorted set of (sequence, abundance) from a `dada`/`dada-pseudo` output JSON.
pub fn asv_set(path: &Path) -> BTreeSet<(String, i64)> {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    v["asvs"]
        .as_array()
        .unwrap_or_else(|| panic!("no asvs array in {}", path.display()))
        .iter()
        .map(|a| {
            (
                a["sequence"].as_str().unwrap().to_ascii_uppercase(),
                a["abundance"].as_i64().unwrap(),
            )
        })
        .collect()
}

/// Set of (uppercased) sequences in a FASTA file.
pub fn fasta_seqs(path: &Path) -> BTreeSet<String> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('>') && !l.trim().is_empty())
        .map(|l| l.trim().to_ascii_uppercase())
        .collect()
}

/// A JSON output with the version tag stripped: it embeds the git hash, so it
/// differs between any two builds without meaning the results differ.
pub fn normalized(path: &Path) -> String {
    let raw = dada2_rs::misc::read_all_maybe_gz(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let text = String::from_utf8_lossy(&raw).into_owned();
    match text.find("\"dada2_rs_version\"") {
        Some(i) => {
            let end = text[i..].find(',').map(|j| i + j).unwrap_or(text.len());
            format!("{}{}", &text[..i], &text[end..])
        }
        None => text,
    }
}

/// Sorted names of the JSON outputs in `dir`.
pub fn json_outputs(dir: &Path) -> Vec<String> {
    let mut names: Vec<_> = std::fs::read_dir(dir)
        .expect("read output dir")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".json") || n.ends_with(".json.gz"))
        .collect();
    names.sort();
    names
}

/// Every JSON output in `a` must equal its namesake in `b`, version tag aside.
pub fn assert_same_outputs(a: &Path, b: &Path, label: &str) {
    let names = json_outputs(a);
    assert!(!names.is_empty(), "{label}: no outputs produced");
    for n in &names {
        assert_eq!(
            normalized(&a.join(n)),
            normalized(&b.join(n)),
            "{label}: {n} differs"
        );
    }
}
