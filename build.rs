use std::path::{Path, PathBuf};
use std::process::Command;

/// Run `git` with `args` in the crate directory; trimmed stdout on success.
fn git(args: &[&str]) -> Option<String> {
    Command::new("git")
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
}

/// Re-run the build script when HEAD or the refs change, so the embedded SHA
/// tracks the working copy.
///
/// The paths come from git rather than being hard-coded under `.git/`: in a
/// worktree `.git` is a file, `HEAD` lives in the per-worktree git dir and the
/// refs in the shared common dir. Only existing paths are watched, because
/// Cargo treats a missing `rerun-if-changed` path as always changed and would
/// then recompile the crate on every build.
fn watch_git_state() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    // Relative output is relative to the cwd, which for a build script is the
    // manifest dir.
    let resolve = |p: String| -> PathBuf { manifest.join(p) };
    let Some(git_dir) = git(&["rev-parse", "--git-dir"]).map(resolve) else {
        return;
    };
    let common = git(&["rev-parse", "--git-common-dir"])
        .map(resolve)
        .unwrap_or_else(|| git_dir.clone());

    let watched = [
        git_dir.join("HEAD"),
        common.join("refs/heads"),
        common.join("refs/tags"),
        common.join("packed-refs"),
    ];
    for p in watched.iter().filter(|p| Path::exists(p)) {
        println!("cargo:rerun-if-changed={}", p.display());
    }
}

/// The `rmath-ppois` diagnostic arm (#277): link R's own libR so calc_pA can
/// call its ppois. Returns the version-string tag that marks such a build.
fn link_r_for_rmath_ppois() -> &'static str {
    if std::env::var_os("CARGO_FEATURE_RMATH_PPOIS").is_none() {
        return "";
    }
    println!("cargo:rerun-if-env-changed=R_HOME");
    let r_home = std::env::var("R_HOME")
        .ok()
        .filter(|s| !s.is_empty())
        .or_else(|| {
            let out = Command::new("R").arg("RHOME").output().ok()?;
            out.status
                .success()
                .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
        });
    let r_home = r_home.expect("feature rmath-ppois: set R_HOME or put `R` on PATH");
    let lib = Path::new(&r_home).join("lib");
    println!("cargo:rustc-link-search=native={}", lib.display());
    println!("cargo:rustc-link-lib=dylib=R");
    // Embed the path so the binary finds libR without LD_LIBRARY_PATH.
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib.display());
    "+rmath-ppois"
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    watch_git_state();
    let arm = link_r_for_rmath_ppois();
    // Allow callers to inject the version (e.g. Docker builds without .git).
    println!("cargo:rerun-if-env-changed=DADA2_RS_VERSION_FULL");

    let cargo_version = env!("CARGO_PKG_VERSION");

    // If the caller pre-set DADA2_RS_VERSION_FULL, trust it verbatim. This is
    // the escape hatch for Docker / CI builds where .git is unavailable.
    if let Ok(injected) = std::env::var("DADA2_RS_VERSION_FULL") {
        let injected = injected.trim();
        if !injected.is_empty() {
            println!("cargo:rustc-env=DADA2_RS_VERSION_FULL={injected}{arm}");
            return;
        }
    }

    // Short SHA of HEAD (8 chars). Empty when not in a git checkout or git is
    // unavailable (e.g. building from a release tarball).
    let sha = git(&["rev-parse", "--short=8", "HEAD"]).unwrap_or_default();

    // True when HEAD is exactly the tag matching the current cargo version
    // (either `vX.Y.Z` or `X.Y.Z`). Tagged release builds emit a clean
    // `X.Y.Z` version; everything else gets the `-<sha>` suffix.
    let is_release_tag = git(&["describe", "--exact-match", "--tags", "HEAD"])
        .map(|tag| tag == format!("v{cargo_version}") || tag == cargo_version)
        .unwrap_or(false);

    let version = if is_release_tag || sha.is_empty() {
        cargo_version.to_string()
    } else {
        format!("{cargo_version}-{sha}")
    };

    println!("cargo:rustc-env=DADA2_RS_VERSION_FULL={version}{arm}");
}
