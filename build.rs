use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // Re-run the build script when the git HEAD or refs change so the embedded
    // SHA tracks the working copy. Missing files are ignored by Cargo.
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");
    println!("cargo:rerun-if-changed=.git/refs/tags");
    println!("cargo:rerun-if-changed=.git/packed-refs");

    let cargo_version = env!("CARGO_PKG_VERSION");

    // Short SHA of HEAD (8 chars). Empty when not in a git checkout or git is
    // unavailable (e.g. building from a release tarball).
    let sha = Command::new("git")
        .args(["rev-parse", "--short=8", "HEAD"])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
        .unwrap_or_default();

    // True when HEAD is exactly the tag matching the current cargo version
    // (either `vX.Y.Z` or `X.Y.Z`). Tagged release builds emit a clean
    // `X.Y.Z` version; everything else gets the `-<sha>` suffix.
    let is_release_tag = Command::new("git")
        .args(["describe", "--exact-match", "--tags", "HEAD"])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
        .map(|tag| tag == format!("v{cargo_version}") || tag == cargo_version)
        .unwrap_or(false);

    let version = if is_release_tag || sha.is_empty() {
        cargo_version.to_string()
    } else {
        format!("{cargo_version}-{sha}")
    };

    println!("cargo:rustc-env=DADA2_RS_VERSION_FULL={version}");
}
