//! Runtime CPU-feature dispatch for the DP kernels (issue #336).
//!
//! The kernels are plain Rust that LLVM auto-vectorizes, so their width is set
//! at compile time. A baseline x86-64 build (`cargo install`, Bioconda, the
//! portable Docker image) gets SSE2 only, and its NW kernel ran ~1.5x slower
//! than a `target-cpu=native` build. Each hot kernel therefore has its body in
//! an `#[inline(always)]` function, compiled twice: once into a baseline entry
//! and once into an `#[target_feature(enable = "avx2")]` entry. [`avx2`]
//! chooses between them, once per process.
//!
//! The kernels are integer-only, so both entries give byte-identical output;
//! the `*_variants_match` tests check that on every host that has AVX2.
//!
//! `DADA2RS_SIMD=baseline` forces the baseline entries, so one binary can A/B
//! the two paths. On aarch64, NEON is baseline and there is nothing to choose.

use std::sync::OnceLock;

/// The instruction set the dispatched kernels use in this process.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    Baseline,
    Avx2,
}

impl Level {
    pub fn label(self) -> &'static str {
        match self {
            Level::Baseline => "baseline",
            Level::Avx2 => "avx2",
        }
    }
}

/// What the host supports, ignoring `DADA2RS_SIMD`.
pub fn detected() -> Level {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            return Level::Avx2;
        }
    }
    Level::Baseline
}

/// The level in effect: [`detected`], unless `DADA2RS_SIMD=baseline` forces
/// the baseline entries. Any other value is ignored with a warning.
pub fn level() -> Level {
    static VALUE: OnceLock<Level> = OnceLock::new();
    *VALUE.get_or_init(|| match std::env::var("DADA2RS_SIMD").as_deref() {
        Err(_) | Ok("") | Ok("auto") => detected(),
        Ok("baseline") => Level::Baseline,
        Ok(other) => {
            eprintln!(
                "[dada2-rs] warning: DADA2RS_SIMD={other} is not recognised \
                 (expected `baseline` or `auto`); using auto-detection"
            );
            detected()
        }
    })
}

/// True when the AVX2 entries are to be used. Callers on other architectures
/// never reach the AVX2 entries, which do not exist there.
#[inline]
pub fn avx2() -> bool {
    level() == Level::Avx2
}

/// `popcount(a & b)` over two equal-length bitmaps: the shared-k-mer count of
/// the #295 screen bound. A baseline x86-64 build has no POPCNT instruction,
/// so an AVX2+POPCNT entry is used when the CPU has both (every AVX2 CPU in
/// practice) and [`level`] allows AVX2.
#[inline]
pub fn and_popcount(a: &[u64], b: &[u64]) -> u32 {
    #[cfg(target_arch = "x86_64")]
    if popcnt() {
        // SAFETY: `popcnt` is true only when the CPU has AVX2 and POPCNT.
        return unsafe { and_popcount_avx2(a, b) };
    }
    and_popcount_body(a, b)
}

#[inline(always)]
fn and_popcount_body(a: &[u64], b: &[u64]) -> u32 {
    a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum()
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
fn and_popcount_avx2(a: &[u64], b: &[u64]) -> u32 {
    and_popcount_body(a, b)
}

#[cfg(target_arch = "x86_64")]
fn popcnt() -> bool {
    static VALUE: OnceLock<bool> = OnceLock::new();
    *VALUE.get_or_init(|| avx2() && std::is_x86_feature_detected!("popcnt"))
}

/// For the `*_variants_match` tests: true when the host can run the AVX2
/// entries. Says so when it cannot, since the test then checks nothing.
#[cfg(all(test, target_arch = "x86_64"))]
pub(crate) fn host_has_avx2(test: &str) -> bool {
    let yes = detected() == Level::Avx2;
    if !yes {
        eprintln!("{test}: skipped, this host has no AVX2");
    }
    yes
}
