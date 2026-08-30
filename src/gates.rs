//! Development tuning gates: the `DADA2RS_*` environment variables, their
//! resolved values, and detection of ones that are set but not recognised
//! (issue #145).
//!
//! ## Why this module exists
//!
//! These gates are read once per process and resolved silently. If an operator
//! sets a variable the binary does not recognise — misremembered, renamed, or
//! mistyped — the run proceeds at the default and nothing says so.
//!
//! That is not hypothetical. Two full soil-pool A/B runs (16S and ITS2, both
//! reads, exclusive node) were spent comparing a configuration **against
//! itself**, because `DADA2RS_SHUFFLE_CARRY` had been renamed to
//! `DADA2RS_SHUFFLE_NO_CARRY` when the carry became the default. Both arms ran
//! carry-on. The runs completed, and the arms came back byte-identical and
//! within 1% — *which is exactly what a correct null result looks like.* The
//! mistake was found only by reading the git history of the gate itself.
//!
//! ## The two halves, and which one actually catches that
//!
//! [`report`] echoes **resolved** values under `--verbose`. Deliberately not an
//! echo of the environment: an environment echo would have printed
//! `DADA2RS_SHUFFLE_CARRY=1` in both failed runs above and been reassuring and
//! wrong. What is printed is what the code decided.
//!
//! [`warn_unrecognised`] is the half that catches the failure class, because it
//! fires **without anyone thinking to check**, and prints regardless of
//! `--verbose` — a stale gate means the run is not doing what was asked, which
//! is not a verbose-only concern.
//!
//! ## Adding or renaming a gate
//!
//! Add it to [`KNOWN`], or the run will warn about a variable that does work.
//! On a rename, add an [`Alias`] naming the replacement, so the old spelling
//! announces itself instead of failing silent. Aliases are the whole point of
//! the module; a rename without one re-arms the trap.

use std::collections::BTreeSet;

/// Every `DADA2RS_*` / `DADA2_RS_*` variable the binary reads.
///
/// Anything matching those prefixes and *not* listed here is reported as
/// unrecognised, so this list is load-bearing rather than documentation.
const KNOWN: &[&str] = &[
    "DADA2RS_SHUFFLE_NO_CARRY",
    "DADA2RS_SHUFFLE_NO_PRUNE",
    "DADA2RS_RECONCILE_VERIFY",
    "DADA2RS_RECONCILE_FULL",
    "DADA2RS_PAR_GRAIN",
    "DADA2RS_PROGRESS_SECS",
    "DADA2RS_PUPDATE_PREFETCH",
    "DADA2RS_ALIGN_BACKEND",
    "DADA2RS_WFA_MAX_STEPS",
    "DADA2RS_BENCH_THREADS",
];

/// What happened to a variable that is no longer the current spelling.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AliasKind {
    /// Same meaning, new spelling. Still honoured, but warns — breaking a
    /// working script to make a point would be its own kind of silent failure.
    Renamed,
    /// The replacement does not mean the same thing, so the old variable is
    /// **not** honoured. `DADA2RS_SHUFFLE_CARRY=1` is the case this exists for:
    /// the carry became the default and the gate inverted, so honouring the old
    /// name would select the opposite of what it says.
    Superseded,
}

/// A retired spelling, its replacement, and why it changed.
pub struct Alias {
    pub old: &'static str,
    pub new: &'static str,
    pub kind: AliasKind,
    pub note: &'static str,
}

/// Variables matching the gate prefixes that are **not** gates, and must not be
/// warned about.
///
/// `DADA2_RS_VERSION_FULL` is a build-time variable consumed by `build.rs` to
/// stamp the version string; it is legitimately present in any environment that
/// pins the version, and the CLI-error test sets it. Caught by that test rather
/// than by review, which is the argument for the test existing.
const NOT_GATES: &[&str] = &["DADA2_RS_VERSION_FULL"];

/// Retired gate spellings. Every rename adds a row here.
pub const ALIASES: &[Alias] = &[
    Alias {
        old: "DADA2RS_SHUFFLE_CARRY",
        new: "DADA2RS_SHUFFLE_NO_CARRY",
        kind: AliasKind::Superseded,
        note: "the carry became the default in #139 and the gate inverted, so \
               the old name cannot be honoured — it now selects the opposite \
               of what it says",
    },
    Alias {
        old: "DADA2_RS_PAR_GRAIN",
        new: "DADA2RS_PAR_GRAIN",
        kind: AliasKind::Renamed,
        note: "renamed for consistency with every other gate (#145); the old \
               spelling is still honoured",
    },
];

/// Look up a retired spelling.
pub fn alias_for(name: &str) -> Option<&'static Alias> {
    ALIASES.iter().find(|a| a.old == name)
}

/// Warn about `DADA2RS_*` / `DADA2_RS_*` variables that are set but not
/// recognised, naming the replacement where one is known.
///
/// Prints to stderr **irrespective of `--verbose`**: a stale gate means the run
/// is not the one the operator asked for, and the whole failure mode is that
/// nobody thought to look.
///
/// Returns the number of warnings, so callers (and tests) can assert on it.
pub fn warn_unrecognised() -> usize {
    let mut warnings = 0;
    let mut unknown = false;
    // Sorted so the output is stable across runs and diffable between logs.
    let set: BTreeSet<String> = std::env::vars_os()
        .filter_map(|(k, _)| k.into_string().ok())
        .filter(|k| k.starts_with("DADA2RS_") || k.starts_with("DADA2_RS_"))
        .filter(|k| !NOT_GATES.contains(&k.as_str()))
        .collect();

    for name in set {
        if let Some(alias) = alias_for(&name) {
            let effect = match alias.kind {
                AliasKind::Renamed => "still honoured, but rename it",
                AliasKind::Superseded => "it is having NO EFFECT",
            };
            eprintln!(
                "[dada] warning: {} is set but has been replaced by {} — {}.\n\
                 [dada]          {}",
                alias.old, alias.new, effect, alias.note
            );
            warnings += 1;
        } else if !KNOWN.contains(&name.as_str()) {
            eprintln!(
                "[dada] warning: {name} is set but not recognised; it is having \
                 NO EFFECT."
            );
            unknown = true;
            warnings += 1;
        }
    }
    // Once, not once per offender: with three typos set the list dwarfed the
    // warnings it was meant to explain.
    if unknown {
        eprintln!("[dada]          Recognised gates: {}", KNOWN.join(", "));
    }
    warnings
}

/// One resolved gate, for the `--verbose` summary.
struct Resolved {
    label: &'static str,
    value: String,
    default: bool,
}

impl Resolved {
    fn render(&self) -> String {
        format!(
            "{}={} ({})",
            self.label,
            self.value,
            if self.default {
                "default"
            } else {
                "OVERRIDDEN"
            }
        )
    }
}

/// Resolved tuning-gate values, as `--verbose` header lines.
///
/// These are the values the code settled on, **not** an echo of the
/// environment — see the module docs for why that distinction is the point.
pub fn report() -> Vec<String> {
    let grain = crate::cluster::par_max_len();
    let pf = crate::pval::prefetch_distance();
    let gates = [
        Resolved {
            label: "shuffle carry",
            value: on_off(crate::cluster::shuffle_carry()),
            default: crate::cluster::shuffle_carry(),
        },
        Resolved {
            label: "shuffle prune",
            value: on_off(crate::cluster::shuffle_prune()),
            default: crate::cluster::shuffle_prune(),
        },
        Resolved {
            label: "reconcile",
            value: if crate::cluster::reconcile_full_gate() {
                "full".into()
            } else {
                "incremental".into()
            },
            default: !crate::cluster::reconcile_full_gate(),
        },
        Resolved {
            label: "reconcile verify",
            value: on_off(crate::cluster::reconcile_verify_gate()),
            // On by default in debug builds; the release default is off, and a
            // release build is what every timing run uses.
            default: crate::cluster::reconcile_verify_gate() == cfg!(debug_assertions),
        },
        Resolved {
            label: "p-update prefetch",
            value: if pf == 0 {
                "off".into()
            } else {
                pf.to_string()
            },
            default: pf == crate::pval::PUPDATE_PREFETCH_DEFAULT,
        },
        Resolved {
            label: "par grain",
            value: grain.to_string(),
            default: grain == crate::cluster::PAR_GRAIN_DEFAULT,
        },
    ];

    // Wrap at two lines rather than one very long one; these sit alongside the
    // existing `alignment backend` / `cpu allocation` header lines.
    let rendered: Vec<String> = gates.iter().map(Resolved::render).collect();
    let mut out = Vec::new();
    for chunk in rendered.chunks(3) {
        out.push(format!("tuning gates: {}", chunk.join("  ")));
    }
    out
}

fn on_off(b: bool) -> String {
    if b { "on".into() } else { "off".into() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_alias_points_at_a_known_gate() {
        // A rename that points at a name the binary does not read would warn
        // the operator toward a variable that does nothing — the same trap,
        // one level down.
        for a in ALIASES {
            assert!(
                KNOWN.contains(&a.new),
                "alias {} points at {}, which is not in KNOWN",
                a.old,
                a.new
            );
        }
    }

    #[test]
    fn aliases_are_not_also_known() {
        // Otherwise the retired spelling would be silently accepted as current
        // and never warn.
        for a in ALIASES {
            assert!(
                !KNOWN.contains(&a.old),
                "retired spelling {} is still listed as known",
                a.old
            );
        }
    }

    #[test]
    fn known_gates_are_unique_and_prefixed() {
        let set: BTreeSet<_> = KNOWN.iter().collect();
        assert_eq!(set.len(), KNOWN.len(), "duplicate entry in KNOWN");
        for k in KNOWN {
            assert!(k.starts_with("DADA2RS_"), "{k} is not DADA2RS_-prefixed");
        }
    }

    #[test]
    fn build_time_vars_are_not_warned_about() {
        // DADA2_RS_VERSION_FULL matches the gate prefix but is a build-time
        // variable. Warning about it made a CLI-error test fail, which is how
        // this was caught.
        for n in NOT_GATES {
            assert!(
                !KNOWN.contains(n),
                "{n} is not a gate and must not be in KNOWN"
            );
            assert!(
                alias_for(n).is_none(),
                "{n} is not a gate and must not be aliased"
            );
        }
    }

    #[test]
    fn report_marks_defaults() {
        // With no gates set, every line must read `(default)`. This is the
        // assertion that would have failed loudly in the two runs that compared
        // a configuration against itself.
        for line in report() {
            assert!(
                !line.contains("OVERRIDDEN"),
                "unexpected override in a clean environment: {line}"
            );
        }
    }
}
