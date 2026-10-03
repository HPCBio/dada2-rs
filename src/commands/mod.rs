//! Subcommand handlers, one module per command group (issue #67). `main.rs`
//! parses the command line and dispatches here.

pub(crate) mod chimera;
mod common;
pub(crate) mod dada;
pub(crate) mod derep;
pub(crate) mod errors;
pub(crate) mod eval;
pub(crate) mod merge;
pub(crate) mod prep;
pub(crate) mod qc;
pub(crate) mod seqtable;
pub(crate) mod taxonomy;
