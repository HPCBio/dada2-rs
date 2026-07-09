//! dada2-rs: a Rust port of the DADA2 amplicon denoising algorithm.
//!
//! This library crate hosts the full module tree; the `dada2-rs` binary
//! (`src/main.rs`) is a thin CLI front-end over these modules.
#![allow(clippy::doc_overindented_list_items)]

pub mod chimera;
pub mod chimera_diagnostics;
pub mod cli;
pub mod cluster;
pub mod cluster_trace;
pub mod containers;
pub mod dada;
pub mod derep;
pub mod error;
pub mod error_models;
pub mod evaluate;
pub mod failed_uniques;
pub mod filter;
pub mod filter_trim;
pub mod kdist_calibrate;
pub mod kmers;
pub mod learn_errors;
pub mod loess;
pub mod merge_pairs;
pub mod misc;
pub mod nwalign;
pub mod pval;
pub mod remove_bimera;
pub mod remove_primers;
pub mod sequence_table;
pub mod summary;
pub mod taxonomy;
pub mod wfa;
