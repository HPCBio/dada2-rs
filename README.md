# dada2-rs

This is an experimental implementation of DADA2 in Rust, using Claude Sonnet 4.6 for the bulk of the work. 

This work will be following [rewrites.bio](https://rewrites.bio), including full credit to the original work by Ben Callahan and other DADA2 contributors. 

## Plans

The short term plan: 

1) Port key code over to Rust: dereplication, error modeling, denoising, merging, chimera removal
2) Decouple any underlying C++ code from R and reimplement key steps as subcommands, R classes as Rust structs/classes, etc.
3) Allow intermediate outputs (in JSON) that can be evaluated for debugging purposes or for plotting in R, Python, etc.
4) Add basic regression tests that follow those within the original dada2 repository and expect results similar to those expected

We do anticipate porting read trimming or taxonomic classification at a later point

## AI Assistance Disclosure
This tool was written with the assistance of AI coding agents, specificall Claude Code, using Sonnet 4.6. All commits using AI are noted.

Correctness is validated by comparing output against DADA2 v1.36 on a suite of real sequencing datasets - not by manual code review alone. 
AI generated the implementation; humans defined the validation criteria, made some key coding updates, and verified results.
