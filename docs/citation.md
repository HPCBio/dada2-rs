# Citation

`dada2-rs` is a **reimplementation** of DADA2. If you use it, please cite the
original work that describes the algorithm. Machine-readable citation metadata
for this project lives in
[`CITATION.cff`](https://github.com/HPCBio/dada2-rs/blob/main/CITATION.cff) at
the repository root.

## Cite the original DADA2 algorithm

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. **DADA2:
  High-resolution sample inference from Illumina amplicon data.** *Nature
  Methods*. 2016;13:581–583.
  doi:[10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- Rosen MJ, Callahan BJ, Fisher DS, Holmes SP. **Denoising PCR-amplified
  metagenome data.** *BMC Bioinformatics*. 2012;13:283.
  doi:[10.1186/1471-2105-13-283](https://doi.org/10.1186/1471-2105-13-283)

## Provenance & guidelines

This project follows the [rewrites.bio](https://rewrites.bio) guidelines as
closely as possible, with the aim of ensuring the original DADA2 developers and
contributors retain clear credit for their work, following the original
implementation's details, including tests and benchmarks, and releasing as open
source under the original licensing.

!!! note "Results may differ slightly from R DADA2"
    Some implementation details (the R/C++ → Rust conversion, including the error
    models) cause small differences from R DADA2. We strive to reproduce results
    as closely as possible, and provide the ability to use R and/or Python for
    custom error-model analysis to more closely emulate the original
    implementation.
