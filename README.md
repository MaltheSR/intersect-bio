# intersect-bio

[![GitHub Actions status](https://github.com/malthesr/intersect-bio/workflows/CI/badge.svg)](https://github.com/malthesr/intersect-bio/actions)

**intersect-bio** is a small Rust crate for efficient iteration over intersecting genomic positions in biological files: support is included for VCF files if the feature flag `rust-htslib` is set (default), and the workload for implementing new file types is designed to be small. See the crate documentation for more extensive details (see "Documentation" section below).

## Usage

**intersect-bio** is not published on [crates.io](https://crates.io/), as it is mainly thought of as a personal helper crate for now.

Nevertheless, to use **angsd-io**, you can depend on this github repo. To do so, add the following to the `[dependencies]` section of your `Cargo.toml`:

```
intersect-bio = { git = "https://github.com/malthesr/intersect-bio.git" }
```

For more information, including on how to depend on a particular commit, see [here](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#specifying-dependencies-from-git-repositories).

## Examples

The [`examples`](examples/) sub-directory contain a runnable example of illustrative basic usage of VCF intersection. The following will intersect an arbitrary number of VCFs, printing the `CHROM` and `POS` fields for sites founds in all files, as well as the genotype of the first sample contained in each:

```
cargo run --release --example intersect_vcfs [PATH_TO_VCFS...]
```

## Documentation

The documentation can be built and viewed locally by running

```
cargo doc --open
```
