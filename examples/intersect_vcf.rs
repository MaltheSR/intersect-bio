//! Intersect multiple VCF genotypes
//!
//! For each site in intersection, write the chromosome and position of the site,
//! and the genotype of the first sample in each VCF at that site.

#![cfg(feature = "rust-htslib")]

use std::io::{self, Write};

use rust_htslib::bcf;

use intersect_bio::{ChromPos, Intersect};

fn main() -> io::Result<()> {
    let paths = std::env::args().skip(1);

    // Create VCF readers from input
    let mut readers = paths
        .map(|p| bcf::Reader::from_path(p))
        .collect::<rust_htslib::errors::Result<Vec<_>>>()
        .expect("cannot open VCF reader");

    // Create intersect iterator over readers
    let intersect = Intersect::vcfs(readers.as_mut_slice());

    // Create writer to stdout
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = io::BufWriter::new(handle);

    // Print chrom, pos, and first genotype in each record
    for site in intersect {
        let site = site?;

        let chrom = site[0].chrom();
        let pos = site[0].pos();
        let gt = site.iter().map(|record| {
            format!("{}", record.genotypes().expect("cannot get record genotypes").get(0))
        }).collect::<Vec<_>>();

        writeln!(writer, "{}\t{}\t{}", chrom, pos, gt.join(";"))?;
    }

    Ok(())
}
