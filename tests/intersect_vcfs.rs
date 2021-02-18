use std::{fs, io, path};

use rust_htslib::bcf::{self, Read};

use intersect_bio::{ChromPos, Intersect};

mod setup;

use setup::{bcftools_intersect, index_vcf, write_vcf};

const VCF_DIR: &str = "tests/data/";
const VCF_NAMES: [&str; 3] = ["test1.vcf.gz", "test2.vcf.gz", "test3.vcf.gz"];
const INTERSECT_VCF_NAME: &str = "intersect.vcf.gz";

/// Creates the full path to the VCF directory.
fn vcf_dir() -> path::PathBuf {
    let mut dir = path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    dir.push(VCF_DIR);
    dir
}

/// Creates a full path to a VCF file from the file name.
fn vcf_path<P>(name: P) -> path::PathBuf
where
    P: AsRef<path::Path>,
{
    let mut dir = vcf_dir();
    dir.push(name);
    dir
}

/// Open a VCF reader.
fn vcf_reader<P>(path: P) -> io::Result<bcf::Reader>
where
    P: AsRef<path::Path>,
{
    bcf::Reader::from_path(path).map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))
}

#[test]
fn intersect_vcfs() -> io::Result<()> {
    fs::create_dir_all(vcf_dir())?;

    // Required file paths
    let vcf_paths = VCF_NAMES
        .iter()
        .map(|name| vcf_path(name))
        .collect::<Vec<_>>();

    let intersect_vcf_path = vcf_path(INTERSECT_VCF_NAME);

    // Create files if they do not exist
    if !(vcf_paths.iter().all(|x| x.exists()) && intersect_vcf_path.exists()) {
        for (i, path) in vcf_paths.iter().enumerate() {
            write_vcf(path, i as u64)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
            index_vcf(path)?;
        }

        // Intersect them using bcftools
        bcftools_intersect(&vcf_paths, intersect_vcf_path.clone())?;
    }

    // Setup iterators
    let mut bcftools_vcf = vcf_reader(intersect_vcf_path.clone())?;
    let bcftools_records = bcftools_vcf
        .records()
        .map(|x| x.map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string())));

    let mut vcfs = vcf_paths
        .iter()
        .map(|p| vcf_reader(p))
        .collect::<io::Result<Vec<_>>>()?;
    let intersect = Intersect::vcfs(&mut vcfs);

    // Check all records match
    let mut counter = 0;
    for (intersected_site, bcftools_site) in intersect.zip(bcftools_records) {
        let intersected_site = intersected_site?;

        // Sanity check that intersecting sites actually intersect
        assert!(intersected_site
            .iter()
            .all(|x| x.intersect(&intersected_site[0])));

        // Check that intersecting sites match bcftools
        assert!(intersected_site[0].intersect(&bcftools_site?));

        counter += 1;
    }

    // Check that iterators had equal length (since zip is not strict)
    let mut bcftools_vcf = vcf_reader(intersect_vcf_path)?;
    let bcftools_records = bcftools_vcf.records();
    let n = bcftools_records.count();

    assert_eq!(counter, n);

    Ok(())
}
