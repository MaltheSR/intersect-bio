use std::{fs, io, path};

mod setup;

use setup::{index_vcf, write_vcf, bcftools_intersect};

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
    P: AsRef<path::Path>
{
    let mut dir = vcf_dir();
    dir.push(name);
    dir
}


#[test]
fn intersect_vcfs() -> io::Result<()> {
    fs::create_dir_all(vcf_dir())?;

    // Required file paths
    let vcf_paths = VCF_NAMES.iter().map(|name| vcf_path(name)).collect::<Vec<_>>();

    let intersect_vcf_path = vcf_path(INTERSECT_VCF_NAME);

    // Create files if they do not exist
    if !(vcf_paths.iter().all(|x| x.exists()) && intersect_vcf_path.exists()) {
        for (i, path) in vcf_paths.iter().enumerate() {
            write_vcf(path, i as u64)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
            index_vcf(path)?;
        }

        // Intersect them using bcftools
        bcftools_intersect(&vcf_paths, intersect_vcf_path)?;
    }

    //
    assert!(false);

    Ok(())
}
