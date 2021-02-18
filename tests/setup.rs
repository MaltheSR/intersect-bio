use std::{fs, io, path, process};

use rand::{prelude::IteratorRandom, Rng, SeedableRng};

use rust_htslib::bcf;

const N_CONTIGS: usize = 15;
const MAX_CONTIG: u8 = 20;
const N_POSITIONS: usize = 500;
const MAX_POSITION: i64 = 1000;

/// Write a random VCF
///
/// The VCF will contain positions in N_CONTIGS random contigs with IDs between 1 and MAX_CONTIG.
/// Each conting will contain N_POSITIONS random positions between 1 and MAX_POSITION. All alleles
/// are set to A/C, an AC=1 flag is set in the INFO column, and a random genotype is chosen.
pub fn write_vcf<P>(path: P, seed: u64) -> rust_htslib::errors::Result<()>
where
    P: AsRef<path::Path>,
{
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    // Setup header with random contigs
    let mut header = bcf::Header::new();

    let possible_contigs: Vec<u8> = (1..MAX_CONTIG).collect();
    let mut contigs = possible_contigs.iter().choose_multiple(&mut rng, N_CONTIGS);
    contigs.sort();

    for contig in contigs.iter() {
        let header_contig_line = format!(r#"##contig=<ID={},length=10>"#, contig);
        header.push_record(header_contig_line.as_bytes());
    }

    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    header.push_record(header_gt_line.as_bytes());

    let header_ac_line = r#"##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">"#;
    header.push_record(header_ac_line.as_bytes());

    header.push_sample(format!("sample{}", seed).as_bytes());

    // Setup VCF
    let mut vcf = bcf::Writer::from_path(path, &header, false, bcf::Format::VCF)?;

    let possible_positions: Vec<i64> = (1..MAX_POSITION as i64).collect();

    // Write records with random positions and genotypes for each contig
    for contig in contigs.into_iter() {
        let mut positions = possible_positions
            .iter()
            .choose_multiple(&mut rng, N_POSITIONS);
        positions.sort();

        let genotypes = (&mut rng)
            .sample_iter(rand::distributions::Uniform::new(0i32, 2i32))
            .take(N_POSITIONS)
            .map(|x| {
                [
                    bcf::record::GenotypeAllele::Unphased(x),
                    bcf::record::GenotypeAllele::Phased(x),
                ]
            })
            .collect::<Vec<_>>();

        for (position, genotype) in positions.into_iter().zip(genotypes.into_iter()) {
            let mut record = vcf.empty_record();

            let rid = vcf.header().name2rid(contig.to_string().as_bytes())?;
            record.set_rid(Some(rid));

            record.set_pos(*position);

            record.set_alleles(&[b"A", b"C"])?;

            record.push_info_integer(b"AC", &[1])?;

            record.push_genotypes(&genotype)?;

            vcf.write(&record)?;
        }
    }

    Ok(())
}

/// Index VCF
///
/// Use bcftools to create a CSI index. This will fail if bcftools is not in the $PATH.
pub fn index_vcf<P>(path: P) -> io::Result<()>
where
    P: AsRef<path::Path>,
{
    let vcf_path = path::PathBuf::from(path.as_ref());

    let status = process::Command::new("bcftools")
        .arg("index")
        .arg(vcf_path)
        .status()?;

    if status.success() {
        Ok(())
    } else {
        Err(io::Error::new(io::ErrorKind::Other, "failed to index VCF"))
    }
}

/// Intersect VCFs
///
/// Use bcftools to intersect VCFs by chaining bcftools merge and filtering any sites with missing
/// genotypes. This will fail if bcftools is not in the $PATH.
pub fn bcftools_intersect<P>(in_paths: &[P], out_path: P) -> io::Result<()>
where
    P: AsRef<path::Path>,
{
    let in_paths = in_paths
        .iter()
        .map(|x| x.as_ref().as_os_str())
        .collect::<Vec<_>>();

    let tmp_path = out_path.as_ref().with_extension("gz.tmp");

    // Create merged VCF
    let tmp_status = process::Command::new("bcftools")
        .arg("merge")
        .args(in_paths)
        .arg("-o")
        .arg(tmp_path.clone())
        .status()?;

    // Filter missing sites, leaving only intersection
    let status = process::Command::new("bcftools")
        .args(&["view", "-e", r#"GT[*] = "mis""#, "-O", "z", "-o"])
        .arg(out_path.as_ref().as_os_str())
        .arg(tmp_path.clone())
        .status()?;

    index_vcf(out_path)?;

    // Clean up
    fs::remove_file(tmp_path)?;

    if tmp_status.success() && status.success() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::Other,
            "failed to intersect VCFs with bcftools",
        ))
    }
}
