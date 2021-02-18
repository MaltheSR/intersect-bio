use std::{convert::TryFrom, io};

use rust_htslib::bcf;

use crate::{ChromDict, ChromPos, Intersect};

impl<'a, R> Intersect<Records<'a, R>>
where
    R: bcf::Read,
{
    /// Create new intersect iterator from VCF readers.
    ///
    /// Chromosome dictionary is automatically created based on header information. VCF files
    /// are assumed to be sorted.
    pub fn vcfs(readers: &'a mut [R]) -> Self {
        let headers = readers.iter().map(|x| x.header()).collect::<Vec<_>>();

        let dict = ChromDict::from(headers.as_slice());

        let iters = readers
            .iter_mut()
            .map(|x| Records(x.records()))
            .collect::<Vec<_>>();

        Self::new(iters, dict)
    }
}

/// VCF record iterator.
///
/// This is a thin wrapper around the [`rust_htslib::bcf::Records`] iterator,
/// transforming the `rust_htslib` errors into `std::io::Error`.
///
/// Users should not need to interact with this struct, but it has to be public
/// since it is exposed as a type argument in the [`Intersect::vcfs`] constructor.
pub struct Records<'a, R>(bcf::Records<'a, R>)
where
    R: bcf::Read;

impl<'a, R> Iterator for Records<'a, R>
where
    R: bcf::Read,
{
    type Item = io::Result<bcf::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0
            .next()
            .map(|x| x.map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string())))
    }
}

impl ChromPos for bcf::Record {
    fn chrom(&self) -> &str {
        let rid = self.rid().expect("VCF record has no rid");

        let bytes = self
            .header()
            .rid2name(rid)
            .expect("cannot get VCF record contig name");

        std::str::from_utf8(bytes).expect("cannot convert VCF record contig name to UTF8")
    }

    fn pos(&self) -> u32 {
        u32::try_from(self.pos()).expect("cannot convert VCF position to u32")
    }
}

impl From<&[&bcf::header::HeaderView]> for ChromDict {
    fn from(headers: &[&bcf::header::HeaderView]) -> Self {
        ChromDict::from_intersection(headers.iter().map(|x| contigs(x)).collect())
    }
}

/// Get contig names from VCF header.
fn contigs(header: &bcf::header::HeaderView) -> Vec<String> {
    header
        .header_records()
        .into_iter()
        .filter_map(|x| match x {
            bcf::header::HeaderRecord::Contig { values, .. } => {
                let id = values
                    .into_iter()
                    .find(|(k, _)| k == "ID")
                    .expect("VCF header contig line did not contain 'ID' field");

                Some(id.1)
            }
            _ => None,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contigs_from_header() -> rust_htslib::errors::Result<()> {
        let ids = vec![1, 2, 4, 7];

        let mut header = bcf::Header::new();

        for id in ids.iter() {
            let header_contig_line = format!(r#"##contig=<ID={},length=10>"#, id);
            header.push_record(header_contig_line.as_bytes());
        }

        let vcf = bcf::Writer::from_path("/dev/null", &header, false, bcf::Format::BCF)?;
        let header = vcf.header();

        let expected = ids.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        assert_eq!(contigs(&header), expected);

        Ok(())
    }
}
