use std::{convert::TryFrom, io};

use rust_htslib::bcf;

use crate::{ChromDict, ChromPos, Merge};

impl<'a, R> Merge<Records<'a, R>, bcf::Record>
where
    R: bcf::Read
{
    pub fn vcfs(readers: &'a mut [R]) -> Self {
        let headers = readers.iter().map(|x| x.header()).collect::<Vec<_>>();

        let dict = ChromDict::from(headers.as_slice());

        let iters = readers.iter_mut().map(|x| Records(x.records())).collect::<Vec<_>>();

        Self::new(iters, dict)
    }
}

struct Records<'a, R>(bcf::Records<'a, R>) where R: bcf::Read;

impl<'a, R> Iterator for Records<'a, R>
where
    R: bcf::Read
{
    type Item = io::Result<bcf::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|x| x.map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string())))
    }
}

impl ChromPos for bcf::Record {
    fn chrom(&self) -> &str {
        let rid = self.rid().expect("VCF record has no rid");

        let bytes = self.header().rid2name(rid).expect("cannot get VCF record contig name");

        std::str::from_utf8(bytes).expect("cannot convert VCF record contig name to UTF8")
    }

    fn pos(&self) -> u32 {
        u32::try_from(self.pos()).expect("cannot convert VCF position to u32")
    }
}

impl From<&[&bcf::header::HeaderView]> for ChromDict {
    fn from(headers: &[&bcf::header::HeaderView]) -> Self {
        ChromDict::from_merged_chromosomes(headers.into_iter().map(|x| contigs(x)).collect())
    }
}

fn contigs(header: &bcf::header::HeaderView) -> Vec<String> {
    header
        .header_records()
        .into_iter()
        .filter_map(|x| {
            match x {
                bcf::header::HeaderRecord::Contig { key, .. } => Some(key),
                _ => None,
            }
        })
        .collect()
}
