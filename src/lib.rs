mod merge;
mod chrom_dict;

#[cfg(feature = "rust-htslib")]
mod rust_htslib;

pub use self::{chrom_dict::ChromDict, merge::Merge};

pub trait ChromPos {
    fn chrom(&self) -> &str;

    fn pos(&self) -> u32;

    fn colocated(&self, other: &Self) -> bool {
        self.chrom() == other.chrom() && self.pos() == other.pos()
    }
}

impl<T> ChromPos for (T, u32)
where
    T: AsRef<str>,
{
    fn chrom(&self) -> &str {
        self.0.as_ref()
    }

    fn pos(&self) -> u32 {
        self.1
    }
}
