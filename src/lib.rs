mod merge;
mod chrom_dict;

#[cfg(feature = "rust-htslib")]
mod rust_htslib;

pub use self::{chrom_dict::ChromDict, merge::Merge};

/// A genomic position.
///
/// Trait for an entity whose location along a genome can be described by an integer coordinate
/// along some chromosome (or similar, e.g. contig).
pub trait ChromPos {
    /// Get the chromosome ID.
    fn chrom(&self) -> &str;

    /// Get the position along the chromosome.
    fn pos(&self) -> u32;

    /// Check whether two position are on the same chromosome with the same position along that
    /// chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// # use merge_bio::ChromPos;
    /// assert!(("1", 1).intersect(&("1", 1)));
    /// ```
    fn intersect(&self, other: &Self) -> bool {
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
