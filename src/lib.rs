#![deny(missing_docs)]

//! This crate provides an iterator over the intersection of an arbitrary number of biological files
//! containing genomic sites pre-sorted by chromosome and position.
//!
//! After calculating the sorted intersection of chromosome IDs in all files (which can typically be
//! done from header information in trivial time), this is achieved with a single pass through each
//! file.  That is, the time complexity is linear in the total number of sites, and the only memory
//! required is to hold a single site per input file in RAM at any given time.
//!
//! If the `rust-htslib` feature flag is set, such intersect iteration comes pre-supported for VCF
//! files for convenience and illustration. However, the goal is also to allow easy implementation
//! for other file types through use of generics. Each of these points is described below.
//!
//! # Implementing a new file format
//!
//! The crate is designed around the situation where the chromosome IDs can be quickly extracted
//! from header information or some other source before iteration. This information should be stored
//! in a [`ChromDict`] (a "chromosome dictionary"), which has helper methods for construction once
//! IDs from each source has been obtained.
//!
//! Apart from this, a (fallible) iterator over each input source must be implemented.  Each
//! iteration must yield a `std::io::Result<T>`, where `T` is [`ChromPos`].
//!
//! Once these requirements are met, intersection is provided by passing any number of iterators
//! and their corresponding chromosome dictionary to the [`Intersect`] iterator.
//!
//! # Intersecting VCFs
//!
//! If the `rust-htslib` feature flag is set, `intersect-bio` comes pre-packaged with support for
//! intersecting VCF files. Note that VCFs are assumed to be pre-sorted; furthermore, the contig
//! records in each VCF header is assumed to be sorted according the same sort order as used for the
//! files themselves.
//!
//! ``` no_run
//! use intersect_bio::{ChromPos, Intersect};
//! use rust_htslib::bcf;
//!
//! let paths = vec!["test1.vcf.gz", "test2.vcf.gz", "test3.vcf.gz"];
//!
//! let mut readers = paths
//!     .iter()
//!     .map(|p| bcf::Reader::from_path(p))
//!     .collect::<rust_htslib::errors::Result<Vec<_>>>()
//!     .expect("cannot open VCF reader");
//!
//! let intersection = Intersect::vcfs(readers.as_mut_slice());
//!
//! for site in intersection {
//!     let site = site.expect("failed to read site");
//!
//!     let first_record = &site[0];
//!
//!     assert!(
//!         site.iter().skip(1).all(|record| {
//!             record.chrom() == first_record.chrom() && record.pos() == first_record.pos()
//!         })
//!     );
//! }
//! ```
//!
//! A similar, runnable example is contained in the `examples/` directory of the repository.

mod chrom_dict;
mod intersect;

#[cfg(feature = "rust-htslib")]
mod rust_htslib;

pub use self::{chrom_dict::ChromDict, intersect::Intersect};

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
    /// # use intersect_bio::ChromPos;
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
