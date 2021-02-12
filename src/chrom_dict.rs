use std::cmp;

use indexmap::IndexSet;

use crate::ChromPos;

/// A chromosome dictionary.
///
/// Contains an ordered set of chromosome identifers. Merging positions from multiple sources
/// requires creating a chromosome dictionary containing the identifiers of all chromosomes
/// that occur in all input sources, in the same order that they occur in those sources.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChromDict(IndexSet<String>);

impl ChromDict {
    /// Order positions according to chromosome dictionary.
    ///
    /// If the chromosome identifiers of both positions are contained in the dict,
    /// returns an ordering relative to the chromosome dictionary. Otherwise, returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::cmp::Ordering;
    /// use merge_bio::ChromDict;
    ///
    /// let chromosomes = vec!["1", "2"];
    /// let dict = ChromDict::from_chromosomes(chromosomes);
    ///
    /// let first = ("1", 34);
    /// let second = ("1", 50);
    /// assert_eq!(dict.compare(&first, &second), Some(Ordering::Less));
    ///
    /// let third = ("2", 5);
    /// assert_eq!(dict.compare(&third, &second), Some(Ordering::Greater));
    ///
    /// let fourth = ("3", 11);
    /// assert_eq!(dict.compare(&first, &fourth), None);
    /// ```
    pub fn compare<T>(&self, first: &T, second: &T) -> Option<cmp::Ordering>
    where
        T: ChromPos,
    {
        // Check that both chromosomes are contained
        if !(self.contains(first) && self.contains(second)) {
            return None;
        }

        // Get ordering
        if first.chrom() == second.chrom() {
            // Same chromosome, order by position within chromosome
            Some(first.pos().cmp(&second.pos()))
        } else {
            // Different chromosomes, order by chromosome order
            Some(Ord::cmp(
                &self.0.get_index_of(first.chrom()).unwrap(),
                &self.0.get_index_of(second.chrom()).unwrap(),
            ))
        }
    }

    /// Checks whether position is located on a chromosome contained in chromosome dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::cmp::Ordering;
    /// use merge_bio::ChromDict;
    ///
    /// let chromosomes = vec!["1", "2"];
    /// let dict = ChromDict::from_chromosomes(chromosomes);
    ///
    /// let first = ("1", 34);
    /// assert!(dict.contains(&first));
    ///
    /// let second = ("3", 11);
    /// assert!(!dict.contains(&second));
    /// ```
    pub fn contains<T>(&self, chrom_pos: &T) -> bool
    where
        T: ChromPos,
    {
        self.0.contains(chrom_pos.chrom())
    }

    /// Create chromosome dictionary from a single set of chromosome identifiers.
    ///
    /// This assumes that chromosomes have already been intersected and ordered.
    /// For most use-cases, the [`from_merged_chromosomes`](Self::from_merged_chromosomes)
    /// constructor will likely be more convenient.
    ///
    /// # Examples
    ///
    /// ```
    /// use merge_bio::ChromDict;
    ///
    /// let chromosomes = vec!["1", "2"];
    ///
    /// let dict = ChromDict::from_chromosomes(chromosomes);
    /// ```
    pub fn from_chromosomes<I, T>(chromosomes: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        let set: IndexSet<String> = chromosomes.into_iter().map(|x| x.to_string()).collect();

        Self::new(set)
    }

    /// Create chromosome dictionary by merging chromosome identifiers from multiple sources.
    ///
    /// Note that this assumes that the subset of chromosome identifiers that occur in all source
    /// occur in the same ordering in each of the sources.
    ///
    /// # Examples
    ///
    /// ```
    /// use merge_bio::ChromDict;
    ///
    /// let first = vec!["1", "2", "3", "4"];
    /// let second = vec!["1", "2", "4"];
    /// let third = vec!["2", "3", "4"];
    ///
    /// let dict = ChromDict::from_merged_chromosomes(vec![first, second, third]);
    ///
    /// assert_eq!(dict, ChromDict::from_chromosomes(vec!["2", "4"]));
    /// ```
    pub fn from_merged_chromosomes<I, T>(chromosomes: Vec<I>) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        let mut set = IndexSet::<String>::default();

        for (i, chrom) in chromosomes.into_iter().enumerate() {
            let new_set: IndexSet<String> = chrom.into_iter().map(|x| x.to_string()).collect();

            if i == 0 {
                set = new_set;
            } else {
                set.retain(|x| new_set.contains(x));
            }
        }

        Self::new(set)
    }

    /// Create new chromosome dictionary.
    fn new(ordering: IndexSet<String>) -> Self {
        Self(ordering)
    }
}
