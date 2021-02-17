use std::{cmp, iter::FromIterator};

use indexmap::IndexSet;

use crate::ChromPos;

/// Ordered chromosome dictionary.
///
/// Efficient merging of positions across multiple ordered files requires pre-computing the subset
/// of chromosomes that occur in all files, as well as the ordering of this subset. We refer to this
/// information as a "chromosome dictionary".
///
/// Typically, the ordered chromosome IDs for each file can be obtained from a header (or similar),
/// and the chromosome dictionary may then be conveniently constructed using
/// [`from_intersection`](Self::from_intersection).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChromDict(IndexSet<String>);

impl ChromDict {
    /// Order positions relative to dictionary.
    ///
    /// If both positions are on chromosomes in the dictionary, returns the ordering of positions.
    /// Otherwise, returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::cmp::Ordering;
    /// # use merge_bio::ChromDict;
    /// let ids = vec!["1", "2"];
    /// let dict = ChromDict::from_ids(ids);
    ///
    /// assert_eq!(dict.compare(&("1", 2), &("2", 1)), Some(Ordering::Less));
    /// assert_eq!(dict.compare(&("2", 5), &("2", 2)), Some(Ordering::Greater));
    /// assert_eq!(dict.compare(&("1", 2), &("3", 2)), None);
    /// ```
    pub fn compare<T>(&self, first: &T, second: &T) -> Option<cmp::Ordering>
    where
        T: ChromPos,
    {
        if !(self.contains(first) && self.contains(second)) {
            return None;
        }

        if first.chrom() == second.chrom() {
            Some(first.pos().cmp(&second.pos()))
        } else {
            Some(Ord::cmp(
                &self.0.get_index_of(first.chrom()).unwrap(),
                &self.0.get_index_of(second.chrom()).unwrap(),
            ))
        }
    }

    /// Checks whether position is on a chromosome in the dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// # use merge_bio::ChromDict;
    ///
    /// let dict = ChromDict::from_ids(vec!["1", "2"]);
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

    /// Create dictionary from chromosome IDs.
    ///
    /// See [`from_intersection`](Self::from_intersection) for creating dictionary from multiple
    /// sources.
    ///
    /// # Examples
    ///
    /// ```
    /// # use merge_bio::ChromDict;
    /// let ids = vec!["1", "2"];
    ///
    /// let dict = ChromDict::from_ids(ids);
    /// ```
    pub fn from_ids<I, T>(ids: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        let set: IndexSet<String> = ids.into_iter().map(|x| x.to_string()).collect();

        Self::new(set)
    }

    /// Intersect dictionaries.
    ///
    /// Subset `self` to only contain entries also found in `other`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use merge_bio::ChromDict;
    /// let mut first_dict = ChromDict::from_ids(vec!["1", "2", "4", "5"]);
    /// let second_dict = ChromDict::from_ids(vec!["2", "3", "4"]);
    ///
    /// first_dict.intersect(&second_dict);
    /// assert_eq!(first_dict, ChromDict::from_ids(vec!["2", "4"]));
    /// ```
    pub fn intersect(&mut self, other: &Self) {
        self.0.retain(|x| other.0.contains(x))
    }

    /// Create dictionary from intersection of chromosome IDs from multiple sources.
    ///
    /// This takes IDs from multiple sources and finds the intersection.
    /// It is assumed that the IDs are sorted the same way in each source.
    ///
    /// # Examples
    ///
    /// ```
    /// # use merge_bio::ChromDict;
    /// let first_ids = vec!["1", "2", "3", "4"];
    /// let second_ids = vec!["1", "2", "4"];
    /// let third_ids = vec!["2", "3", "4"];
    ///
    /// let dict = ChromDict::from_intersection(vec![first_ids, second_ids, third_ids]);
    ///
    /// assert_eq!(dict, ChromDict::from_ids(vec!["2", "4"]));
    /// ```
    pub fn from_intersection<I, T>(mut id_sources: Vec<I>) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        if id_sources.is_empty() {
            return Self::default();
        }

        let mut dict = Self::from_iter(id_sources.pop().unwrap());

        id_sources.into_iter().for_each(|src| {
            dict.intersect(&Self::from_iter(src))
        });

        dict
    }

    /// Create new dictionary.
    fn new(ordering: IndexSet<String>) -> Self {
        Self(ordering)
    }
}

impl Default for ChromDict {
    fn default() -> Self {
        ChromDict::new(IndexSet::<String>::default())
    }
}

impl<T> FromIterator<T> for ChromDict
where
    T: ToString
{
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = T>
    {
        Self::new(iter.into_iter().map(|x| x.to_string()).collect())
    }
}
