use std::cmp;

use indexmap::IndexSet;

use crate::ChromPos;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SequenceDict(IndexSet<String>);

impl SequenceDict {
    /// Order [`ChromPos`] instances according to sequence dictionary.
    ///
    /// If both instances are contained in the dict, returns an ordering relative to
    /// the sequence dictionary. Otherwise, returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::cmp::Ordering;
    /// use merge_bio::SequenceDict;
    ///
    /// let sequence = vec!["1", "2"];
    /// let dict = SequenceDict::from_sequence(sequence);
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

    /// Checks whether [`ChromPos`] is on chromosome contained in sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::cmp::Ordering;
    /// use merge_bio::SequenceDict;
    ///
    /// let sequence = vec!["1", "2"];
    /// let dict = SequenceDict::from_sequence(sequence);
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

    /// Create sequence dict from a single sequence of chromosomes.
    ///
    /// This assumes that sequence chromosomes have already been merged.
    /// For most use-cases, the [`from_merged_sequence`](#method.from_merged_sequence) constructor
    /// will likely be more convenient.
    ///
    /// # Examples
    ///
    /// ```
    /// use merge_bio::SequenceDict;
    ///
    /// let sequence = vec!["1", "2"];
    ///
    /// let dict = SequenceDict::from_sequence(sequence);
    /// ```
    pub fn from_sequence<I, T>(sequence: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        let set: IndexSet<String> = sequence.into_iter().map(|x| x.to_string()).collect();

        Self::new(set)
    }

    /// Create sequence dict by merging from multiple sequences of chromosomes.
    ///
    /// # Examples
    ///
    /// ```
    /// use merge_bio::SequenceDict;
    ///
    /// let first = vec!["1", "2", "3", "4"];
    /// let second = vec!["1", "2", "4"];
    /// let third = vec!["2", "3", "4"];
    ///
    /// let dict = SequenceDict::from_merged_sequences(vec![first, second, third]);
    ///
    /// assert_eq!(dict, SequenceDict::from_sequence(vec!["2", "4"]));
    /// ```
    pub fn from_merged_sequences<I, T>(sequences: Vec<I>) -> Self
    where
        I: IntoIterator<Item = T>,
        T: ToString,
    {
        let mut set = IndexSet::<String>::default();

        for (i, seq) in sequences.into_iter().enumerate() {
            let new_set: IndexSet<String> = seq.into_iter().map(|x| x.to_string()).collect();

            if i == 0 {
                set = new_set;
            } else {
                set.retain(|x| new_set.contains(x));
            }
        }

        Self::new(set)
    }

    /// Create new sequence dict.
    fn new(ordering: IndexSet<String>) -> Self {
        Self(ordering)
    }
}
