use std::{cmp, io};

use crate::{ChromPos, ChromDict};

/// Iterator over intersecting positions from different sources.
///
/// Merging assumes that positions from all sources follow the same ordering. That is,
/// the subset of chromosomes that appear in all sources must occur in the same order
/// in each source; within each chromosome, positions must occur in ascending order.
///
/// The intersection of chromosomes and their ordering must be pre-calculated from
/// from header information or otherwise and stored in a [`ChromDict`]. See in particular
/// [`ChromDict::from_merged_chromosomes`](ChromDict::from_merged_chromosomes)
/// for help.
pub struct Merge<I, T>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos
{
    iters: Vec<Search<I, T>>,
    dict: ChromDict
}

impl<I, T> Merge<I, T>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos
{
    /// Advance all iterators to the next position on a chromosome in the chromosome dict.
    fn all_next_in_dict(&mut self) -> Option<io::Result<MultiChromPos<T>>> {
        let dict = &self.dict;

        self.iters
            .iter_mut()
            .map(|x| x.next_in_dict(dict))
            .collect::<Option<io::Result<Vec<T>>>>()
            .map(|x| x.map(MultiChromPos))
    }

    /// Create a new merge iterator.
    ///
    /// The intersection of chromosomes and their ordering must be pre-calculated
    /// and stored in a [`ChromDict`].
    /// See in particular [`ChromDict::from_merged_chromosomes`](ChromDict::from_merged_chromosomes)
    /// for help.
    pub fn new(input: Vec<I>, dict: ChromDict) -> Self {
        Self {
            iters: input.into_iter().map(|x| Search::new(x)).collect(),
            dict,
        }
    }
}

impl<I, T> Iterator for Merge<I, T>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos
{
    type Item = io::Result<Vec<T>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut positions = match self.all_next_in_dict()? {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        };

        while !positions.is_intersection() {
            // Find the greatest position
            let argmax = positions.argmax(&self.dict)?;

            // Find the first position in each iterator at least as great as the current greatest
            // (The awkward indexing here is required to appease the borrow checker)
            for i in 0..positions.0.len() {
                if !positions.0[i].colocated(&positions.0[argmax]) {
                    let target = &positions.0[argmax];
                    positions.0[i] = match self.iters[i].search(target, &self.dict)? {
                        Ok(v) => v,
                        Err(e) => return Some(Err(e)),
                    };
                }
            }
        }

        Some(Ok(positions.0))
    }
}

/// Multiple positions that may or may not be intersecting.
struct MultiChromPos<T>(Vec<T>) where T: ChromPos;

impl<T> MultiChromPos<T>
where
    T: ChromPos,
{
    /// Check if all positions are colocated.
    fn is_intersection(&self) -> bool {
        let first = &self.0[0];

        self.0.iter().skip(1).all(|x| x.colocated(first))
    }

    /// Returns the index of the greatest position.
    pub fn argmax(&self, dict: &ChromDict) -> Option<usize> {
        let mut argmax = 0;

        for (i, position) in self.0.iter().enumerate().skip(1) {
            match dict.compare(position, &self.0[argmax]) {
                Some(cmp::Ordering::Greater) => argmax = i,
                Some(cmp::Ordering::Equal) => (),
                Some(cmp::Ordering::Less) => (),
                None => return None,
            }
        }

        Some(argmax)
    }
}

/// Search struct for finding positions in ordered iterators.
struct Search<I, T>(I)
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos
;

impl<I, T> Search<I, T>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos
{
    /// Create new position search iterator.
    pub fn new(inner: I) -> Self {
        Self(inner)
    }

    /// Search to next position on a chromosome contained in chromosome dictionary.
    ///
    /// If iterator is exhausted before such a position is found, returns None.
    fn next_in_dict(&mut self, dict: &ChromDict) -> Option<io::Result<T>> {
        while let Some(v) = self.0.next() {
            match v {
                Ok(v) => {
                    if dict.contains(&v) {
                        return Some(Ok(v));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    /// Search for a specific position in iterator.
    ///
    /// Returns target position if found, otherwise returns the first position that is greater than
    /// the target position (relative to chromosome dictionary). If iterator is exhaused before
    /// finding a position equal to or greater than the target, returns None.
    pub fn search(&mut self, target: &T, dict: &ChromDict) -> Option<io::Result<T>> {
        while let Some(v) = self.next_in_dict(dict) {
            match v {
                Ok(v) => match dict.compare(&v, target) {
                    Some(cmp::Ordering::Equal) | Some(cmp::Ordering::Greater) => {
                        return Some(Ok(v))
                    }
                    Some(cmp::Ordering::Less) => continue,
                    None => return None,
                },
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn next_in_dict() {
        let positions = vec![("1", 1), ("1", 2), ("2", 1), ("2", 3), ("4", 2), ("5", 1)];

        let dict = ChromDict::from_chromosomes(vec!["2", "4"]);

        let mut search = Search::new(positions.into_iter().map(|x| Ok(x)));

        assert!(matches!(search.next_in_dict(&dict), Some(Ok(("2", 1)))));
        assert!(matches!(search.next_in_dict(&dict), Some(Ok(("2", 3)))));
        assert!(matches!(search.next_in_dict(&dict), Some(Ok(("4", 2)))));
        assert!(matches!(search.next_in_dict(&dict), None));
    }

    #[test]
    fn search() {
        let positions = vec![("1", 1), ("1", 2), ("2", 1), ("2", 3), ("4", 2), ("5", 1)];

        let dict = ChromDict::from_chromosomes(vec!["2", "4"]);

        let mut iter = Search::new(positions.into_iter().map(|x| Ok(x)));

        assert!(matches!(iter.search(&("2", 1), &dict), Some(Ok(("2", 1)))));
        assert!(matches!(iter.search(&("2", 2), &dict), Some(Ok(("2", 3)))));
        assert!(matches!(iter.search(&("4", 1), &dict), Some(Ok(("4", 2)))));
        assert!(matches!(iter.search(&("4", 3), &dict), None));
    }
}
