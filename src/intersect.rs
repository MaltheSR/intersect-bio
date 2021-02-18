use std::{
    cmp, io,
    ops::{Index, IndexMut},
};

use crate::{ChromDict, ChromPos};

/// Intersect iterator.
///
/// An iterator over the intersection of positions in pre-sorted files, where a position
/// is anything that implements [`ChromPos`]. Merging requires that a chromosome dictionary
/// is computed ahead of time. See [`ChromDict`] for details.
pub struct Intersect<I> {
    iters: Vec<Search<I>>,
    dict: ChromDict,
}

impl<I> Intersect<I> {
    /// Create new intersect iterator.
    pub fn new(input: Vec<I>, dict: ChromDict) -> Self {
        Self {
            iters: input.into_iter().map(|x| Search::new(x)).collect(),
            dict,
        }
    }
}

impl<I, T> Intersect<I>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos,
{
    /// Find next candidate positions.
    ///
    /// A candidate position is any position located on any of the chromosomes contained
    /// in the current chromosome dictionary; if a position is not on such a chromosome,
    /// it cannot be part of an intersection.
    fn next_candidates(&mut self) -> Option<io::Result<Positions<T>>> {
        let dict = &self.dict;

        self.iters
            .iter_mut()
            .map(|x| x.next_candidate(dict))
            .collect::<Option<io::Result<Vec<T>>>>()
            .map(|x| x.map(Positions))
    }
}

impl<I, T> Iterator for Intersect<I>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos,
{
    type Item = io::Result<Vec<T>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut positions = match self.next_candidates()? {
            Ok(v) => v,
            Err(e) => return Some(Err(e)),
        };

        let n = positions.len();

        while !positions.is_intersection() {
            // Find the max position, and forward all iterators currently at a position less than or
            // equal to max to the first position greater than or equal to max (awkward indexing is
            // required to appease borrow checker)
            let argmax = positions.argmax(&self.dict)?;

            for i in (0..argmax).chain(argmax + 1..n) {
                let max = &positions[argmax];

                if !positions[i].intersect(max) {
                    positions[i] = match self.iters[i].search(max, &self.dict)? {
                        Ok(v) => v,
                        Err(e) => return Some(Err(e)),
                    };
                }
            }
        }

        Some(Ok(positions.0))
    }
}

/// Multiple positions.
///
/// Helper newtype for a collection of positions that may or may not be intersecting.
struct Positions<T>(Vec<T>);

impl<T> Positions<T>
where
    T: ChromPos,
{
    /// Get number of positions.
    fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if all positions intersect.
    fn is_intersection(&self) -> bool {
        let first = &self.0[0];

        self.0.iter().skip(1).all(|x| x.intersect(first))
    }

    /// Get index of the greatest position.
    ///
    /// If all positions are located on chromosomes contained in chromosome dictionary,
    /// returns the index of the positions with the greatest position. Otherwise, returns
    /// `None`. If multiple positions are tied for greatest, returns the first of these.
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

impl<T> Index<usize> for Positions<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> IndexMut<usize> for Positions<T> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

/// Search iterator.
///
/// Helper newtype for position iterators to search forward for positions meeting particular
/// criteria.
struct Search<I>(I);

impl<I> Search<I> {
    /// Create new search iterator.
    pub fn new(inner: I) -> Self {
        Self(inner)
    }
}

impl<I, T> Search<I>
where
    I: Iterator<Item = io::Result<T>>,
    T: ChromPos,
{
    /// Find next candidate position.
    ///
    /// A candidate position, relative to some chromosome dictionary, is any position located on
    /// a chromosome contained in the dictionary. If the iterator is exhausted before such a
    /// position is found, returns None.
    fn next_candidate(&mut self, dict: &ChromDict) -> Option<io::Result<T>> {
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

    /// Search for target position.
    ///
    /// Returns target position if found, otherwise returns the first position that is greater than
    /// the target position, relative to chromosome dictionary. If iterator is exhausted before
    /// finding a position equal to or greater than the target, returns None.
    pub fn search(&mut self, target: &T, dict: &ChromDict) -> Option<io::Result<T>> {
        while let Some(v) = self.next_candidate(dict) {
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

    fn mock_source<'a>(v: Vec<(&'a str, u32)>) -> impl Iterator<Item = io::Result<(&'a str, u32)>> {
        v.into_iter().map(|x| Ok(x))
    }

    fn mock_input<'a>(
        vs: Vec<Vec<(&'a str, u32)>>,
    ) -> Vec<impl Iterator<Item = io::Result<(&'a str, u32)>>> {
        vs.into_iter().map(|x| mock_source(x)).collect()
    }

    #[test]
    fn intersect() {
        let dict = ChromDict::from_ids(vec!["2", "4"]);

        let input = mock_input(vec![
            vec![("1", 1), ("1", 2), ("2", 1), ("2", 3), ("4", 1)],
            vec![
                ("1", 1),
                ("1", 2),
                ("2", 2),
                ("2", 3),
                ("4", 1),
                ("4", 5),
                ("5", 1),
            ],
            vec![("2", 1), ("2", 2), ("2", 3), ("3", 1), ("4", 1), ("4", 7)],
        ]);

        let mut intersect = Intersect::new(input, dict);

        assert_eq!(
            intersect.next().unwrap().unwrap(),
            vec![("2", 3), ("2", 3), ("2", 3)]
        );
        assert_eq!(
            intersect.next().unwrap().unwrap(),
            vec![("4", 1), ("4", 1), ("4", 1)]
        );
        assert!(matches!(intersect.next(), None));
    }

    #[test]
    fn positions_intersect() {
        let mut positions = Positions(vec![("1", 1), ("1", 1), ("1", 1), ("1", 1), ("1", 1)]);
        assert!(positions.is_intersection());

        positions.0[0] = ("1", 2);
        assert!(!positions.is_intersection());

        positions.0[0] = ("2", 1);
        assert!(!positions.is_intersection());
    }

    #[test]
    fn positions_argmax() {
        let dict = ChromDict::from_ids(vec!["1", "2"]);

        let mut positions = Positions(vec![("1", 1), ("1", 2), ("1", 5), ("1", 1), ("1", 3)]);
        assert_eq!(positions.argmax(&dict), Some(2));

        positions.0[1] = ("1", 5);
        assert_eq!(positions.argmax(&dict), Some(1));

        positions.0[4] = ("2", 1);
        assert_eq!(positions.argmax(&dict), Some(4));

        positions.0[4] = ("3", 1);
        assert_eq!(positions.argmax(&dict), None);
    }

    #[test]
    fn search_candidate() {
        let positions = vec![("1", 1), ("1", 2), ("2", 1), ("2", 3), ("4", 2), ("5", 1)];

        let dict = ChromDict::from_ids(vec!["2", "4"]);

        let mut search = Search::new(positions.into_iter().map(|x| Ok(x)));

        assert_eq!(search.next_candidate(&dict).unwrap().unwrap(), ("2", 1));
        assert_eq!(search.next_candidate(&dict).unwrap().unwrap(), ("2", 3));
        assert_eq!(search.next_candidate(&dict).unwrap().unwrap(), ("4", 2));
        assert!(matches!(search.next_candidate(&dict), None));
    }

    #[test]
    fn search_position() {
        let positions = vec![("1", 1), ("1", 2), ("2", 1), ("2", 3), ("4", 2), ("5", 1)];

        let dict = ChromDict::from_ids(vec!["2", "4"]);

        let mut iter = Search::new(positions.into_iter().map(|x| Ok(x)));

        assert_eq!(iter.search(&("2", 1), &dict).unwrap().unwrap(), ("2", 1));
        assert_eq!(iter.search(&("2", 2), &dict).unwrap().unwrap(), ("2", 3));
        assert_eq!(iter.search(&("4", 1), &dict).unwrap().unwrap(), ("4", 2));
        assert!(matches!(iter.search(&("4", 3), &dict), None));
    }
}
