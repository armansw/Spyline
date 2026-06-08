//! A column-oriented collection of weighted observations.
//!
//! Direct port of the Java `spyline.Points` container. The three input vectors
//! are truncated to the length of the shortest one, exactly like the original.

use crate::point::Point;

/// A set of weighted data points stored as parallel `x`, `y` and `w` vectors.
#[derive(Debug, Clone, Default)]
pub struct Points {
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
}

impl Points {
    /// Build a point set, truncating all three series to the shortest length.
    pub fn new(x: Vec<f64>, y: Vec<f64>, w: Vec<f64>) -> Self {
        let min_len = x.len().min(y.len()).min(w.len());
        Points {
            x: x[..min_len].to_vec(),
            y: y[..min_len].to_vec(),
            w: w[..min_len].to_vec(),
        }
    }

    /// Number of stored points.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Returns `true` when there are no points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Materialise the point at `index` as a [`Point`].
    pub fn item(&self, index: usize) -> Point {
        Point::new(self.x[index], self.y[index], self.w[index])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncates_to_shortest_series() {
        let p = Points::new(vec![1.0, 2.0, 3.0], vec![4.0, 5.0], vec![6.0, 7.0, 8.0]);
        assert_eq!(p.len(), 2);
        assert_eq!(p.item(0), Point::new(1.0, 4.0, 6.0));
        assert_eq!(p.item(1), Point::new(2.0, 5.0, 7.0));
    }

    #[test]
    fn empty_when_any_series_empty() {
        let p = Points::new(vec![], vec![1.0], vec![1.0]);
        assert!(p.is_empty());
    }
}
