//! A single weighted observation `(x, y)` with weight `w`.
//!
//! Direct port of the Java `spyline.Point` value type.

/// A weighted data point used by the curve fitter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    x: f64,
    y: f64,
    w: f64,
}

impl Point {
    /// Create a new point with abscissa `x`, ordinate `y` and weight `w`.
    pub fn new(x: f64, y: f64, w: f64) -> Self {
        Point { x, y, w }
    }

    /// The weight of the point.
    pub fn w(&self) -> f64 {
        self.w
    }

    /// The abscissa (x coordinate).
    pub fn x(&self) -> f64 {
        self.x
    }

    /// The ordinate (y coordinate).
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Override the weight.
    pub fn set_w(&mut self, w: f64) {
        self.w = w;
    }

    /// Override the abscissa.
    pub fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    /// Override the ordinate.
    pub fn set_y(&mut self, y: f64) {
        self.y = y;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getters_return_constructed_values() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
        assert_eq!(p.w(), 3.0);
    }

    #[test]
    fn setters_mutate_in_place() {
        let mut p = Point::new(0.0, 0.0, 0.0);
        p.set_x(4.0);
        p.set_y(5.0);
        p.set_w(6.0);
        assert_eq!(p, Point::new(4.0, 5.0, 6.0));
    }
}
