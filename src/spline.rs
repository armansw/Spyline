//! B-spline model: evaluation, derivatives and knot manipulation.
//!
//! Faithful port of the Java `spyline.Spline` class. The knot vector uses the
//! usual clamped layout: the degree is repeated `k + 1` times at each boundary
//! with `g` interior knots in between. Integer index arithmetic mirrors the
//! original exactly, so a handful of helpers ([`Spline::kn`], [`Spline::co`],
//! …) wrap the signed-index access used throughout the algorithms.

/// Recursive factorial matching the Java implementation: `factorial(n) == n`
/// for `n <= 2`, so `factorial(0) == 0` and `factorial(3) == 6`.
fn factorial(n: i32) -> i32 {
    if n <= 2 {
        n
    } else {
        n * factorial(n - 1)
    }
}

fn sort_ascending(v: &mut [f64]) {
    v.sort_by(|a, b| {
        a.partial_cmp(b)
            .expect("knot coordinates must be comparable")
    });
}

/// A weighted B-spline of fixed degree with a mutable knot grid.
#[derive(Debug, Clone)]
pub struct Spline {
    k: i32,
    k_fact: i32,
    coefficients: Vec<f64>,
    g: i32,
    a: f64,
    b: f64,
    knots: Vec<f64>,
}

impl Spline {
    // --- signed-index accessors (mirror Java list `get`/`set`) -------------

    #[inline]
    fn kn(&self, i: i32) -> f64 {
        self.knots[i as usize]
    }

    #[inline]
    fn set_kn(&mut self, i: i32, v: f64) {
        self.knots[i as usize] = v;
    }

    #[inline]
    fn co(&self, i: i32) -> f64 {
        self.coefficients[i as usize]
    }

    #[inline]
    fn set_co(&mut self, i: i32, v: f64) {
        self.coefficients[i as usize] = v;
    }

    // --- construction ------------------------------------------------------

    /// Build a spline of degree `spline_degree` from `coefficients` and a list
    /// of `horizontal_knots` (boundary + interior knots, any order).
    pub fn new(coefficients: Vec<f64>, mut horizontal_knots: Vec<f64>, spline_degree: i32) -> Self {
        let k = spline_degree;
        let k_fact = factorial(k);
        sort_ascending(&mut horizontal_knots);
        let g = horizontal_knots.len() as i32 - 2;
        let a = horizontal_knots[0];
        let b = horizontal_knots[(g + 1) as usize];

        let mut knots = Vec::new();
        for _ in 0..k + 1 {
            knots.push(a);
        }
        for _ in 0..g {
            knots.push(0.0);
        }
        for _ in 0..k + 1 {
            knots.push(b);
        }

        let mut spline = Spline {
            k,
            k_fact,
            coefficients,
            g,
            a,
            b,
            knots,
        };
        for i in 0..g {
            spline.set_kn(i + k + 1, horizontal_knots[(i + 1) as usize]);
        }
        spline
    }

    // --- queries -----------------------------------------------------------

    /// Index of the knot interval containing `point`, searching from `min_id`.
    /// Returns `-1` when `point` lies outside `[a, b]`.
    pub fn left_node_index(&self, point: f64, min_id: i32) -> i32 {
        if point < self.a || point > self.b {
            return -1;
        }
        let mut l = min_id;
        while l < self.g + self.k && (self.kn(l) > point || self.kn(l + 1) <= point) {
            l += 1;
        }
        l
    }

    /// Number of interior knots `g`.
    pub fn internal_knots_num(&self) -> i32 {
        self.g
    }

    /// Spline degree `k`.
    pub fn degree(&self) -> i32 {
        self.k
    }

    /// Left boundary `a`.
    pub fn left_bound(&self) -> f64 {
        self.a
    }

    /// Right boundary `b`.
    pub fn right_bound(&self) -> f64 {
        self.b
    }

    /// `k!` precomputed at construction.
    pub fn degree_factorial(&self) -> i32 {
        self.k_fact
    }

    /// Borrow the control coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Replace the control coefficients.
    pub fn set_coefficients(&mut self, coefficients: Vec<f64>) {
        self.coefficients = coefficients;
    }

    /// Override the left edge.
    pub fn set_left_edge(&mut self, left_edge: f64) {
        self.a = left_edge;
    }

    /// Override the right edge.
    pub fn set_right_edge(&mut self, right_edge: f64) {
        self.b = right_edge;
    }

    /// Return the boundary + interior knots (the "horizontal" knot grid).
    pub fn knots(&self) -> Vec<f64> {
        let mut sliced = Vec::new();
        for i in self.k..self.k + self.g + 2 {
            sliced.push(self.kn(i));
        }
        sliced
    }

    /// Rebuild the internal clamped knot vector from a horizontal knot grid.
    pub fn set_knots(&mut self, mut horizontal_knots: Vec<f64>) {
        sort_ascending(&mut horizontal_knots);
        self.g = horizontal_knots.len() as i32 - 2;
        self.a = horizontal_knots[0];
        self.b = horizontal_knots[(self.g + 1) as usize];

        let mut local_knots = Vec::new();
        for _ in 0..self.k + 1 {
            local_knots.push(self.a);
        }
        for _ in 0..self.g {
            local_knots.push(0.0);
        }
        for _ in 0..self.k + 1 {
            local_knots.push(self.b);
        }
        self.knots = local_knots;
        for i in 0..self.g {
            self.set_kn(i + self.k + 1, horizontal_knots[(i + 1) as usize]);
        }
    }

    // --- basis functions ---------------------------------------------------

    /// Evaluate a single B-spline basis function of degree `deg` anchored at
    /// `knot_id`.
    pub fn b_spline(&self, point: f64, deg: i32, knot_id: i32) -> f64 {
        if point < self.kn(knot_id) || point > self.kn(knot_id + deg + 1) {
            return 0.0;
        }

        if deg == 0 {
            return if point != self.kn(knot_id + deg + 1) {
                1.0
            } else {
                0.0
            };
        }

        if self.kn(knot_id + deg) < self.kn(knot_id + deg + 1) {
            let mut j = 0;
            while j < deg && self.kn(knot_id + j) == self.kn(knot_id + j + 1) {
                j += 1;
            }
            if j == deg {
                return ((self.kn(knot_id + deg + 1) - point)
                    / (self.kn(knot_id + deg + 1) - self.kn(knot_id)))
                .powi(deg);
            }
        }

        let l = self.left_node_index(point, knot_id);
        let mut buff = vec![0.0; (deg + 1) as usize];
        buff[(knot_id - 1 + deg) as usize] = 1.0;

        for j in 1..deg + 1 {
            let mut i = l;
            while i >= l - deg + j {
                let alpha = (point - self.kn(i)) / (self.kn(i + 1 + deg - j) - self.kn(i));
                buff[(i - l + deg) as usize] = alpha * buff[(i - l + deg) as usize]
                    + (1.0 - alpha) * buff[(i - 1 - l + deg) as usize];
                i -= 1;
            }
        }

        buff[deg as usize]
    }

    /// Evaluate all non-zero B-spline basis functions of degree `deg` at
    /// `point` (stable De Boor recursion).
    pub fn b_splines(&self, point: f64, deg: i32) -> Vec<f64> {
        let mut buff = vec![0.0; (self.k + 1) as usize];
        if deg > self.k {
            return buff;
        }

        let l = self.left_node_index(point, 0);
        buff[deg as usize] = 1.0;

        for r in 1..deg + 1 {
            let mut v = l - r + 1;
            let mut w2 = (self.kn(v + r) - point) / (self.kn(v + r) - self.kn(v));
            buff[(deg - r) as usize] = w2 * buff[(deg - r + 1) as usize];
            for i in deg - r + 1..deg {
                let w1 = w2;
                v += 1;
                w2 = (self.kn(v + r) - point) / (self.kn(v + r) - self.kn(v));
                buff[i as usize] = (1.0 - w1) * buff[i as usize] + w2 * buff[(i + 1) as usize];
            }
            buff[deg as usize] = (1.0 - w2) * buff[deg as usize];
        }
        buff
    }

    /// Insert a new interior knot at integer `coordinate`, updating coefficients.
    pub fn insert_node(&mut self, coordinate: i32) {
        let j = self.left_node_index(coordinate as f64, 0);

        if j < 0 || self.kn(j) == coordinate as f64 {
            return;
        }

        self.knots.push(coordinate as f64);
        sort_ascending(&mut self.knots);
        self.coefficients.push(self.co(self.g + self.k));

        let mut i = self.g + self.k;
        while i >= j + 1 {
            self.set_co(i, self.co(i - 1));
            i -= 1;
        }

        let mut i = j;
        while i >= j - self.k + 1 {
            let ri = (coordinate as f64 - self.kn(i)) / (self.kn(i + self.k + 1) - self.kn(i));
            self.set_co(i, ri * self.co(i) + (1.0 - ri) * self.co(i - 1));
            i -= 1;
        }

        self.g += 1;
    }

    /// `der_degree`-th derivative of a single basis function (recursive).
    pub fn b_spline_derivative(&self, point: f64, l: i32, i: i32, der_degree: i32) -> f64 {
        if der_degree == 0 {
            return self.b_spline(point, l, i);
        }

        if l == 0 {
            return 0.0;
        }

        let mut spline = 0.0;
        let c1 = self.kn(i + l) - self.kn(i);
        let c2 = self.kn(i + l + 1) - self.kn(i + 1);
        if c1 != 0.0 {
            spline += self.b_spline_derivative(point, l - 1, i, der_degree - 1) / c1;
        }
        if c2 != 0.0 {
            spline -= self.b_spline_derivative(point, l - 1, i + 1, der_degree - 1) / c2;
        }
        l as f64 * spline
    }

    /// Evaluate the spline at `point` via De Boor's algorithm.
    pub fn value(&self, point: f64) -> f64 {
        if point < self.a || point > self.b {
            return 0.0;
        }

        let l = self.left_node_index(point, 0);
        if l < 0 {
            return 0.0;
        }

        let mut buff = vec![0.0; (self.k + 1) as usize];
        for i in 0..self.k + 1 {
            buff[i as usize] = self.co(i + l - self.k);
        }

        for j in 1..self.k + 1 {
            let mut i = l;
            while i >= l - self.k + j {
                let alpha = (point - self.kn(i)) / (self.kn(i + 1 + self.k - j) - self.kn(i));
                buff[(i - l + self.k) as usize] = alpha * buff[(i - l + self.k) as usize]
                    + (1.0 - alpha) * buff[(i - 1 - l + self.k) as usize];
                i -= 1;
            }
        }

        buff[self.k as usize]
    }

    /// `der_degree`-th derivative of the spline value at `point`.
    pub fn value_derivative(&self, point: f64, der_degree: i32) -> f64 {
        if der_degree == 0 {
            return self.value(point);
        }

        if der_degree > self.k {
            return 0.0;
        }

        let l = self.left_node_index(point, 0);

        if l < 0 {
            if der_degree > 1 {
                return 0.0;
            }
            return -1.0;
        }

        let mut alpha = 1.0;
        let mut spline = 0.0;

        for i in 0..der_degree {
            alpha *= (self.k - i) as f64;
        }

        let mut buff = vec![0.0; (self.k + 1) as usize];
        for i in 0..self.k + 1 {
            buff[i as usize] = self.co(i + l - self.k);
        }

        for j in 1..der_degree {
            let mut i = j;
            while i >= l - self.k + j {
                buff[(i - l + self.k) as usize] = (buff[(i - l + self.k) as usize]
                    - buff[(i - l - 1 + self.k) as usize])
                    / (self.kn(i + 1 + self.k - j) - self.kn(i));
                i -= 1;
            }
        }

        for i in der_degree..self.k + 1 {
            spline += buff[i as usize] * self.b_spline(point, self.k - der_degree, l + i - self.k);
        }

        alpha * spline
    }

    /// Leading derivative jump of basis `i` at knot `q`.
    pub fn lead_derivative_difference(&self, i: i32, q: i32) -> f64 {
        if i < q - self.k - 1 || i > q {
            return 0.0;
        }
        let numerator = (2 * (self.k % 2) - 1) as f64
            * self.k_fact as f64
            * (self.kn(i + self.k + 1) - self.kn(i));
        let mut denominator = 1.0;
        for j in i..i + self.k + 2 {
            if j != q {
                denominator *= self.kn(q) - self.kn(j);
            }
        }
        numerator / denominator
    }

    /// Derivative of the leading-derivative jump w.r.t. knot `l`, given the
    /// precomputed jump `lddk`.
    pub fn lead_der_diff_der_knot_with(&self, lddk: f64, i: i32, q: i32, l: i32) -> f64 {
        if l < i || l > i + self.k + 1 {
            return 0.0;
        }

        if l != i && l != q && l != i + self.k + 1 {
            return lddk / (self.kn(q) - self.kn(l));
        }

        let c = (2 * (self.k % 2) - 1) as f64 * self.k_fact as f64;
        let mut product = c / lddk;
        let mut total_sum = 0.0;

        if q != i && q != i + self.k + 1 {
            if l == i {
                let temp =
                    (self.kn(i + self.k + 1) - self.kn(i)) * (self.kn(i + self.k + 1) - self.kn(q));
                return lddk / (self.kn(q) - self.kn(i) / temp);
            }
            if l == i + self.k + 1 {
                let temp = (self.kn(i + self.k + 1) - self.kn(i)) * (self.kn(q) - self.kn(i));
                return lddk / (self.kn(q) - self.kn(i + self.k + 1)) / temp;
            }
            if q == l {
                product *= self.kn(i + self.k + 1) - self.kn(i);
                for j in i..i + self.k + 2 {
                    if j != q {
                        total_sum += product / (self.kn(q) - self.kn(j));
                    }
                }
                product *= product;
                return -c * (self.kn(i + self.k + 1) - self.kn(i)) * total_sum / product;
            }
        } else {
            for j in i + 1..i + self.k + 1 {
                total_sum += product / (self.kn(q) - self.kn(j));
            }
            return -c * total_sum / (product * product);
        }
        0.0
    }

    /// Convenience wrapper computing the leading-derivative jump internally.
    pub fn lead_der_diff_der_knot(&self, i: i32, q: i32, l: i32) -> f64 {
        if l < i || l > i + self.k + 1 {
            return 0.0;
        }
        self.lead_der_diff_der_knot_with(self.lead_derivative_difference(i, q), i, q, l)
    }

    /// Derivative of the spline value at `point` w.r.t. the knot `knot_id`.
    pub fn value_derivative_knot(&self, point: f64, knot_id: i32) -> f64 {
        if point < self.a {
            if knot_id == self.k + 1 {
                let numerator = -(self.k as f64) * (self.co(1) - self.co(0)) * (point - self.a);
                let denominator = (self.kn(knot_id) - self.a) * (self.kn(knot_id) - self.a);
                return numerator / denominator;
            }
            return 0.0;
        }

        if point > self.b {
            if knot_id == self.g + self.k {
                let numerator = -(self.k as f64)
                    * (self.co(self.g + self.k) - self.co(self.g + self.k - 1))
                    * (point - self.b);
                let denominator = (self.kn(knot_id) - self.b) * (self.kn(knot_id) - self.b);
                return numerator / denominator;
            }
            return 0.0;
        }

        if point <= self.kn(knot_id - self.k) || point >= self.kn(knot_id + self.k) {
            return 0.0;
        }

        let mut l = self.left_node_index(point, 0);
        if l < 0 {
            return 0.0;
        }

        if l >= knot_id {
            l += 1;
        }

        let mut buff = vec![0.0; (self.k + 1) as usize];

        for i in 0..self.k + 1 {
            if i < knot_id - l || i > knot_id - l + self.k {
                buff[i as usize] = 0.0;
            } else {
                let mut val = self.co(i + l - self.k - 1) - self.co(i + l - self.k);
                if i + l + 1 <= knot_id {
                    val /= self.kn(i + l + 1) - self.kn(i + l - self.k);
                } else if i <= knot_id {
                    val /= self.kn(i + l) - self.kn(i + l - self.k);
                } else {
                    val /= self.kn(i + l) - self.kn(i + l - self.k - 1);
                }
                buff[i as usize] = val;
            }
        }

        for j in 1..self.k + 1 {
            let mut i = l;
            while i >= l - self.k + j {
                let alpha = if i + 1 + self.k - j <= knot_id {
                    (point - self.kn(i)) / (self.kn(i + 1 + self.k - j) - self.kn(i))
                } else if i <= knot_id {
                    (point - self.kn(i)) / (self.kn(i + self.k - j) - self.kn(i))
                } else {
                    (point - self.kn(i - 1)) / (self.kn(i + self.k - j) - self.kn(i - 1))
                };
                buff[(i - l + self.k) as usize] = alpha * buff[(i - l + self.k) as usize]
                    + (1.0 - alpha) * buff[(i - 1 - l + self.k) as usize];
                i -= 1;
            }
        }
        buff[self.k as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_spline() -> Spline {
        let knots = vec![1.0, 4.0, 8.0, 12.0, 16.0, 39039.0, 17.0];
        let coefficients = vec![0.0; knots.len() - 2];
        Spline::new(coefficients, knots, 3)
    }

    #[test]
    fn factorial_matches_java_semantics() {
        assert_eq!(factorial(0), 0);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(2), 2);
        assert_eq!(factorial(3), 6);
        assert_eq!(factorial(4), 24);
    }

    #[test]
    fn bounds_and_degree_are_exposed() {
        let s = sample_spline();
        assert_eq!(s.degree(), 3);
        assert_eq!(s.degree_factorial(), 6);
        assert_eq!(s.left_bound(), 1.0);
        assert_eq!(s.right_bound(), 39039.0);
    }

    #[test]
    fn outside_bounds_returns_zero() {
        let s = sample_spline();
        assert_eq!(s.value(-5.0), 0.0);
        assert_eq!(s.value(50_000.0), 0.0);
        assert_eq!(s.left_node_index(-5.0, 0), -1);
    }
}
