//! Weighted least-squares B-spline fitting with optimal knot placement.
//!
//! Faithful port of the Java `spyline.CurveFitter`. The original relied on
//! Apache Commons Math for the Cholesky solve; here that is replaced by the
//! small self-contained [`cholesky_lower`] routine, which returns the same
//! lower-triangular factor `L` such that `A = L * Lᵀ`.

use crate::points::Points;
use crate::spline::Spline;

/// Cholesky decomposition of a symmetric positive-definite matrix `a`.
///
/// Returns the lower-triangular factor `L` (zeros above the diagonal) with
/// `a == L * Lᵀ`, matching `org.apache.commons.math3.linear.CholeskyDecomposition`.
fn cholesky_lower(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    l
}

/// Stateful fitter tracking the data misfit, penalty and combined error.
#[derive(Debug, Clone, Default)]
pub struct CurveFitter {
    m_delta: f64,
    p: f64,
    m_error: f64,
    m_penalty: f64,
}

impl CurveFitter {
    /// Create a fitter seeded with the penalty of `spline`.
    pub fn new(spline: &Spline) -> Self {
        let mut fitter = CurveFitter {
            m_delta: 0.0,
            p: 0.0,
            m_error: 0.0,
            m_penalty: 0.0,
        };
        fitter.penalty(spline);
        fitter
    }

    /// Smoothness penalty: sum of reciprocal interior interval lengths.
    pub fn penalty(&mut self, spline: &Spline) -> f64 {
        let knots = spline.knots();
        let g = spline.internal_knots_num();
        self.m_penalty = 0.0;
        for i in 0..g + 1 {
            self.m_penalty += 1.0 / (knots[(i + 1) as usize] - knots[i as usize]);
        }
        self.m_error = self.m_delta + self.p * self.m_penalty;
        self.m_penalty
    }

    /// Weighted residual sum of squares plus optional smoothing term.
    pub fn delta(&mut self, spline: &Spline, points: &Points, smoothing_weight: f64) -> f64 {
        self.m_delta = 0.0;
        for i in 0..points.len() {
            let pt = points.item(i);
            let e = pt.w() * (pt.y() - spline.value(pt.x()));
            self.m_delta += e * e;
        }
        let mut nu = 0.0;
        if smoothing_weight > 0.0 {
            let k = spline.degree();
            let g = spline.internal_knots_num();
            let coefficients = spline.coefficients();
            for q in k + 1..g + k + 1 {
                let mut e = 0.0;
                for i in q - k - 1..q + 1 {
                    e += coefficients[i as usize] * spline.lead_derivative_difference(i, q);
                    nu += e * e;
                }
            }
            nu *= smoothing_weight;
        }
        self.m_delta += nu;
        self.m_error = self.m_delta + self.p * self.m_penalty;
        self.m_delta
    }

    /// Derivative of the penalty with respect to interior knot `knot_id`.
    pub fn penalty_derivative(spline: &Spline, knot_id: i32) -> f64 {
        let knots = spline.knots();
        let a = knots[(knot_id - 1) as usize];
        let b = knots[knot_id as usize];
        let c = knots[(knot_id + 1) as usize];
        let bma = b - a;
        let cmb = c - b;
        1.0 / (cmb * cmb) - 1.0 / (bma * bma)
    }

    /// Combined error `delta + p * penalty`.
    pub fn error(&mut self, spline: &Spline, points: &Points, sw: f64) -> f64 {
        let d = self.delta(spline, points, sw);
        let pen = self.penalty(spline);
        d + self.p * pen
    }

    /// Gradient of the error with respect to interior knot `knot_id`.
    pub fn error_gradient(
        &self,
        spline: &Spline,
        points: &Points,
        smoothing_weight: f64,
        knot_id: i32,
    ) -> f64 {
        let mut grad_error = 0.0;
        for i in 0..points.len() {
            let pt = points.item(i);
            let w_sq = pt.w() * pt.w();
            let diff = pt.y() - spline.value(pt.x());
            grad_error -=
                spline.value_derivative_knot(pt.x(), knot_id + spline.degree()) * w_sq * diff;
        }

        if self.p > 0.0 {
            grad_error += 0.5 * self.p * Self::penalty_derivative(spline, knot_id);
        }

        if smoothing_weight > 0.0 {
            let mut sm_error = 0.0;
            let k = spline.degree();
            let g = spline.internal_knots_num();
            let coefficients = spline.coefficients();
            for q in k + 1..g + k + 1 {
                let mut sum1 = 0.0;
                let mut sum2 = 0.0;
                for i in q - k - 1..q + 1 {
                    let ci = coefficients[i as usize];
                    let lead_der_diff = spline.lead_derivative_difference(i, q);
                    sum1 += ci * lead_der_diff;
                    sum2 += ci
                        * spline.lead_der_diff_der_knot_with(
                            lead_der_diff,
                            i,
                            q,
                            knot_id + spline.degree(),
                        );
                }
                sm_error += sum1 * sum2;
            }
            grad_error += smoothing_weight * sm_error;
        }

        2.0 * grad_error
    }

    /// Line-search objective: place knots at `fixed + alpha * direction`,
    /// re-approximate, and return the resulting error (or `-1.0` on failure).
    pub fn theta(
        &mut self,
        spline: &mut Spline,
        points: &Points,
        sw: f64,
        alpha: f64,
        direction: &[f64],
        fixed_knots: &[f64],
    ) -> f64 {
        let g = spline.internal_knots_num();
        let mut knots = vec![0.0; (g + 2) as usize];
        knots[0] = spline.left_bound();
        knots[(g + 1) as usize] = spline.right_bound();
        for i in 0..g {
            knots[(i + 1) as usize] = fixed_knots[(i + 1) as usize] + alpha * direction[i as usize];
        }
        spline.set_knots(knots);

        if Self::approximate(spline, points, sw) {
            self.error(spline, points, sw)
        } else {
            -1.0
        }
    }

    /// Squared Euclidean norm of `v`.
    pub fn norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum()
    }

    /// Solve the weighted normal equations for the control coefficients.
    ///
    /// Returns `false` if a sample falls outside the spline support.
    pub fn approximate(spline: &mut Spline, points: &Points, smoothing_weight: f64) -> bool {
        let k = spline.degree();
        let g = spline.internal_knots_num();
        let n = points.len() as i32;
        let dim = (g + k + 1) as usize;
        let mut coefficients = vec![0.0; dim];
        let mut a_mat = vec![vec![0.0; dim]; dim];

        let mut l = 0;
        for r in 0..n {
            let pt = points.item(r as usize);
            l = spline.left_node_index(pt.x(), l);
            if l < 0 {
                return false;
            }
            let b_splines = spline.b_splines(pt.x(), k);
            for i in 0..k + 1 {
                let w_sq = pt.w() * pt.w();
                for j in 0..i + 1 {
                    a_mat[(i + l - k) as usize][(j + l - k) as usize] +=
                        w_sq * b_splines[i as usize] * b_splines[j as usize];
                }
                coefficients[(i + l - k) as usize] += w_sq * pt.y() * b_splines[i as usize];
            }
        }

        if smoothing_weight > 0.0 {
            for q in 0..g {
                for i in q..q + k + 2 {
                    let ai = spline.lead_derivative_difference(i, q + k + 1);
                    for j in q..i + 1 {
                        a_mat[i as usize][j as usize] +=
                            smoothing_weight * ai * spline.lead_derivative_difference(j, q + k + 1);
                    }
                }
            }
        }

        for i in 0..g + k + 1 {
            for j in 0..i {
                a_mat[j as usize][i as usize] = a_mat[i as usize][j as usize];
            }
        }

        let l_mat = cholesky_lower(&a_mat);

        // Forward substitution: solve L y = b.
        for i in 0..g + k + 1 {
            for j in 0..i {
                coefficients[i as usize] -=
                    l_mat[i as usize][j as usize] * coefficients[j as usize];
            }
            coefficients[i as usize] /= l_mat[i as usize][i as usize];
        }

        // Back substitution: solve Lᵀ x = y.
        let mut i = g + k;
        while i >= 0 {
            let mut j = g + k;
            while j >= i + 1 {
                coefficients[i as usize] -=
                    l_mat[j as usize][i as usize] * coefficients[j as usize];
                j -= 1;
            }
            coefficients[i as usize] /= l_mat[i as usize][i as usize];
            i -= 1;
        }

        spline.set_coefficients(coefficients);
        true
    }

    /// One-dimensional minimisation along `direction` (the heart of the
    /// optimal-grid search). Mutates `spline` to the located minimiser.
    pub fn spec_dimensional_minimization(
        &mut self,
        spline: &mut Spline,
        points: &Points,
        sw: f64,
        direction: &[f64],
        error_derivative: &[f64],
        fixed_knots: &[f64],
    ) -> bool {
        let g = spline.internal_knots_num();
        let knots = spline.knots();
        let mut alpha_max = f64::INFINITY;
        let a = spline.left_bound();
        let b = spline.right_bound();

        if direction[0] < 0.0 {
            alpha_max = (a - knots[1]) / direction[0];
        }
        for i in 0..g - 1 {
            if direction[i as usize] > direction[(i + 1) as usize] {
                alpha_max = alpha_max.min(
                    (knots[(i + 2) as usize] - knots[(i + 1) as usize])
                        / (direction[i as usize] - direction[(i + 1) as usize]),
                );
            }
        }
        if direction[(g - 1) as usize] > 0.0 {
            alpha_max = alpha_max.min((b - knots[g as usize]) / direction[(g - 1) as usize]);
        }

        let theta0 = self.m_error;
        let mut theta0_der = 0.0;
        for (d, e) in direction.iter().zip(error_derivative.iter()) {
            theta0_der += d * e;
        }
        let mut alpha0 = 0.0;
        let mut alpha2 = alpha_max / (1.0 - theta0 / alpha_max / theta0_der);
        let mut alpha1 = 0.5 * alpha2;
        let mut q0 = self.m_delta;
        let mut r0 = self.m_penalty;
        let mut theta1 = self.theta(spline, points, sw, alpha1, direction, fixed_knots);
        if theta1 < 0.0 {
            return false;
        }
        let mut q1 = self.m_delta;
        let mut r1 = self.m_penalty;

        let mut iteration = 0;
        let max_num_of_iterations = 10;
        while theta1 >= theta0 && iteration < max_num_of_iterations {
            let alpha_tilde =
                -0.5 * theta0_der * alpha1 * alpha1 / (theta1 - theta0 - theta0_der * alpha1);
            alpha1 = (0.1 * alpha1).max(alpha_tilde);
            theta1 = self.theta(spline, points, sw, alpha1, direction, fixed_knots);
            if theta1 < 0.0 {
                return false;
            }
            q1 = self.m_delta;
            r1 = self.m_penalty;
            iteration += 1;
        }

        if iteration > 0 {
            if theta1 > theta0 {
                self.theta(spline, points, sw, alpha0, direction, fixed_knots);
            }
            return true;
        }

        let mut theta2 = self.theta(spline, points, sw, alpha2, direction, fixed_knots);
        if theta2 < 0.0 {
            return false;
        }
        let mut q2 = self.m_delta;
        let mut r2 = self.m_penalty;

        while theta2 < theta1 {
            alpha0 = alpha1;
            q0 = q1;
            r0 = r1;
            alpha1 = alpha2;
            theta1 = theta2;
            q1 = q2;
            r1 = r2;
            alpha2 = (2.0 * alpha1).min(0.5 * (alpha_max + alpha1));
            theta2 = self.theta(spline, points, sw, alpha2, direction, fixed_knots);
            if theta2 < 0.0 {
                return false;
            }
            q2 = self.m_delta;
            r2 = self.m_penalty;
        }

        let a0 = q0;
        let diff1 = alpha1 - alpha0;
        let diff2 = alpha2 - alpha0;
        let mut a2 = (q1 - q0) / diff1;
        a2 -= (q2 - q0) / diff2;
        a2 /= alpha1 - alpha2;
        let mut a1 = (q1 - a0) / diff1;
        a1 -= a2 * diff1;

        let fraction = diff1 / diff2;
        let numerator = r1 - r0 - fraction * (r2 - r0);
        let temp = ((alpha_max - alpha1) / (alpha_max - alpha0)).ln();
        let denominator = temp - fraction * ((alpha_max - alpha2) / (alpha_max - alpha0)).ln();
        let b2 = numerator / denominator;
        let b1 = (r1 - r0 - b2 * temp) / diff1;

        let a_coef = -2.0 * a2;
        let b_coef = -a_coef * (alpha_max + alpha0) - self.p * b1 - a1;
        let c = (self.p * b1 + a1 + a_coef * alpha0) * alpha_max - self.p * b2;

        let root1 = -0.5 * (b_coef + (b_coef * b_coef - 4.0 * a_coef * c).sqrt()) / a_coef;
        let root2 = -b_coef / a_coef - root1;

        let mut alpha_res = 0.0;
        if 0.0 < root1 && root1 < alpha_max {
            alpha_res = root1;
        } else if 0.0 < root2 && root2 < alpha_max {
            alpha_res = root2;
        }

        let theta_res = self.theta(spline, points, sw, alpha_res, direction, fixed_knots);

        theta_res >= 0.0
    }

    /// Seed the interior knot grid from the data distribution.
    pub fn initiate_grid(spline: &mut Spline, points: &Points) -> bool {
        let k = spline.degree();
        let g = spline.internal_knots_num();
        let mut knots = vec![0.0; (g + 2) as usize];
        knots[0] = spline.left_bound();
        knots[(g + 1) as usize] = spline.right_bound();
        let n = points.len() as i32;

        let mut unique_size = 0;
        let mut index = 0;

        while index < n && points.item(index as usize).x() < knots[(g + 1) as usize] {
            if index != 0
                && points.item(index as usize).x() != points.item((index - 1) as usize).x()
            {
                unique_size += 1;
            }
            index += 1;
        }

        if unique_size <= 0 {
            return false;
        }

        if unique_size < g + k + 1 {
            return false;
        }

        let points_per_knot = unique_size as f64 / (g + 1) as f64;
        let mut knot_index = 1;
        let mut i = 1;
        let mut counter = 0;

        while knot_index < g + 1 {
            while (counter as f64) < knot_index as f64 * points_per_knot
                || points.item(i as usize).x() == points.item((i - 1) as usize).x()
            {
                if points.item(i as usize).x() != points.item((i - 1) as usize).x() {
                    counter += 1;
                }
                i += 1;
            }
            knots[knot_index as usize] =
                0.5 * (points.item(i as usize).x() + points.item((i - 1) as usize).x());
            knot_index += 1;
        }

        spline.set_knots(knots);
        true
    }

    /// Full optimal-grid curve fit: seed the grid, then iterate the
    /// conjugate-gradient knot optimisation until the convergence criteria are
    /// met.
    ///
    /// Note: the original Java never increments its loop counter, so the
    /// search always takes the steepest-descent branch and terminates purely
    /// on the `eps1` / `eps2` criteria. That behaviour is preserved verbatim.
    pub fn approximate_with_optimal_grid(
        &mut self,
        spline: &mut Spline,
        points: &Points,
        smoothing_weight: f64,
        eps1: f64,
        eps2: f64,
    ) -> bool {
        if !Self::initiate_grid(spline, points) {
            return false;
        }

        let g = spline.internal_knots_num();

        if !Self::approximate(spline, points, smoothing_weight) {
            return false;
        }

        let mut direction = vec![0.0; g as usize];
        let mut error_derivative = vec![0.0; g as usize];

        self.delta(spline, points, smoothing_weight);
        self.p = eps1 * self.m_delta * (spline.right_bound() - spline.left_bound())
            / (g + 1) as f64
            / (g + 1) as f64;
        self.m_error = self.m_delta + self.p * self.penalty(spline);

        for i in 0..g {
            error_derivative[i as usize] =
                self.error_gradient(spline, points, smoothing_weight, i + 1);
            direction[i as usize] = -error_derivative[i as usize];
        }

        let mut old_norm = Self::norm(&direction);
        let mut criteria1 = eps1 + eps2;
        let mut criteria2 = criteria1;
        let max_num_of_iter = 1000;

        // Faithful to the Java source: `iteration` is never incremented.
        let iteration = 0;
        let eps2_seq = eps2 * eps2;

        while (criteria1 >= eps1 || criteria2 >= eps2_seq) && iteration < max_num_of_iter {
            let fixed_knots = spline.knots();
            let old_error = self.m_error;
            if !self.spec_dimensional_minimization(
                spline,
                points,
                smoothing_weight,
                &direction,
                &error_derivative,
                &fixed_knots,
            ) {
                spline.set_knots(fixed_knots);
                return Self::approximate(spline, points, smoothing_weight);
            }

            for i in 0..g {
                error_derivative[i as usize] =
                    self.error_gradient(spline, points, smoothing_weight, i + 1);
            }

            let new_norm = Self::norm(&error_derivative);

            if iteration % g == 0 {
                for i in 0..g {
                    direction[i as usize] = -error_derivative[i as usize];
                }
            } else {
                let temp = new_norm / old_norm;
                for i in 0..g {
                    direction[i as usize] *= temp;
                    direction[i as usize] -= error_derivative[i as usize];
                }
            }

            let mut numerator = 0.0;
            let mut denominator = 0.0;

            let knots = spline.knots();
            for i in 1..knots.len() - 1 {
                let temp = knots[i] - fixed_knots[i];
                numerator += temp * temp;
                denominator += fixed_knots[i] * fixed_knots[i];
            }

            criteria1 = (old_error - self.m_error).abs() / old_error;
            criteria2 = numerator / denominator;

            old_norm = new_norm;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_reconstructs_matrix() {
        // A = [[4,12,-16],[12,37,-43],[-16,-43,98]] (classic SPD example).
        let a = vec![
            vec![4.0, 12.0, -16.0],
            vec![12.0, 37.0, -43.0],
            vec![-16.0, -43.0, 98.0],
        ];
        let l = cholesky_lower(&a);
        let n = a.len();
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += l[i][k] * l[j][k];
                }
                assert!((s - a[i][j]).abs() < 1e-9, "mismatch at ({i},{j})");
            }
        }
    }

    #[test]
    fn norm_is_sum_of_squares() {
        assert_eq!(CurveFitter::norm(&[3.0, 4.0]), 25.0);
    }
}
