//! End-to-end tests exercising the full fitting pipeline on the reference
//! dataset used by the command-line driver.

use spyline::curve_fitter::CurveFitter;
use spyline::points::Points;
use spyline::spline::Spline;

fn reference_dataset() -> Points {
    let x = vec![
        1.0, 3.0, 7.0, 15.0, 31.0, 63.0, 127.0, 255.0, 383.0, 639.0, 895.0, 1151.0, 1663.0,
        2175.0, 3199.0, 4223.0, 6271.0, 10367.0, 14463.0, 22655.0, 39039.0,
    ];
    let y = vec![
        900.0, 779.0, 692.0, 617.0, 551.0, 476.0, 387.0, 283.0, 255.0, 221.0, 200.0, 185.0,
        164.0, 150.0, 81.0, 65.0, 46.0, 26.0, 16.0, 7.0, 1.0,
    ];
    let w = vec![1.0; x.len()];
    Points::new(x, y, w)
}

fn build_spline() -> Spline {
    let mut knots = vec![1.0, 4.0, 8.0, 12.0, 16.0, 39039.0];
    let coefficients = vec![0.0; knots.len() - 2];
    let penultimate = knots[knots.len() - 2];
    knots.push(penultimate + 1.0);
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Spline::new(coefficients, knots, 3)
}

#[test]
fn uniform_grid_fit_is_finite_and_bounded() {
    let points = reference_dataset();
    let mut spline = build_spline();

    assert!(CurveFitter::initiate_grid(&mut spline, &points));
    assert!(CurveFitter::approximate(&mut spline, &points, 1e-9));

    for i in 0..points.len() {
        let v = spline.value(points.item(i).x());
        assert!(v.is_finite(), "spline value must be finite");
        assert!((-100.0..=1100.0).contains(&v), "value {v} out of range");
    }
}

#[test]
fn optimal_grid_keeps_knots_sorted_and_reduces_error() {
    let points = reference_dataset();
    let mut spline = build_spline();
    let mut fitter = CurveFitter::new(&spline);

    CurveFitter::initiate_grid(&mut spline, &points);
    CurveFitter::approximate(&mut spline, &points, 1e-9);
    let uniform_error = fitter.error(&spline, &points, 1e-9);

    assert!(fitter.approximate_with_optimal_grid(&mut spline, &points, 1e-9, 1e-3, 1e-3));
    let optimal_error = fitter.error(&spline, &points, 1e-9);

    // The optimal grid should not be worse than the uniform one.
    assert!(
        optimal_error <= uniform_error * 1.0001,
        "optimal error {optimal_error} should not exceed uniform error {uniform_error}"
    );

    // Knots must remain in non-decreasing order.
    let knots = spline.knots();
    for w in knots.windows(2) {
        assert!(w[0] <= w[1], "knots must stay sorted: {:?}", knots);
    }
}

#[test]
fn fitted_curve_is_monotonically_decreasing_in_trend() {
    // The reference data is decreasing; the smoothed fit should broadly follow.
    let points = reference_dataset();
    let mut spline = build_spline();
    let mut fitter = CurveFitter::new(&spline);

    CurveFitter::initiate_grid(&mut spline, &points);
    CurveFitter::approximate(&mut spline, &points, 1e-9);
    fitter.approximate_with_optimal_grid(&mut spline, &points, 1e-9, 1e-3, 1e-3);

    let left = spline.value(spline.left_bound());
    let right = spline.value(spline.right_bound());
    assert!(left > right, "fit should trend downward: {left} -> {right}");
}
