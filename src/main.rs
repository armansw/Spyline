//! Command-line driver, ported from the Java `spyline.App`.
//!
//! Fits a cubic smoothing spline to a fixed sample dataset twice — once on a
//! uniform interior grid and once with the optimal-grid optimisation — then
//! writes a sampled curve and the knot locations to `result.txt`.

use std::fs::File;
use std::io::{self, Write};

use spyline::curve_fitter::CurveFitter;
use spyline::points::Points;
use spyline::spline::Spline;

/// Format a float the way Java's `Double.toString` does: integral values keep a
/// trailing `.0` so the output matches the reference reproduction.
fn fmt(v: f64) -> String {
    if v.is_nan() {
        return "NaN".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let s = format!("{v}");
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{s}.0")
    }
}

fn write_series(out: &mut impl Write, label: &str, values: &[f64]) -> io::Result<()> {
    write!(out, "{label}")?;
    for v in values {
        write!(out, "{} ", fmt(*v))?;
    }
    // Java wrote three platform line separators between sections.
    write!(out, "\n\n\n")?;
    Ok(())
}

fn main() -> io::Result<()> {
    let x = vec![
        1.0, 3.0, 7.0, 15.0, 31.0, 63.0, 127.0, 255.0, 383.0, 639.0, 895.0, 1151.0, 1663.0, 2175.0,
        3199.0, 4223.0, 6271.0, 10367.0, 14463.0, 22655.0, 39039.0,
    ];
    let y = vec![
        900.0, 779.0, 692.0, 617.0, 551.0, 476.0, 387.0, 283.0, 255.0, 221.0, 200.0, 185.0, 164.0,
        150.0, 81.0, 65.0, 46.0, 26.0, 16.0, 7.0, 1.0,
    ];
    let w = vec![1.0; x.len()];
    let points = Points::new(x.clone(), y, w);

    // Sample 1000 abscissae uniformly across the data range.
    let n_curve = 1000;
    let first = points.item(0).x();
    let last = points.item(x.len() - 1).x();
    let d = (last - first) / (n_curve - 1) as f64;
    let x_curve: Vec<f64> = (0..n_curve).map(|i| first + d * i as f64).collect();

    // Initial knot grid (boundary + interior knots), matching the Java driver:
    // an extra knot is appended next to the right boundary before sorting.
    let mut knots = vec![first, 4.0, 8.0, 12.0, 16.0, last];
    let coefficients = vec![0.0; knots.len() - 2];
    let penultimate = knots[knots.len() - 2];
    knots.push(penultimate + 1.0);
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut spline = Spline::new(coefficients, knots, 3);
    let mut fitter = CurveFitter::new(&spline);

    let q = 1e-9;
    CurveFitter::initiate_grid(&mut spline, &points);
    CurveFitter::approximate(&mut spline, &points, q);

    let y_curve: Vec<f64> = x_curve.iter().map(|&p| spline.value(p)).collect();

    let uniform_knots = spline.knots();
    let knots_y_uni: Vec<f64> = uniform_knots.iter().map(|&p| spline.value(p)).collect();

    fitter.approximate_with_optimal_grid(&mut spline, &points, q, 1e-3, 1e-3);

    let y_curve_opt: Vec<f64> = x_curve.iter().map(|&p| spline.value(p)).collect();

    let optimal_knots = spline.knots();
    let knots_y: Vec<f64> = optimal_knots.iter().map(|&p| spline.value(p)).collect();

    let mut file = File::create("result.txt")?;
    write_series(&mut file, "x_curve:  ", &x_curve)?;
    write_series(&mut file, "y_curve:  ", &y_curve)?;
    write_series(&mut file, "y_curve_opt:  ", &y_curve_opt)?;
    write_series(&mut file, "knots:  ", &uniform_knots)?;
    write_series(&mut file, "knots_y_uni:  ", &knots_y_uni)?;
    write_series(&mut file, "knots2:  ", &optimal_knots)?;
    write_series(&mut file, "knots_y:  ", &knots_y)?;

    Ok(())
}
