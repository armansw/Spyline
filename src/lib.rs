//! # Spyline
//!
//! A weighted least-squares **B-spline curve fitting** library, ported from the
//! original Java implementation to safe, dependency-free Rust.
//!
//! The crate is organised into four small modules that mirror the original
//! class structure:
//!
//! * [`point`] — a single weighted observation.
//! * [`points`] — a column-oriented collection of observations.
//! * [`spline`] — B-spline evaluation, derivatives and knot manipulation.
//! * [`curve_fitter`] — least-squares fitting with optimal knot placement.
//!
//! ## Example
//!
//! ```
//! use spyline::points::Points;
//! use spyline::spline::Spline;
//! use spyline::curve_fitter::CurveFitter;
//!
//! let x = vec![1.0, 3.0, 7.0, 15.0, 31.0, 63.0, 127.0, 255.0];
//! let y = vec![900.0, 779.0, 692.0, 617.0, 551.0, 476.0, 387.0, 283.0];
//! let w = vec![1.0; x.len()];
//! let points = Points::new(x.clone(), y, w);
//!
//! let knots = vec![x[0], 4.0, 8.0, 12.0, 16.0, *x.last().unwrap(), 17.0];
//! let coefficients = vec![0.0; knots.len() - 2];
//! let mut spline = Spline::new(coefficients, knots, 3);
//!
//! CurveFitter::initiate_grid(&mut spline, &points);
//! assert!(CurveFitter::approximate(&mut spline, &points, 1e-9));
//! let _y = spline.value(10.0);
//! ```

// The numerical routines are a faithful port that relies on explicit index
// arithmetic; the lints below would obscure the correspondence with the source.
#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::same_item_push,
    clippy::int_plus_one,
    clippy::too_many_arguments
)]

pub mod curve_fitter;
pub mod point;
pub mod points;
pub mod spline;

pub use curve_fitter::CurveFitter;
pub use point::Point;
pub use points::Points;
pub use spline::Spline;
