//! # Spyline
//!
//! A weighted least-squares **B-spline curve fitting** library, ported from the
//! original Java implementation to safe, dependency-free Rust.
//!
//! This module currently exposes the core data types; the spline model and
//! fitter are introduced in later modules.

#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::same_item_push,
    clippy::int_plus_one,
    clippy::too_many_arguments
)]

pub mod point;
pub mod points;

pub use point::Point;
pub use points::Points;
