# Spyline

[![CI](https://github.com/armansw/spyline/actions/workflows/ci.yml/badge.svg)](https://github.com/armansw/spyline/actions/workflows/ci.yml)
[![Release](https://github.com/armansw/spyline/actions/workflows/release.yml/badge.svg)](https://github.com/armansw/spyline/actions/workflows/release.yml)
[![Docs](https://github.com/armansw/spyline/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/armansw/spyline/actions/workflows/deploy-docs.yml)

Spyline performs **weighted least-squares B-spline curve fitting** with
automatic, optimal knot placement. Given a set of weighted `(x, y)`
observations it builds a smooth polynomial spline approximation of the
underlying function.

This is a **Rust** rewrite of the original Java implementation. The numerical
results match the Java reference to within floating-point round-off
(`< 1e-11` relative error), and the Apache Commons Math `CholeskyDecomposition`
dependency has been replaced by a small, self-contained Cholesky solver — so
the crate has **zero runtime dependencies**.

## Features

- Cubic (and arbitrary-degree) clamped B-splines with De Boor evaluation.
- Spline value and knot derivatives.
- Weighted least-squares fitting via the normal equations + Cholesky solve.
- Optional second-derivative smoothing term.
- Optimal interior-knot placement through a penalised line search.
- No external crates required.

## Project layout

| Module                | Responsibility                                   |
| --------------------- | ------------------------------------------------ |
| `point`               | A single weighted observation `(x, y, w)`.       |
| `points`              | A column-oriented collection of observations.    |
| `spline`              | B-spline evaluation, derivatives, knot edits.    |
| `curve_fitter`        | Least-squares fitting + optimal knot placement.  |

## Build & run

```sh
# Build the library and CLI driver.
cargo build --release

# Run the driver: fits the bundled dataset and writes result.txt.
cargo run --release

# Run the full test suite (unit + integration + doctests).
cargo test
```

The CLI driver reproduces the original program: it fits the bundled sample
dataset on both a uniform interior grid and an optimised grid, then writes the
sampled curves and knot locations to `result.txt`.

## Library usage

```rust
use spyline::{Points, Spline, CurveFitter};

let x = vec![1.0, 3.0, 7.0, 15.0, 31.0, 63.0, 127.0, 255.0];
let y = vec![900.0, 779.0, 692.0, 617.0, 551.0, 476.0, 387.0, 283.0];
let w = vec![1.0; x.len()];
let points = Points::new(x.clone(), y, w);

let knots = vec![x[0], 4.0, 8.0, 12.0, 16.0, *x.last().unwrap(), 17.0];
let coefficients = vec![0.0; knots.len() - 2];
let mut spline = Spline::new(coefficients, knots, 3);

let mut fitter = CurveFitter::new(&spline);
CurveFitter::initiate_grid(&mut spline, &points);
fitter.approximate_with_optimal_grid(&mut spline, &points, 1e-9, 1e-3, 1e-3);

println!("fit at x=10: {}", spline.value(10.0));
```

## Docker

```sh
docker build -t spyline .
docker run --rm spyline
```

Images are also published to the GitHub Container Registry on every push to
the default branch and on tagged releases (`ghcr.io/armansw/spyline`).

## Continuous integration & deployment

- **CI** (`ci.yml`) — formatting, clippy, and a cross-platform test matrix
  (Linux / macOS / Windows, stable + beta).
- **Release** (`release.yml`) — builds platform binaries and publishes a
  GitHub Release on every `v*.*.*` tag.
- **Docs** (`deploy-docs.yml`) — builds `rustdoc` and deploys it to GitHub
  Pages.
- **Docker** (`docker.yml`) — builds and pushes the container image to GHCR.

## License

Licensed under the MIT License — see [LICENSE](LICENSE).
