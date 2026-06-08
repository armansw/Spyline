# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Complete Rust rewrite** of the original Java implementation:
  - `point` / `points` value types.
  - `spline` — B-spline evaluation, derivatives and knot manipulation.
  - `curve_fitter` — weighted least-squares fitting with optimal knot placement.
- Self-contained Cholesky solver replacing the Apache Commons Math dependency,
  giving the crate zero runtime dependencies.
- Command-line driver that reproduces `result.txt` from the bundled dataset.
- Unit, integration and doc tests.
- CI matrix (Linux / macOS / Windows, stable + beta) with formatting and clippy gates.
- Release pipeline producing cross-platform binaries on tagged versions.
- GitHub Pages documentation deployment.
- Multi-stage Dockerfile and GHCR container publishing.

### Changed

- Output number formatting mirrors Java's `Double.toString` so generated
  `result.txt` matches the reference to floating-point round-off.

### Notes

- The optimal-grid search preserves the original control-flow exactly,
  including the quirk where the loop counter is never incremented.
