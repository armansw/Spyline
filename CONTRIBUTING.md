# Contributing to Spyline

Thanks for your interest in improving Spyline! This document describes the
workflow and the checks your changes must pass.

## Development workflow

1. Create a feature branch off the default branch:
   ```sh
   git switch -c feature/my-change
   ```
2. Make your changes with accompanying tests.
3. Run the local quality gates (see below).
4. Open a pull request describing the motivation and the approach.

## Quality gates

All of the following must pass locally and in CI before a change can be merged:

```sh
cargo fmt --all -- --check      # formatting
cargo clippy --all-targets -- -D warnings   # lints
cargo test --all-features       # unit + integration + doc tests
```

## Coding guidelines

- The numerical modules (`spline`, `curve_fitter`) are a faithful port of the
  original Java algorithms. When changing them, preserve the numerical
  behaviour and update the integration tests accordingly.
- Public items must carry doc comments; examples in doc comments are run as
  tests.
- Keep the crate dependency-free unless there is a compelling reason otherwise.

## Commit messages

Use conventional, imperative commit subjects, e.g.:

```
feat(spline): add periodic knot support
fix(curve_fitter): guard against zero-length intervals
docs(readme): document the Docker workflow
```
