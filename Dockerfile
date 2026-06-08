# ---- build stage ----
FROM rust:1-slim AS builder
WORKDIR /usr/src/spyline

# Cache dependencies first.
COPY Cargo.toml ./
RUN mkdir src \
    && echo "fn main() {}" > src/main.rs \
    && echo "" > src/lib.rs \
    && cargo build --release || true
RUN rm -rf src

# Build the real sources.
COPY . .
RUN cargo build --release

# ---- runtime stage ----
FROM debian:bookworm-slim AS runtime
LABEL org.opencontainers.image.title="spyline" \
      org.opencontainers.image.description="Weighted least-squares B-spline curve fitting" \
      org.opencontainers.image.source="https://github.com/armansw/spyline"
WORKDIR /app
COPY --from=builder /usr/src/spyline/target/release/spyline /usr/local/bin/spyline
ENTRYPOINT ["spyline"]
