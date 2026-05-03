# ---- Build stage ----
FROM rust:1-bookworm AS builder

# Optional pre-computed version string (e.g. "0.1.0-deadbeef" or "0.1.0").
# When unset, build.rs falls back to the bare Cargo.toml version because
# .git is not part of the Docker context.
ARG DADA2_RS_VERSION_FULL
ENV DADA2_RS_VERSION_FULL=${DADA2_RS_VERSION_FULL}

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY src/ src/

RUN cargo build --release \
    && strip target/release/dada2-rs

# ---- Runtime stage ----
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/dada2-rs /usr/local/bin/dada2-rs

ENTRYPOINT ["dada2-rs"]
