# syntax=docker/dockerfile:1

FROM rust:1.79-bullseye AS build
WORKDIR /app

# cache dependencies
COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo 'fn main() {}' > src/main.rs && cargo build --release && rm -rf src

# copy actual sources
COPY src ./src
COPY config.example.toml ./
COPY book ./book

RUN cargo build --release --bin AnnoQC

FROM debian:bookworm-slim AS runtime

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends diamond mafft ca-certificates \ 
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /app/target/release/AnnoQC /usr/local/bin/annoqc

ENV DIAMOND_BIN=/usr/bin/diamond \
    MAFFT_BIN=/usr/bin/mafft

VOLUME ["/data"]
WORKDIR /data

ENTRYPOINT ["annoqc"]
CMD ["--help"]
