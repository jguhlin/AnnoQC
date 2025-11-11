#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v pixi >/dev/null 2>&1; then
  echo "pixi is required to build release binaries; install it or run via 'pixi run'." >&2
  exit 1
fi

pixi run cargo build --release

VERSION=$(grep -m1 '^version' Cargo.toml | sed -E 's/version\s*=\s*"([^"]+)"/\1/')
if [[ -z "$VERSION" ]]; then
  echo "Unable to determine version from Cargo.toml" >&2
  exit 1
fi

TARGET_TRIPLE=$(rustc -vV | awk '/host:/ {print $2}')
BIN_PATH="target/release/AnnoQC"
if [[ ! -x "$BIN_PATH" ]]; then
  BIN_PATH="target/release/annoqc"
fi
if [[ ! -x "$BIN_PATH" ]]; then
  echo "Release binary not found at $BIN_PATH" >&2
  exit 1
fi

DIST_DIR="$ROOT_DIR/dist"
PKG_NAME="annoqc-${VERSION}-${TARGET_TRIPLE}"
STAGE_DIR="$DIST_DIR/$PKG_NAME"
rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR"

cp "$BIN_PATH" "$STAGE_DIR/annoqc"
chmod 755 "$STAGE_DIR/annoqc"

for asset in LICENSE config.example.toml config.weights_demo.toml; do
  if [[ -f "$asset" ]]; then
    cp "$asset" "$STAGE_DIR/"
  fi
done

echo "Built artifacts staged at $STAGE_DIR"
mkdir -p "$DIST_DIR"
tar -C "$DIST_DIR" -czf "$DIST_DIR/${PKG_NAME}.tar.gz" "$PKG_NAME"
sha256sum "$DIST_DIR/${PKG_NAME}.tar.gz" > "$DIST_DIR/${PKG_NAME}.tar.gz.sha256"
rm -rf "$STAGE_DIR"

echo "Release package written to $DIST_DIR/${PKG_NAME}.tar.gz"
