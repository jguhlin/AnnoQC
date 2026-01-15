#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/package_release.sh [--target <triple>] [--docker] [--docker-tag <name>]

Options:
  --target <triple>   Build for the given Rust target triple (defaults to host).
  --docker            Build a Docker image using the local Dockerfile (tagged with the version).
  --docker-tag <tag>  Custom Docker tag (defaults to annqoc:<version>). Implies --docker.
  -h, --help          Show this help.
USAGE
}

TARGET_OVERRIDE=""
DOCKER_BUILD=0
DOCKER_TAG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target)
      TARGET_OVERRIDE="$2"
      shift 2
      ;;
    --docker)
      DOCKER_BUILD=1
      shift
      ;;
    --docker-tag)
      DOCKER_BUILD=1
      DOCKER_TAG="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v pixi >/dev/null 2>&1; then
  echo "pixi is required to build release binaries; install it or run via 'pixi run'." >&2
  exit 1
fi

TARGET_TRIPLE=${TARGET_OVERRIDE:-$(rustc -vV | awk '/host:/ {print $2}')}
CARGO_TARGET_ARGS=()
BIN_SUBDIR="release"
if [[ -n "$TARGET_OVERRIDE" ]]; then
  CARGO_TARGET_ARGS=(--target "$TARGET_TRIPLE")
  BIN_SUBDIR="$TARGET_TRIPLE/release"
fi

pixi run cargo build --release "${CARGO_TARGET_ARGS[@]}"

VERSION=$(grep -m1 '^version' Cargo.toml | sed -E 's/version\s*=\s*"([^"]+)"/\1/')
if [[ -z "$VERSION" ]]; then
  echo "Unable to determine version from Cargo.toml" >&2
  exit 1
fi

BIN_BASE="target/${BIN_SUBDIR}/AnnoQC"
if [[ ! -x "$BIN_BASE" ]]; then
  BIN_BASE="target/${BIN_SUBDIR}/annoqc"
fi
if [[ ! -x "$BIN_BASE" && -x "${BIN_BASE}.exe" ]]; then
  BIN_BASE="${BIN_BASE}.exe"
fi
if [[ ! -x "$BIN_BASE" ]]; then
  echo "Release binary not found at $BIN_BASE" >&2
  exit 1
fi

DIST_DIR="$ROOT_DIR/dist"
PKG_NAME="annoqc-${VERSION}-${TARGET_TRIPLE}"
STAGE_DIR="$DIST_DIR/$PKG_NAME"
rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR"

BIN_NAME="annoqc"
if [[ "$BIN_BASE" == *.exe ]]; then
  BIN_NAME="annoqc.exe"
fi
cp "$BIN_BASE" "$STAGE_DIR/$BIN_NAME"
chmod 755 "$STAGE_DIR/$BIN_NAME"

for asset in LICENSE config.example.toml config.weights_demo.toml; do
  if [[ -f "$asset" ]]; then
    cp "$asset" "$STAGE_DIR/"
  fi
done

echo "Built artifacts staged at $STAGE_DIR"
mkdir -p "$DIST_DIR"
tar -C "$DIST_DIR" -czf "$DIST_DIR/${PKG_NAME}.tar.gz" "$PKG_NAME"
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum "$DIST_DIR/${PKG_NAME}.tar.gz" > "$DIST_DIR/${PKG_NAME}.tar.gz.sha256"
fi
rm -rf "$STAGE_DIR"

echo "Release package written to $DIST_DIR/${PKG_NAME}.tar.gz"

if [[ $DOCKER_BUILD -eq 1 ]]; then
  if ! command -v docker >/dev/null 2>&1; then
    echo "Docker not found; skipping image build" >&2
  else
    TAG=${DOCKER_TAG:-"annoqc:${VERSION}"}
    echo "Building Docker image $TAG"
    docker build --build-arg VERSION="${VERSION}" -t "$TAG" .
    echo "Docker image $TAG built"
  fi
fi
