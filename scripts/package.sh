#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DIST_DIR="${ROOT_DIR}/dist"
DEFAULT_TARGETS=(
  x86_64-unknown-linux-gnu
  aarch64-apple-darwin
  x86_64-apple-darwin
  x86_64-pc-windows-gnu
)

if [[ $# -gt 0 ]]; then
  TARGETS=("$@")
else
  TARGETS=("${DEFAULT_TARGETS[@]}")
fi

mkdir -p "${DIST_DIR}"

run_builder() {
  local target="$1"
  if command -v cross >/dev/null 2>&1; then
    echo "[package] building ${target} with cross"
    cross build --release --target "${target}"
  else
    echo "[package] building ${target} with cargo"
    cargo build --release --target "${target}"
  fi
}

package_binary() {
  local target="$1"
  local bin_name="AnnoQC"
  local ext=""
  if [[ "${target}" == *"windows"* ]]; then
    ext=".exe"
  fi
  local bin_path="${ROOT_DIR}/target/${target}/release/${bin_name}${ext}"
  if [[ ! -f "${bin_path}" ]]; then
    echo "[package] missing binary for ${target}: ${bin_path}" >&2
    return 1
  fi
  local archive_name="${bin_name}-${target}"
  local archive_path
  if [[ "${target}" == *"windows"* ]]; then
    if command -v zip >/dev/null 2>&1; then
      archive_path="${DIST_DIR}/${archive_name}.zip"
      (cd "$(dirname "${bin_path}")" && zip -q "${archive_path}" "${bin_name}${ext}")
    else
      archive_path="${DIST_DIR}/${archive_name}.tar.gz"
      (cd "$(dirname "${bin_path}")" && tar -czf "${archive_path}" "${bin_name}${ext}")
    fi
  else
    archive_path="${DIST_DIR}/${archive_name}.tar.gz"
    (cd "$(dirname "${bin_path}")" && tar -czf "${archive_path}" "${bin_name}${ext}")
  fi
  echo "[package] packaged ${archive_path}"
}

for target in "${TARGETS[@]}"; do
  run_builder "${target}"
  package_binary "${target}"
done
