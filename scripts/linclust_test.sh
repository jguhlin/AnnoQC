#!/usr/bin/env bash
# Reproduce/diagnose DIAMOND linclust crashes with controlled runs.
# Usage:
#   scripts/linclust_test.sh [FASTA] [OUTDIR]
# Defaults:
#   FASTA = share/uniprot/uniprot_sprot.fasta.gz
#   OUTDIR = linclust_diagnose_$(date +%Y%m%d_%H%M%S)

set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FASTA_IN=${1:-"${ROOT_DIR}/share/uniprot/uniprot_sprot.fasta.gz"}
OUTROOT=${2:-"${ROOT_DIR}/linclust_diagnose_$(date +%Y%m%d_%H%M%S)"}

mkdir -p "${OUTROOT}" "${OUTROOT}/logs" "${OUTROOT}/tmp"

log() { echo "[linclust_test] $*" | tee -a "${OUTROOT}/logs/run.log"; }

cmd_exists() { command -v "$1" >/dev/null 2>&1; }

if ! cmd_exists pixi; then
  echo "pixi not found in PATH. Please install/activate Pixi env first." >&2
  exit 1
fi

if ! pixi run diamond --version >/dev/null 2>&1; then
  echo "diamond not available via pixi. Try: pixi run diamond --version" >&2
  exit 1
fi

log "diamond version: $(pixi run diamond --version | head -n1)"
log "input FASTA: ${FASTA_IN}"

# Prepare an uncompressed FASTA copy (some environments behave better without gzip streaming)
FASTA_UN="${OUTROOT}/uniprot_sprot.uncompressed.fasta"
if [[ "${FASTA_IN}" == *.gz ]]; then
  log "Decompressing input to ${FASTA_UN}"
  gunzip -c "${FASTA_IN}" > "${FASTA_UN}"
else
  log "Copying input to ${FASTA_UN}"
  cp -f "${FASTA_IN}" "${FASTA_UN}"
fi

# Helper: take first N records from FASTA (header lines start with '>')
take_fasta_n() {
  local n=$1; shift
  awk -v N="$n" 'BEGIN{c=0} /^>/{c++} c<=N{print}' "$@"
}

# Create a 50k-record subset for rapid iteration
FASTA_50K="${OUTROOT}/subset_50k.fasta"
if [[ ! -f "${FASTA_50K}" ]]; then
  log "Creating subset of first 50,000 records -> ${FASTA_50K}"
  take_fasta_n 50000 "${FASTA_UN}" > "${FASTA_50K}"
fi

run_case() {
  local label=$1; shift
  local fasta=$1; shift
  local threads=$1; shift
  local extra=("$@")
  local outdir="${OUTROOT}/${label}"
  mkdir -p "${outdir}" "${outdir}/tmp"
  log "Running case '${label}' (threads=${threads}) on $(basename "${fasta}")"
  # Ensure per-run diamond.log is unique
  rm -f diamond.log
  # Use a bounded tmpdir and verbose logging
  ( 
    cd "${ROOT_DIR}" && \
    TMPDIR="${outdir}/tmp" pixi run diamond linclust \
      -d "${fasta}" \
      -o "${outdir}/clusters" \
      --approx-id 50 \
      --threads "${threads}" \
      --verbose --log "${extra[@]}"
  )
  local status=$?
  if [[ -f diamond.log ]]; then
    mv -f diamond.log "${outdir}/diamond.log"
  fi
  echo "${status}" > "${outdir}/exit.code"
  if [[ ${status} -ne 0 ]]; then
    log "Case '${label}' FAILED (exit ${status}). See ${outdir}/diamond.log"
  else
    log "Case '${label}' OK"
  fi
}

# Cases to probe stability/resource sensitivity
# Note: DIAMOND linclust does not accept some blastp flags (e.g., --block-size, --index-chunks).
# Keep the test minimal and portable: only flags known to be supported by 'linclust'.
run_case full_t8        "${FASTA_UN}" 8
run_case full_t4        "${FASTA_UN}" 4
run_case subset50k_t8   "${FASTA_50K}" 8
run_case subset50k_t4   "${FASTA_50K}" 4

log "Done. Summaries:"
for d in full_t8 full_t4 subset50k_t8 subset50k_t4; do
  printf "  %-14s exit=%s\n" "$d" "$(cat "${OUTROOT}/${d}/exit.code" 2>/dev/null || echo NA)" | tee -a "${OUTROOT}/logs/run.log"
done

log "Inspect logs under: ${OUTROOT}/<case>/diamond.log"
log "You can share the failing case dir for debugging."
