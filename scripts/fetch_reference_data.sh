#!/usr/bin/env bash
# Fetches reference datasets (SwissProt, NCBI taxonomy, Pfam) for offline analysis.
# Usage: scripts/fetch_reference_data.sh [output_dir]
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR=${1:-"${ROOT_DIR}/share"}
UNIPROT_DIR="${OUT_DIR}/uniprot"
TAXONOMY_DIR="${OUT_DIR}/taxonomy"
PFAM_DIR="${OUT_DIR}/pfam"
REFPROT_DIR="${OUT_DIR}/uniprot/reference_proteomes"

mkdir -p "${UNIPROT_DIR}" "${TAXONOMY_DIR}" "${PFAM_DIR}" "${REFPROT_DIR}"

curl_fetch() {
  local url="$1"
  local dest="$2"
  if [[ -f "${dest}" ]]; then
    echo "[fetch] Skipping existing ${dest##*/}"
    return 0
  fi
  echo "[fetch] Downloading ${url}"
  curl -L --fail --retry 4 --retry-delay 5 -o "${dest}.partial" "${url}"
  mv "${dest}.partial" "${dest}"
}

# --- SwissProt FASTA ---
SP_FASTA="${UNIPROT_DIR}/uniprot_sprot.fasta.gz"
curl_fetch "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" "${SP_FASTA}"

# --- NCBI Taxonomy dump ---
TAX_ARCHIVE="${TAXONOMY_DIR}/new_taxdump.tar.gz"
curl_fetch "https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" "${TAX_ARCHIVE}"
if [[ ! -d "${TAXONOMY_DIR}/new_taxdump" ]]; then
  echo "[extract] Unpacking NCBI taxdump"
  mkdir -p "${TAXONOMY_DIR}/new_taxdump"
  tar -C "${TAXONOMY_DIR}/new_taxdump" -xf "${TAX_ARCHIVE}"
fi

# --- Pfam models ---
PFAM_BASE="http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_HMM_GZ="${PFAM_DIR}/Pfam-A.hmm.gz"
PFAM_HMM="${PFAM_DIR}/Pfam-A.hmm"
PFAM_DAT_GZ="${PFAM_DIR}/Pfam-A.hmm.dat.gz"
PFAM_CLANS_GZ="${PFAM_DIR}/Pfam-A.clans.tsv.gz"
PFAM_SEED_GZ="${PFAM_DIR}/Pfam-A.seed.gz"

curl_fetch "${PFAM_BASE}/Pfam-A.hmm.gz" "${PFAM_HMM_GZ}"
curl_fetch "${PFAM_BASE}/Pfam-A.hmm.dat.gz" "${PFAM_DAT_GZ}"
curl_fetch "${PFAM_BASE}/Pfam-A.clans.tsv.gz" "${PFAM_CLANS_GZ}"
curl_fetch "${PFAM_BASE}/Pfam-A.seed.gz" "${PFAM_SEED_GZ}"

if [[ ! -f "${PFAM_HMM}" ]]; then
  echo "[extract] Decompressing Pfam-A.hmm"
  gunzip -c "${PFAM_HMM_GZ}" > "${PFAM_HMM}.partial"
  mv "${PFAM_HMM}.partial" "${PFAM_HMM}"
fi

for gz in "${PFAM_DAT_GZ}" "${PFAM_CLANS_GZ}"; do
  dest="${gz%.gz}"
  if [[ ! -f "${dest}" ]]; then
    echo "[extract] Decompressing ${gz##*/}"
    gunzip -c "${gz}" > "${dest}.partial"
    mv "${dest}.partial" "${dest}"
  fi
done

if command -v hmmpress >/dev/null 2>&1; then
  if [[ ! -f "${PFAM_HMM}.h3f" ]]; then
    echo "[hmmpress] Generating HMMER indices"
    (cd "${PFAM_DIR}" && hmmpress "${PFAM_HMM}")
  else
    echo "[hmmpress] Pfam indices already present"
  fi
else
  echo "[warn] hmmpress not found in PATH; skipping index generation" >&2
fi

echo "[done] Reference data ready under ${OUT_DIR}"

# --- UniProt Reference Proteomes README (mapping proteome IDs → taxa) ---
# Downloaded last so previous steps aren't affected. Skips if present.
REFPROT_README_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"
REFPROT_README_PATH="${REFPROT_DIR}/README"
curl_fetch "${REFPROT_README_URL}" "${REFPROT_README_PATH}"
