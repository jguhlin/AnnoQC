#!/usr/bin/env bash
set -euo pipefail

# Compare MAFFT vs SPOA on head200 with clean outputs and timing.

MAFFT_OUT="results_a9_head200_mafft_bench"
SPOA_OUT="results_a9_head200_spoa_bench"
FASTA="a9_head200.faa"
DB="uniprot_sprot.dmnd"
REF="share/uniprot/uniprot_sprot.fasta.gz"
PFAM_DB="share/pfam/Pfam-A.hmm"
PFAM_CLANS="share/pfam/Pfam-A.clans.tsv"
REFPROT_DB="share/refprot/aves/aves_refprot.dmnd"
THREADS=12

run_one() {
  local ALIGNER=$1
  local OUT=$2
  rm -rf "$OUT"
  mkdir -p "$OUT"
  echo "[$(date -Is)] running $ALIGNER -> $OUT"
  CARGO_HOME=.cargo /usr/bin/time -f "wall %e\nuser %U\nsys %S\nmaxrss %M" \
    pixi run cargo run -- \
      --config config.weights_demo.toml \
      analyze \
      --fasta "$FASTA" \
      --db "$DB" \
      --out "$OUT" \
      --threads "$THREADS" \
      --log-format json \
      --diamond-mode single \
      --aligner "$ALIGNER" \
      --mafft-bin mafft \
      --mafft-fast \
      --mafft-threads-per-job 2 \
      --mafft-max-jobs 6 \
      --reference-fasta "$REF" \
      --hmmscan-bin hmmscan \
      --pfam-db "$PFAM_DB" \
      --pfam-clans "$PFAM_CLANS" \
      --enable-taxonomy \
      --taxonomy-min-consensus 3 \
      --taxonomy-top-hits 20 \
      --csv-verbose \
      --export-high \
      --export-high-path "$OUT/high_scoring.faa" \
      --refprot-db "$REFPROT_DB" \
      >"$OUT/run.log" 2>&1
  echo "[$(date -Is)] finished $ALIGNER"
}

run_one mafft "$MAFFT_OUT"
run_one spoa "$SPOA_OUT"
