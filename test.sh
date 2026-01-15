#!/usr/bin/env bash
set -euo pipefail

# Test launcher: prepare (databases) + 200/1000 analyses with informative logging.
# - Prepare builds SwissProt DMND and Aves ref-proteome DB (idempotent)
# - Each analyze run logs to <out>/run.log for inspection.
# - Set ALIGNER=spoa to run the SPOA backend (MAFFT flags still passed for fallback).

prepare_refs() {
  echo "[test.sh] Starting prepare at $(date -Is)" | tee -a prepare_run.log
  # Ensure reference data is present (SwissProt, taxdump, Pfam, README)
  pixi run bash scripts/fetch_reference_data.sh share 2>&1 | tee -a prepare_run.log
  # Build SwissProt DMND and run checkpointed clustering + Aves refprot makedb (if README present)
  RUST_LOG=info pixi run cargo run -- \
    prepare \
    --fasta share/uniprot/uniprot_sprot.fasta.gz \
    --db-out uniprot_sprot.dmnd \
    --threads 16 \
    --resume \
    --log-format json 2>&1 | tee -a prepare_run.log
  echo "[test.sh] Finished prepare at $(date -Is)" | tee -a prepare_run.log
}

run_analysis() {
  local FASTA="$1"
  local OUT="$2"
  local THREADS="$3"
  shift 3 || true
  local EXTRA_ARGS=("$@")
  local ALIGNER=${ALIGNER:-mafft}
  mkdir -p "$OUT"
  local LOG="$OUT/run.log"
  echo "[test.sh] Starting $(basename "$FASTA") run into $OUT at $(date -Is)" | tee -a "$LOG"
  local START_TS=$(date +%s)

  local REFPROT_DB="share/refprot/aves/aves_refprot.dmnd"
  local REFPROT_FLAGS=()
  if [[ -f "$REFPROT_DB" ]]; then
    REFPROT_FLAGS=(--refprot-db "$REFPROT_DB")
  fi

  RUST_LOG=info pixi run cargo run -- \
    --config config.weights_demo.toml \
    analyze \
    --fasta "$FASTA" \
    --db uniprot_sprot.dmnd \
    --out "$OUT" \
    --threads "$THREADS" \
    --log-format json \
    --diamond-mode single \
    --aligner "$ALIGNER" \
    --mafft-bin mafft \
    --pfam-db share/pfam/Pfam-A.hmm \
    --pfam-clans share/pfam/Pfam-A.clans.tsv \
    --mafft-fast \
    --mafft-threads-per-job 2 \
    --mafft-max-jobs 6 \
    --reference-fasta share/uniprot/uniprot_sprot.fasta.gz \
    --hmmscan-bin hmmscan \
    --enable-taxonomy \
    --taxonomy-min-consensus 3 \
    --taxonomy-top-hits 20 \
    --csv-verbose \
    --export-high \
    --export-high-path "$OUT/high_scoring.faa" \
    "${REFPROT_FLAGS[@]}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee -a "$LOG"

  local END_TS=$(date +%s)
  local ELAPSED=$((END_TS - START_TS))
  printf '[test.sh] Finished %s run at %s (elapsed %ds)\n' "$(basename "$FASTA")" "$(date -Is)" "$ELAPSED" | tee -a "$LOG"
}

run_200() {
  run_analysis a9_head200.faa results_a9_head200_refprot 12
}

run_1000() {
  run_analysis a9_head1000.faa results_a9_head1000_refprot 12
}

run_5000() {
  run_analysis a9_head5000.faa results_a9_head5000_refprot 16
}

run_full() {
  run_analysis a9.faa results_a9_full_refprot 24
}

case "${1:-}" in
  prepare)
    prepare_refs
    ;;
  200)
    run_200
    ;;
  1000)
    run_1000
    ;;
  5000)
    run_5000
    ;;
  full)
    run_full
    ;;
  both)
    run_200
    run_1000
    ;;
  all|"")
    prepare_refs
    run_200
    run_1000
    run_5000
    run_full
    ;;
  *)
    echo "Usage: $0 [prepare|200|1000|5000|full|both|all]" >&2
    exit 1
    ;;
esac
