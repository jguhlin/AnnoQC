#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ./scripts/bench_aligners.sh --fasta <path> --db <path> --out-root <path> [options]

Options:
  --threads <n>       Thread count (default: 8)
  --config <path>     Optional config.toml
  --mafft-fast        Use MAFFT fast mode (default: enabled)
  --no-mafft-fast     Disable MAFFT fast mode
USAGE
}

fasta=""
db=""
out_root=""
threads=8
config=""
mafft_fast=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta)
      fasta="$2"; shift 2;;
    --db)
      db="$2"; shift 2;;
    --out-root)
      out_root="$2"; shift 2;;
    --threads)
      threads="$2"; shift 2;;
    --config)
      config="$2"; shift 2;;
    --mafft-fast)
      mafft_fast=1; shift;;
    --no-mafft-fast)
      mafft_fast=0; shift;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown arg: $1" >&2
      usage; exit 1;;
  esac
 done

if [[ -z "$fasta" || -z "$db" || -z "$out_root" ]]; then
  usage
  exit 1
fi

run_bench() {
  local aligner="$1"
  local out_dir="$out_root/$aligner"
  if [[ -e "$out_dir" ]]; then
    echo "Output already exists: $out_dir" >&2
    exit 1
  fi
  local args=(pixi run cargo run -- analyze --fasta "$fasta" --db "$db" --out "$out_dir" --threads "$threads" --aligner "$aligner")
  if [[ -n "$config" ]]; then
    args=(pixi run cargo run -- --config "$config" analyze --fasta "$fasta" --db "$db" --out "$out_dir" --threads "$threads" --aligner "$aligner")
  fi
  if [[ "$mafft_fast" -eq 1 ]]; then
    args+=(--mafft-fast)
  fi
  "${args[@]}"
  python - <<PY "$out_dir/run_metrics.json"
import json,sys
path = sys.argv[1]
with open(path) as f:
    data = json.load(f)
steps = {s["name"]: float(s["seconds"]) for s in data.get("steps", [])}
msa = steps.get("mafft", 0.0)
total = sum(steps.values())
print(f"{path}: total={total:.3f}s mafft_step={msa:.3f}s")
PY
}

mkdir -p "$out_root"
run_bench mafft
run_bench spoa
