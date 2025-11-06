# Analyze

AnnoQC’s `analyze` subcommand combines DIAMOND homology searches, MAFFT conserved-region alignments, and intrinsic sequence checks to produce a per-gene scorecard. The JSONL output exposes detailed evidence, while the CSV mirrors the headline metrics for downstream dashboards.

## Prerequisites

- DIAMOND installed and the SwissProt database prepared via `cargo run -- prepare`.
- (Optional) MAFFT installed when you want conserved-region analysis (`--mafft-bin`).
- A reference FASTA that contains the DIAMOND subjects (defaults to `uniprot_sprot.fasta.gz`).

## Example Run

1. Copy and adjust `config.example.toml`:

```toml
fasta = "examples/demo.faa"
db = "uniprot_sprot.dmnd"
out = "results"
top = 5
diamond_bin = "diamond"
mafft_bin = "/usr/bin/mafft"
reference_fasta = "uniprot_sprot.fasta.gz"
[scoring.weights]
homology = 0.6
intrinsic = 0.4
```

2. Run the analyzer:

```bash
cargo run -- analyze --config config.example.toml
```

This produces `results/qc_report.jsonl`, `results/qc_summary.csv`, and `results/run.json` (manifest with tool versions and resolved configuration).

## JSONL Scorecard Snapshot

Each JSON line records the final score, per-pillar components, evidence blocks, and warnings. New fields include coverage metrics and a fusion/split flag, plus intrinsic metrics. When MAFFT is enabled, an alignment block is present.

```json
{
  "gene_id": "Gene00042",
  "final_score": 0.87,
  "classification": "High",
  "score_components": { "homology": 0.92, "intrinsic": 0.81 },
  "homology": {
    "hits_count": 12,
    "top_hit": "sp|P12345|REF_HUMAN",
    "top_bitscore": 212.6,
    "top_evalue": 3.2e-68,
    "top_qcov": 0.94,
    "top_scov": 0.91,
    "bitscore_density": 3.21,
    "coverage_ratio": 1.03,
    "coverage_delta": 0.03,
    "fusion_split_flag": false
  },
  "intrinsic": {
    "ambiguous_fraction": 0.0,
    "max_homopolymer": 5,
    "low_complexity_fraction": 0.08,
    "low_complexity_windows": 2,
    "orf_has_start": true,
    "orf_has_stop": true,
    "orf_internal_stop": false
  },
  "alignment": {
    "mafft_enabled": true,
    "strategy_used": "Add",
    "conserved_fraction": 0.89,
    "pairwise_identity": 0.86,
    "sequences_aligned": 6,
    "top_hits_considered": 5,
    "query_gap_fraction": 0.04,
    "gap_run_count": 1,
    "max_gap_run": 4,
    "motif_mismatch_fraction": 0.11
  },
  "warnings": []
}
```

Rows with potential issues advertise warnings such as missing start/stop codons, internal stops, excessive low-complexity content, or suspicious fusion/split coverage.

## CSV Columns

`qc_summary.csv` contains headline homology numbers, coverage metrics, and a `warnings` column for quick filtering:

```
gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,taxonomy_score,taxonomy_status,warnings
```

The CSV is suitable for spreadsheets or dashboards, while the JSONL is richer for downstream pipelines.
