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

Each JSON line records the final score, per-pillar components, evidence blocks, and warnings:

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

`qc_summary.csv` contains the same headline numbers plus a `warnings` column for quick filtering:

```
gene_id,hits_count,top_hit,top_bitscore,...,mafft_motif_mismatch_fraction,low_complexity_fraction,low_complexity_windows,orf_has_start,orf_has_stop,orf_internal_stop,warnings
Gene00042,12,sp|P12345|REF_HUMAN,212.600,3.20e-68,0.940,0.910,3.210000,1.030,0.030,0,0.0000,5,0.9200,0.8100,0.8700,high,0.890,0.860,6,0.040,1,4,0.1100,0.0800,2,1,1,0,
Gene00077,4,sp|Q8XXX7|ALT_MOUSE,98.300,2.50e-18,0.540,0.310,1.820000,1.742,0.230,1,0.1200,14,0.4200,0.5100,0.4660,medium,0.420,0.380,5,0.120,3,9,0.5800,0.2700,9,0,0,1,Internal stop codon detected;High low-complexity content (27.0% windows)
```

The CSV is suitable for spreadsheets or dashboards, while the JSONL is richer for downstream pipelines.
