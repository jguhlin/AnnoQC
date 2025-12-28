# Checklist: src/taxonomy.rs

## Severity counts
- HIGH: 4
- LOW: 1

## Issues

- [ ] [HIGH] Unbounded lineage reconstruction may cause stack overflow or OOM (performance) (function: reconstruct_lineage) lines [230, 250] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Unvalidated string parsing may panic on malformed input (code_safety) (function: parse_taxid) lines [300, 305] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Non-deterministic taxon selection due to HashMap iteration order; results vary across runs. (statistical_integrity) (function: summarize_panel) lines [260, 285] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Repeated HashMap inserts per lineage element cause O(N*L) overhead; can be streamlined. (performance) (function: fallback_coarse_consensus) lines [210, 235] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Struct fields lack explanatory doc comments; unclear defaults. (docstrings) (function: TaxonomyConsensusConfig) lines [30, 70] | model: nemotron-3-nano:30b | block: module