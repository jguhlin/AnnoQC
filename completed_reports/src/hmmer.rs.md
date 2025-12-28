# Checklist: src/hmmer.rs

## Severity counts
- HIGH: 6
- MEDIUM: 6

## Issues

- [ ] [HIGH] Mutex contention in task queue and results; use channels or concurrent data structures (performance) (function: run_hmmscan_batch_opts) lines [400, 415] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Missing error handling for stdin write failures (code_safety) (function: run_hmmscan) lines [210, 220] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Assumes clan mappings are exhaustive; no fallback for unmapped accessions (statistical_integrity) (function: collapse_by_clan) lines [600, 610] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Allocates Vec per line causing O(n) heap allocs; use slice iteration (performance) (function: parse_domtblout) lines [138, 152] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc comment for parameters and return value (docstrings) (function: analyze_orphan_domains) lines [71, 92] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Arc<Mutex<>> serialization bottleneck under high concurrency (performance) (function: run_hmmscan_batch_opts) lines [320, 350] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] unwrap_or_default() masks potential errors in thread closure (code_safety) (function: run_hmmscan_batch_opts) lines [340, 350] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Hardcoded core/accessory thresholds lack validation (statistical_integrity) (function: domains_architecture_score) lines [520, 550] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] unwrap_or on parsing can silently drop hits (code_safety) (function: parse_domtblout) lines [240, 260] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Score calculation may return 1.0 when status is None, potentially masking edge cases. (statistical_integrity) (function: analyze_orphan_domains) lines [70, 80] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Contention on Arc<Mutex> results can bottleneck with many threads; suggest per-thread aggregation. (performance) (function: run_hmmscan_batch_opts) lines [260, 300] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Missing doc comment and parameter/return descriptions; unclear contract. (docstrings) (function: domains_architecture_score) lines [145, 155] | model: nemotron-3-nano:30b | block: module