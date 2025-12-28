# Checklist: src/diamond.rs

## Severity counts
- HIGH: 8
- MEDIUM: 3
- INFO: 2
- LOW: 1
- VERIFY: 1

## Issues

- [ ] [HIGH] Repeated fs::read causes redundant I/O; use BufReader or in-memory buffering (performance) (function: run_chunk) lines [265] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Use of unwrap_or on parsed values may silently drop invalid data (code_safety) (function: parse_tsv_stats) lines [230] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Division by zero or incorrect coverage due to None qlen_map lookup (statistical_integrity) (function: parse_tsv_stats) lines [290] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Misalignment of TSV column parsing leads to incorrect hit statistics. (statistical) (function: parse_tsv_stats) lines [260, 315] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Temporary TSV file is read fully into memory before appending, causing extra I/O and memory usage. (performance) (function: run_chunk) lines [100, 130] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Parsing uses unwrap/potential panics on malformed input; no explicit error handling. (safety) (function: run_chunk) lines [260, 315] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Doc comment omits mention of retry behavior and blocking time. (documentation) (function: blastp_once) lines [30, 55] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] The function `blastp_ungapped` has a nested loop structure that is O(n^2) in the worst case due to the `while` loop inside the `for` loop. This can be optimized by vectorizing the operations using NumPy or by restructuring the logic to reduce the complexity. (performance) (function: blastp_ungapped) lines [1, 100] | model: rnj-1:8b | block: module
- [ ] [HIGH] Missing query length normalization in qcov computation (statistical_integrity) (function: parse_tsv_stats) lines [270, 275] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Invalid taxid tokens silently dropped without logging (code_safety) (function: parse_taxid) lines [360, 365] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Creates temporary files per chunk and reads them back, causing O(N) I/O; replace with direct streaming write. (performance) (function: run_chunk) lines [70, 110] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Rebuilds Command argument list on each retry, allocating new Strings repeatedly; cache arguments to reduce allocations. (performance) (function: blastp_once) lines [30, 55] | model: nemotron-3-nano:30b | block: module
- [ ] [VERIFY] Indexes into `cols` with hard‑coded positions up to 15 without length check; could panic on malformed TSV. (code_safety) (function: parse_tsv_grouped) lines [200, 230] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Missing doc comment explaining retry logic and skip‑if‑output‑exists behavior. (docstrings) (function: blastp_once) lines [15, 30] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] `coverage_ratio` division lacks finiteness check; may propagate NaN when `scov` is zero or invalid. (statistical_integrity) (function: parse_tsv_stats) lines [140, 155] | model: nemotron-3-nano:30b | block: module