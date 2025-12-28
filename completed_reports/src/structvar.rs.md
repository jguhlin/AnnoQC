# Checklist: src/structvar.rs

## Severity counts
- HIGH: 6
- MEDIUM: 2
- LOW: 1

## Issues

- [ ] [HIGH] Inconsistent coverage calculation leads to leakage in fusion/duplication detection (statistical_integrity) (function: analyze) lines [110, 125] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Nested loop in fusion detection causes O(n^2) complexity (performance) (function: analyze) lines [140, 160] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Potential panic on zero-length subject length (code_safety) (function: analyze) lines [105] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Split detection lacks significance testing; may misclassify low-coverage subjects. (statistical) (function: analyze) lines [210, 215] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Fusion pair loop O(k^2) can dominate with many subjects; replace nested loops with sorted gap scan. (performance) (function: analyze) lines [140, 176] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] `unwrap_or(1)` may panic on missing sseqid; replace with expect and error handling. (safety) lines [30, 35] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Potential data leakage in training process - No validation split in training loop (statistical_integrity) (function: train_model) lines [42, 67] | model: rnj-1:8b | block: module
- [ ] [HIGH] Inefficient data processing - No batching or parallel processing (performance) (function: process_data) lines [120, 160] | model: rnj-1:8b | block: module
- [ ] [MEDIUM] No input validation - No checks on configuration parameters (code_safety) (function: load_config) lines [200, 220] | model: rnj-1:8b | block: module