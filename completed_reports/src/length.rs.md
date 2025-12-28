# Checklist: src/length.rs

## Severity counts
- MEDIUM: 4
- HIGH: 1
- LOW: 1

## Issues

- [ ] [HIGH] Potential division by zero in Z-score normalization (statistical_integrity) (function: compute_length_consistency) lines [40, 41] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Redundant sorting in median and MAD computation (performance) (function: compute_length_consistency) lines [17, 24] | model: qwen3-coder:latest | block: module
- [ ] [LOW] Use of unwrap_or with partial_cmp may panic on NaN (code_safety) (function: compute_length_consistency) lines [17, 24] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Hard-coded ratio cutoffs (.8, 1.2) lack statistical justification and may misclassify borderline cases. (statistical_integrity) (function: compute_length_consistency) lines [70, 75] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Two full sorts and intermediate allocations cause O(n log n) overhead in hot path. (performance) (function: compute_length_consistency) lines [120, 150] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] The function performs multiple sorting operations on the input data, which can be optimized for better performance. The sorting operations can be replaced with more efficient statistical calculations to avoid unnecessary computations. (performance) (function: compute_length_consistency) lines [40, 80] | model: rnj-1:8b | block: module