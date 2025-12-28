# Checklist: src/metrics.rs

## Severity counts
- HIGH: 2
- MEDIUM: 2

## Issues

- [ ] [HIGH] O(n^2) sliding window uniqueness check; use hash set or bit manipulation (performance) (function: compute_intrinsic) lines [70, 85] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Redundant to_ascii_uppercase() calls in inner loop (performance) (function: compute_intrinsic) lines [73, 78] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Sliding-window low-complexity calculation is O(n·w) with per-window allocation, causing high CPU for large sequences. (performance) (function: compute_intrinsic) lines [38, 57] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc comment describing metric semantics and assumptions. (documentation) (function: compute_intrinsic) lines [1, 5] | model: nemotron-3-nano:30b | block: module