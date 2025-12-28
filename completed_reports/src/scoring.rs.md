# Checklist: src/scoring.rs

## Severity counts
- HIGH: 5
- MEDIUM: 5
- LOW: 3
- INFO: 1

## Issues

- [ ] [HIGH] Early min(1.0) operations distort weighted contributions of components (statistical_integrity) (function: compute_homology_score) lines [12, 22] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Missing seed control for deterministic behavior (statistical_integrity) (function: compute_intrinsic_score) lines [24, 29] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Returns NaN when congruence_score is NaN; should clamp or validate input. (statistical_integrity) (function: compute_taxonomy_score) lines [70, 75] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Borrow of Option<&DiamondHitStats> has no lifetime constraint; ensure data outlives call. (code_safety) (function: compute_homology_score) lines [1, 20] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Missing rustdoc documentation for parameters and return value. (docstrings) (function: compute_intrinsic_score) lines [45, 60] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Consider caching the result of `s.coverage_delta` if it's expensive to compute, as it's used multiple times in the function. (performance) (function: compute_homology_score) lines [1, 30] | model: rnj-1:8b | block: module
- [ ] [HIGH] Final score not normalized or calibrated across inputs (statistical_integrity) (function: compute_homology_score) lines [12, 22] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] No handling for class imbalance in congruence score (statistical_integrity) (function: compute_taxonomy_score) lines [30, 32] | model: qwen3-coder:latest | block: module
- [ ] [LOW] Unnecessary final clamp after already clamped values (performance) (function: compute_genomic_score) lines [70, 72] | model: qwen3-coder:latest | block: module
- [ ] [LOW] Potential division by zero if introns_total is not initialized (code_safety) (function: compute_genomic_score) lines [66, 67] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Multiplicative penalties may overpenalize due to non-linear interaction (statistical_integrity) (function: compute_genomic_score_with_cfg) lines [130, 150] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Hardcoded linear ramp assumes divergence_ratio is bounded in [0.5, 0.8] (statistical_integrity) (function: compute_divergence_score) lines [100, 108] | model: qwen3-coder:latest | block: module
- [ ] [LOW] Redundant floating-point divisions and casts can be hoisted; pre‑compute denominators. (performance) (function: compute_homology_score) lines [14, 23] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc comment describing parameters and return value. (documentation) (function: compute_taxonomy_score) lines [27, 30] | model: nemotron-3-nano:30b | block: module