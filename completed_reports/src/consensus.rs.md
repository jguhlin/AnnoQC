# Checklist: src/consensus.rs

## Severity counts
- HIGH: 6
- MEDIUM: 4
- LOW: 2

## Issues

- [ ] [HIGH] Uses unwrap_or on max() of lengths, which may silently miscompute coverage if all lengths are zero (statistical_integrity) (function: AggregatedHitBuilder::finalize) lines [110, 117] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Performs redundant sorting and span merging on potentially large hit sets (performance) (function: union_span_len) lines [158, 175] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Assumes qlen and slen are always positive, but filter(|&l| l > 0) may miss edge cases (code_safety) (function: AggregatedHitBuilder::finalize) lines [122, 123] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Unchecked cfg.len_ratio_tolerance may cause invalid windows leading to overflow or incorrect filtering. (statistical_integrity) (function: select_panel_with_result) lines [130, 150] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] median function sorts slice without handling NaN, may panic on malformed len_ratio values. (statistical_integrity) (function: compute_len_stats) lines [235, 245] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Clones each DiamondHitRow into a Vec inside builder, causing O(n) allocations and copies. (performance) (function: aggregate_hits_by_subject) lines [70, 90] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Sorts ratios vector twice for median and deviations, incurring extra O(m log m) overhead. (performance) (function: compute_len_stats) lines [235, 250] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Panics if sseqid is empty due to unconditional parts.next().unwrap(); can be triggered by malformed input. (code_safety) (function: subject_root) lines [340, 345] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Parses h.evalue with unwrap_or(1.0); panics if evalue string is not a valid float. (code_safety) (function: run_phase) lines [260, 265] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Missing detailed doc comment describing cfg parameters and backfill behavior. (docstrings) (function: select_panel_with_result) lines [73, 95] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Doc comment lacks explanation of cap semantics and effect when cap==0. (docstrings) (function: cap_refprot_hits) lines [115, 125] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Doc comment omits description of len_shift and extra parameters effect on window boundaries. (docstrings) (function: adjust_phase_window) lines [124, 135] | model: nemotron-3-nano:30b | block: module