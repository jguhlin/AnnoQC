# Checklist: src/mafft.rs

## Severity counts
- MEDIUM: 12
- HIGH: 11
- LOW: 5
- INFO: 3
- VERIFY: 2

## Issues

- [ ] [HIGH] Arbitrary truncation thresholds without validation (statistical_integrity) (function: compute_alignment_metrics) lines [270, 290] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Re-allocates Vec<u8> per line (performance) (function: parse_fasta_sequences) lines [170, 180] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Assumes valid FASTA format without validation (code_safety) (function: parse_fasta_sequences) lines [170, 180] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Pairwise identity only compares first sequence to others, biasing results. (statistical) (function: compute_alignment_metrics) lines [95, 108] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Hard‑coded concordance thresholds may misclassify panels with larger variation. (statistical) (function: compute_alignment_metrics) lines [124, 147] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Builds entire FASTA into a Vec before piping; avoid extra allocation. (performance) (function: run_mafft_alignment) lines [45, 61] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Repeated per‑column filters cause O(cols * n) work; could be streamlined. (performance) (function: compute_consensus_gap_run) lines [207, 224] | model: nemotron-3-nano:30b | block: module
- [ ] [VERIFY] Unwraps position/rposition without guard; safety depends on external invariants. (code-safety) (function: parse_fasta_sequences) lines [71, 84] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] AlignmentEngine construction uses positional args that may not match documented order. (code-safety) (function: run_spoa_alignment) lines [165, 169] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing documentation for mafft_fast and mafft_max_jobs fields. (docstrings) (function: AlignerConfig) lines [30, 42] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Fallback behavior is only mentioned in comments, not in doc comment. (docstrings) (function: run_alignment_for_panel) lines [176, 176] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] The MAFFT binary path is hardcoded and may not be portable across systems. This can lead to failures if the binary is not in the expected location. Use a more flexible approach to locate the MAFFT binary. (performance) (function: run_mafft_alignment) lines [50, 60] | model: rnj-1:8b | block: module
- [ ] [MEDIUM] The MAFFT command construction uses string concatenation, which can be inefficient and error-prone. Consider using a more robust method to build the command. (performance) (function: run_mafft_alignment) lines [50, 60] | model: rnj-1:8b | block: module
- [ ] [LOW] The function lacks detailed documentation about the expected behavior, error handling, and return values. Add docstrings to improve maintainability. (docstrings) (function: run_mafft_alignment) lines [50, 60] | model: rnj-1:8b | block: module
- [ ] [HIGH] Divergence ratio may be misleading when panel pairwise identity is near zero (statistical_integrity) (function: compute_alignment_metrics) lines [275, 280] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Redundant ratios computed per column (performance) (function: compute_consensus_gap_run) lines [303, 315] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Bypasses mafft_bin check in tests (code_safety) (function: run_mafft_alignment) lines [60, 68] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Divergence ratio may yield inf/nan when panel identity is zero; interpretation questionable. (statistical) (function: compute_alignment_metrics) lines [140, 150] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] O(N^2) pairwise length checks per column can be costly for large panels. (performance) (function: validate_spoa_alignment) lines [78, 94] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc comment describing query_gap semantics and return behavior. (documentation) (function: compute_consensus_gap_run) lines [260, 270] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] The MAFFT binary path is hardcoded and may not be portable across systems. Use environment variables or configuration files to specify the path. (performance) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module
- [ ] [MEDIUM] The MAFFT binary is executed synchronously, which can block the main thread. Consider using asynchronous execution or a thread pool for non-blocking I/O. (performance) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module
- [ ] [LOW] The function lacks documentation for the `mafft_bin` parameter and return value. Add docstrings to clarify expected behavior and error handling. (docstrings) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module
- [ ] [HIGH] O(N^2 * L) pairwise identity computation is inefficient (performance) (function: compute_alignment_metrics) lines [222, 233] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] No validation of mafft_bin path before spawning process (code_safety) (function: run_mafft_alignment) lines [100, 105] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] divergence_ratio may divide by zero when panel size ≤2, producing undefined scores. (statistical) (function: compute_alignment_metrics) lines [70, 85] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Repeated small Vec allocations per FASTA line cause allocation overhead. (performance) (function: parse_fasta_sequences) lines [462, 479] | model: nemotron-3-nano:30b | block: module
- [ ] [VERIFY] Potential usize overflow in non_gap count calculation on very long sequences. (safety) (function: validate_spoa_alignment) lines [340, 357] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc explaining strategy_used overwrite and returned metric meanings. (documentation) (function: metrics_from_alignment) lines [260] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] No doc describing query_gap flag semantics and 70% threshold rule. (documentation) (function: compute_consensus_gap_run) lines [384, 398] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] The MAFFT binary path is hardcoded and may not exist in the environment, leading to potential failures. Use a more robust method to locate the binary or provide a configurable path. (performance) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module
- [ ] [MEDIUM] The MAFFT binary is executed for each alignment job, which can be slow for large datasets. Consider caching results or using a more efficient alignment method. (performance) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module
- [ ] [LOW] The function lacks documentation for the `mafft_bin` parameter, which is critical for understanding how to configure the binary path. (docstrings) (function: run_mafft_alignment) lines [45, 55] | model: rnj-1:8b | block: module