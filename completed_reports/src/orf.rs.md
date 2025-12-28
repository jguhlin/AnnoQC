# Checklist: src/orf.rs

## Severity counts
- MEDIUM: 3
- HIGH: 2

## Issues

- [ ] [HIGH] Potential panic on empty or malformed ORF input due to unwrap_or (code_safety) (function: translate_longest_orf) lines [54] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Inefficient string conversion and repeated allocations (performance) (function: translate_longest_orf) lines [47, 52] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unsafe indexing without bounds check (code_safety) (function: translate_longest_orf) lines [46] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Builds intermediate AA string and splits on '*', leading to O(N) allocations; replace with direct longest‑segment extraction. (performance) (function: translate_longest_orf) lines [70, 95] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Doc comment omits that the returned path is a temporary file requiring explicit cleanup. (documentation) (function: translate_nt_fasta_to_protein) lines [6, 15] | model: nemotron-3-nano:30b | block: module