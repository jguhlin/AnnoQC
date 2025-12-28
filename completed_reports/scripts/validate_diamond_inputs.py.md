# Checklist: scripts/validate_diamond_inputs.py

## Severity counts
- HIGH: 3
- MEDIUM: 3
- LOW: 1

## Issues

- [ ] [HIGH] Invalid residue detection logic misclassifies valid IUPAC codes (statistical_integrity) (function: analyze_fasta) lines [60, 65] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Redundant line.upper() calls in loop (performance) (function: analyze_fasta) lines [60, 65] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Unsafe subprocess call allows command injection (code_safety) (function: check_db) lines [80, 85] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Missing gzip integrity validation (code_safety) (function: open_text) lines [20, 25] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Per-character Counter updates in hot loop cause Python overhead (L73-84) (performance) (function: analyze_fasta) lines [73, 84] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] subprocess returncode not raised; may hide diamond errors (L95-103) (safety) (function: check_db) lines [95, 103] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Missing docstrings for validation logic (L30-57, L98-107) (documentation) (function: ['analyze_fasta', 'check_db']) lines [30, 57, 98, 107] | model: nemotron-3-nano:30b | block: module