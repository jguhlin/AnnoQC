# Checklist: src/explain.rs

## Severity counts
- HIGH: 2
- MEDIUM: 2
- INFO: 1

## Issues

- [ ] [HIGH] Potential panic on unwrap_or with unknown or NA values (code_safety) (function: print_breakdown) lines [27, 32, 37, 42, 47, 52] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Inefficient line-by-line JSON parsing with redundant allocations (performance) (function: explain_gene) lines [14, 25] | model: qwen3-coder:latest | block: module
- [ ] [INFO] Collects warnings into a temporary Vec before joining; could avoid allocation by directly iterating and building the string. (performance) (function: print_breakdown) lines [124, 138] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing docstring describing parameters, return value, and behavior of early exit; also missing description of printed breakdown format. (documentation) (function: explain_gene) lines [1, 20] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Potential panic on unwrap_or with unknown/NA fallbacks (code_safety) (function: explain_gene) lines [35, 50, 57, 63, 69] | model: qwen3-coder:latest | block: module