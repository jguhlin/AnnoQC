# Checklist: tests/analyze_smoke.rs

## Severity counts
- MEDIUM: 5
- HIGH: 3
- INFO: 1
- LOW: 1

## Issues

- [ ] [HIGH] Missing statistical validation of QC metrics (statistical_integrity) (function: analyze_smoke_produces_outputs) lines [45, 75] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unsafe unwrap on env! macro (code_safety) (function: analyze_smoke_produces_outputs) lines [50, 50] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Inefficient JSON parsing loop (performance) (function: analyze_smoke_produces_outputs) lines [79, 90] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Missing validation of QC metrics for data leakage or improper CV (statistical_integrity) (function: analyze_smoke_produces_outputs) lines [100, 140] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Inefficient JSON parsing in test loop without batching (performance) (function: analyze_smoke_produces_outputs) lines [120, 135] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Unsafe unwrap() on env var and binary path (code_safety) (function: analyze_smoke_produces_outputs) lines [100, 105] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Hard-coded expectation of exactly two scorecards and specific gene IDs; test brittle to fixture changes. (statistical_integrity) (function: analyze_smoke_produces_outputs) lines [70, 95] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Collects all JSONL records into a Vec<Value> causing unnecessary memory allocation; could stream. (performance) (function: analyze_smoke_produces_outputs) lines [84, 102] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] `unwrap()` on Option values may panic on malformed input; replace with proper error handling. (code_safety) (function: analyze_smoke_produces_outputs) lines [115, 130] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing documentation for parameters, returns, and side effects. (docstrings) (function: compile_diamond_stub) lines [15, 45] | model: nemotron-3-nano:30b | block: module