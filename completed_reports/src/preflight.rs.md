# Checklist: src/preflight.rs

## Severity counts
- MEDIUM: 7
- HIGH: 3
- INFO: 2
- LOW: 1

## Issues

- [ ] [HIGH] Unvalidated binary paths may lead to command injection (code_safety) (function: check_tool_version) lines [12, 20] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Silent failure in version checks may lead to incorrect tool state (code_safety) (function: preflight) lines [27, 37] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Redundant string conversions may cause unnecessary allocations (performance) (function: check_tool_version) lines [12, 20] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Potential injection risk from unvalidated binary paths (safety) (function: check_tool_version) lines [40, 45] | model: rnj-1:8b | block: module
- [ ] [MEDIUM] Missing UTF-8 error handling for tool output (safety) (function: check_tool_version) lines [50, 60] | model: rnj-1:8b | block: module
- [ ] [INFO] Missing documentation for struct fields (docstrings) (function: ToolVersions) lines [10, 20] | model: rnj-1:8b | block: module
- [ ] [HIGH] Assumes UTF-8 output; may panic or corrupt data on non-UTF-8 stdout/stderr (code_safety) (function: check_tool_version) lines [13, 22] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Does not validate diamond_bin path is executable before command execution (code_safety) (function: check_diamond_db) lines [50, 51] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] No documentation for arguments or return value; makes usage unclear. (docstrings) (function: check_tool_version) lines [9, 20] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Assumes UTF-8 output; may panic on non-UTF-8 binaries (code_safety) (function: check_tool_version) lines [12, 20] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Does not validate binary path is executable (code_safety) (function: check_diamond_db) lines [33, 40] | model: qwen3-coder:latest | block: module
- [ ] [INFO] Repeated allocation of error string may cause heap churn on large stderr; replace with bounded write. (performance) (function: check_diamond_db) lines [61, 63] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Missing doc comments for return semantics and side effects, risking misuse. (documentation) (function: *multiple*) lines [9, 30] | model: nemotron-3-nano:30b | block: module