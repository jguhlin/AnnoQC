# Checklist: src/rhai_rules.rs

## Severity counts
- MEDIUM: 5
- HIGH: 3
- LOW: 2
- VERIFY: 1

## Issues

- [ ] [HIGH] Potential panic on file_stem() usage (code_safety) (function: load) lines [44, 47] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Repeated Dynamic clone in loop (performance) (function: run) lines [75, 85] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Arc<Engine> clone in hot path (code_safety) (function: run) lines [75, 85] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Cloning input_dyn inside rule loop allocates per iteration causing unnecessary heap pressure; reuse reference or move clone outside. (performance) (function: RhaiRuntime::run) lines [68, 85] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Missing doc comment describing parameters and return semantics; unclear contract for callers. (documentation) (function: RhaiRuntime::load) lines [20, 38] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Potential panic on file_stem() usage (code_safety) (function: RhaiRuntime::run) lines [35] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Repeated cloning of Dynamic input increases memory pressure (performance) (function: RhaiRuntime::run) lines [45, 50] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Missing seed control for stochastic scripts may cause nondeterministic outputs. (statistical_integrity) (function: RhaiRuntime::run) lines [70, 85] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Creates a new Scope per rule causing repeated allocations; can be optimized. (performance) (function: RhaiRuntime::run) lines [61, 73] | model: nemotron-3-nano:30b | block: module
- [ ] [VERIFY] Engine wrapped in Arc may not be Send+Sync; concurrent use unproven. (code_safety) (function: RhaiRuntime struct) lines [15, 20] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Missing doc comments for parameters, returns, and side effects. (docstrings) (function: RhaiRule,RhaiRuntime public methods) lines [18, 45] | model: nemotron-3-nano:30b | block: module