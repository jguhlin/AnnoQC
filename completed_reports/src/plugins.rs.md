# Checklist: src/plugins.rs

## Severity counts
- MEDIUM: 5
- HIGH: 4
- INFO: 1

## Issues

- [ ] [HIGH] Potential panic on invalid file path (code_safety) (function: PluginDefinition::load) lines [169, 175] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unnecessary clone of Arc<Vec<u8>> (performance) (function: run_plugin) lines [190, 190] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] No error handling for PluginOutput deserialization (code_safety) (function: run_plugin) lines [195, 195] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Clones wasm_bytes via .clone() inside Wasm::data, causing extra allocation per call; use slice directly if possible. (performance) (function: PluginDefinition::load) lines [24, 30] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Missing doc comment describing input schema and returned PluginOutput fields; adds maintenance risk. (documentation) (function: run_plugin) lines [55, 68] | model: nemotron-3-nano:30b | block: module
- [ ] [HIGH] Potential panic on invalid file path (code_safety) (function: PluginDefinition::load) lines [75, 79] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unnecessary clone of Arc<Vec<u8>> (performance) (function: run_plugin) lines [89, 90] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unsafe deserialization without input limits (code_safety) (function: run_plugin) lines [94, 94] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Repeated Manifest/Plugin construction in hot path causing per-call allocation; cache or reuse to avoid O(N) overhead. (performance) (function: run_plugin) lines [45, 53] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Potential aliasing of wasm_bytes across threads not examined; requires evidence of disjoint mutable access. (code_safety) | model: nemotron-3-nano:30b | block: module