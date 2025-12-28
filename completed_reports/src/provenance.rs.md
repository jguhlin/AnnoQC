# Checklist: src/provenance.rs

## Severity counts
- MEDIUM: 3
- INFO: 3
- LOW: 1
- HIGH: 1

## Issues

- [ ] [MEDIUM] 64KB buffer may not align with OS page size for optimal I/O (performance) (function: filehash_xx64) lines [13, 18] | model: qwen3-coder:latest | block: module
- [ ] [LOW] Repeated `format!` call in hot path causes heap allocation (performance) (function: filehash_xx64) lines [19] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] No explicit error handling for `file.read()` returning `None` (code_safety) (function: filehash_xx64) lines [15] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Function and parameters lack documentation; add a doc comment clarifying path, seed, and hex output. (docstrings) (function: filehash_xx64) lines [5, 30] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Consider using BufReader to reduce syscalls in the read loop; replace manual buffer with BufReader::new(File::open(path)?)?.read_to_end(&mut hasher). (performance) (function: filehash_xx64) lines [14, 20] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] Potential resource leak — File handle not explicitly closed in error case (Code Safety) (function: filehash_xx64) lines [40, 55] | model: rnj-1:8b | block: module
- [ ] [INFO] Unnecessary buffer size — Reduce buffer size for better cache efficiency (Performance) (function: filehash_xx64) lines [50, 60] | model: rnj-1:8b | block: module
- [ ] [INFO] Missing documentation — Add docstring explaining purpose, input/output, and error handling (Docstrings) (function: filehash_xx64) lines [1, 10] | model: rnj-1:8b | block: module