# Checklist: src/refprot.rs

## Severity counts
- MEDIUM: 3
- HIGH: 2
- INFO: 1
- LOW: 1

## Issues

- [ ] [HIGH] Unchecked indexing in line splitting may panic (code_safety) (function: parse_readme) lines [40] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Repeated string cloning and joining may cause allocations (performance) (function: parse_readme) lines [30, 50] | model: qwen3-coder:latest | block: module
- [ ] [MEDIUM] Unbounded loop may cause performance degradation (performance) (function: select_by_taxon) lines [65, 75] | model: qwen3-coder:latest | block: module
- [ ] [HIGH] Fallback division default may misclassify organisms, leading to biased proteome sets. (statistical_integrity) (function: parse_readme) lines [30, 45] | model: nemotron-3-nano:30b | block: module
- [ ] [MEDIUM] O(N·D) taxonomy climb per entry; can be slow for large inputs. (performance) (function: select_by_taxon) lines [84, 95] | model: nemotron-3-nano:30b | block: module
- [ ] [INFO] Potential shared mutable access in resolver.parent_of across threads. (code_safety) (function: select_by_taxon) lines [70, 90] | model: nemotron-3-nano:30b | block: module
- [ ] [LOW] Docstring does not note truncation at max_proteomes. (docstrings) (function: select_by_taxon) lines [70, 80] | model: nemotron-3-nano:30b | block: module