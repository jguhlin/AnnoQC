# DIAMOND & Taxonomy Refactor Migration

## Goal
Modernize the `analyze` pipeline to (1) support a single-pass DIAMOND workflow with ECS-only post-processing, and (2) integrate Pfam-driven `hmmscan` evidence into the taxonomy pillar.

## Checklist

1. [ ] Augment CLI/config (`AnalyzeArgs`, `AnalyzeConfig`, `config.example.toml`) with `diamond_mode` and taxonomy file paths; stage utility functions (canonical accession reuse).
2. [ ] Refactor pipeline loading: pre-parse FASTA into `Vec<GeneJob>`, compute taxonomy resolver once, introduce single-shot execution path (write/query temp FASTA, parse output) and shared output writer.
3. [ ] Introduce new DIAMOND result parsing helpers reused by batch mode; adjust batch executor to consume shared helpers and remove stdin streaming.
4. [ ] Wire optional `hmmscan` execution using Pfam metadata: run per gene when enabled, parse domain hits, feed into taxonomy scoring (support config paths).
5. [ ] Update docs/tests: `config.example.toml`, `COMPREHENSIVE_PLAN.md`, smoke test to exercise single-shot path, ensure `cargo fmt`, `cargo check`, `cargo test` pass.

Mark items as completed while implementing. Delete this file after all tasks are checked.
