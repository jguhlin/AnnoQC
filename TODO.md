# AnnoQC TODO

Converted from COMPREHENSIVE_PLAN.md. All items start unchecked.

## Status Snapshot (2025-11-06)

- Core CLI operational; run.json, qc_report.jsonl, qc_summary.csv emitted.
- DIAMOND pre-run caching; retries; optional chunked mode; homology metrics incl. coverage_delta/ratio and fusion flag.
- ECS minimal scheduler with JSON/text progress.
- Intrinsic metrics; MAFFT metrics (optional) in JSONL and CSV.
- Scoring: homology+intrinsic weights; thresholds; CSV final_score/classification; taxonomy weight plumbed (default 0.0).
- Prepare: checkpointed makedb/linclust/cluster/recluster placeholder; `--resume`; JSON logs.
- Manifest: schema_version, tool versions, file hashes, config snapshot.
- HMMER/Pfam scaffolding (gated) with domtblout parser; JSONL domains summary.

## TODO Tracker

- [x] Stabilize `prepare` flow (download, validation, clustering artifacts).
- [x] Implement initial homology + intrinsic scoring with JSONL/CSV outputs.
- [x] Expose scoring thresholds/weights via TOML config (`scoring.*`).
- [x] Add reciprocal coverage delta + fusion/split heuristic to scorecards.
- [x] Add conserved-region analysis (MAFFT conserved fractions & pairwise identity).
- [x] Extend intrinsic metrics (low complexity windows, ORF checks, entropy stats).
- [x] Integrate Bevy ECS scheduler for workload balancing and task orchestration.
- [x] Support batch-mode DIAMOND invocations + Auto heuristic + chunk logs.
- [x] Implement taxonomy/domain evidence pillar (optional hmmscan integration).
- [x] Build golden fixtures + CLI smoke tests in CI.
- [x] Document metrics & scoring in mdBook (`book/`) with examples.
- [ ] Package releases (Linux/macOS/Windows) and optional Docker image.
- [x] Build offline UniProt taxonomy cache + lineage resolver.
- [ ] Integrate hmmscan + taxonomy congruence heuristics into scoring model (future).

## Goals

- [ ] Fast, scalable QC for gene annotations with evidence-based scoring.
- [ ] ECS scheduling to balance heterogeneous workloads.
- [ ] Transparent per-gene scorecards with metrics for downstream pipelines.

## Architecture & Data Flow

- [ ] CLI with `prepare` and `analyze` subcommands.
- [x] Implemented `prepare` (makedb only) and `analyze`.
- [ ] External tools wired: DIAMOND, optional MAFFT/HMMER.
- [x] DIAMOND wired: `--version`, `blastp` pre-run, `makedb` in `prepare`.
- [x] DIAMOND TSV parsed for top-hit stats (bitscore, evalue, coverage, density).
- [ ] FASTA handling via `needletail` (gz supported).
- [x] FASTA parse for ids/lengths via needletail.
- [ ] `prepare`: makedb → cluster → realign → recluster artifacts.
- [x] Implement makedb + linclust + cluster with `.done` checkpoints; recluster placeholder.
- [x] Add `--resume` support and progress logs for all prepare steps (checkpointed).
- [ ] `analyze`: read input → schedule tasks → parse DIAMOND → compute metrics → emit outputs.
- [x] Read input + schedule tasks + emit minimal outputs.
- [x] Emit enriched homology fields in JSONL and CSV.
- [x] Add final_score and classification to CSV; JSONL contains per-pillar scores and final_score.

## Testing & CI

- [x] Unit tests for parsing/metrics/report formatting.
- [x] Fixture-based tests for pipelines and scoring.
- [x] Golden tests for summaries, property tests as needed.
- [x] CI: fmt, clippy -D warnings, tests, tiny integration with DIAMOND stub.

## Performance & Reliability

- [ ] Streaming I/O and bounded memory usage.
- [ ] Expose threads/approx-id/member-cover; sensible defaults.
- [ ] Retry transient DIAMOND errors; concise diagnostics.

## Reproducibility

- [ ] Run manifest (tool versions, config snapshot, checksums) at `results/run.json`.
- [x] Basic run manifest with tool versions (diamond/mafft/hmmscan) and inputs.
- [x] Add file hashes for FASTA and DB (xx64).
- [x] Manifest includes schema_version and a config snapshot (resolved settings & weights).
- [ ] Schema versioning for outputs; CSV headers include tool versions.

## Caching & Resume

- [ ] Cache downloads (ETag/Last-Modified) and support resume.
- [ ] Step checkpoints (`.done` files) for idempotent prepare.

## Observability

- [x] Human-readable logs; optional JSON log format.
- [x] Progress and basic counters (genes/sec, queue sizes).
- [x] Step-level timing; analyze emits step_start/step_finish JSON; metrics sidecar.

## Completed in this iteration

- DIAMOND pre-run caching; retries; optional chunked mode.
- ECS scheduler scaffold consuming FASTA + DIAMOND TSV; progress logs.
- Enriched homology metrics; coverage delta/ratio; fusion/split flag.
- Prepare checkpoints (`makedb`, `linclust`, `cluster`) with `--resume`.
- Scoring weights + thresholds; CSV final_score/classification.
- Manifest schema_version + config snapshot; file hashes.
- MAFFT metrics (JSONL and CSV); HMMER/Pfam scaffolding and domains JSONL.


## Next Up (prioritized)

- [x] Schema versioning (run.json) and include basic config weights/hash info.
- [x] JSON log mode and progress counters in analyze loop.
- [x] Reciprocal coverage delta, fusion/split heuristics, and initial scoring weights.
- [x] Optional MAFFT conserved-region metrics for top-N hits integrated into CSV.
- [x] Taxonomy pillar (Phase 2): lineage resolution + taxonomy_score and outputs
  - Wire `TaxonomyResolver::from_sources(cache, reference_fasta, taxdump_dir)` and resolve top-hit accessions.
  - JSONL: taxonomy {status, taxid, name, lineage[]} + scalar `taxonomy_score` (presence now, congruence later).
  - CSV: add taxonomy_score when enabled; taxonomy_status already switches enabled/disabled.
  - Scoring: include in `score_components` and weighted `final_score` when `[scoring.weights].taxonomy > 0`.
  - Tests: add resolver unit tests with tiny headers/fixtures; smoke test unchanged when disabled.

- [x] Pfam/HMMER (Phase 2): domains_score + architecture summary
  - Batch/threaded `hmmscan`; parse domtblout (parser present).
  - JSONL: expand domains list (capped) and `domains_score`; config gates execution.
  - Tests: tiny domtblout fixtures; no external calls in CI.

- [x] DIAMOND batch mode: Auto heuristics and chunk progress logs
  - Implement `DiamondMode::Auto` decision by input size; JSON log per chunk (genes, secs, rate).
  - Improve error reporting; capture tool stderr in logs.

- [x] Observability & Logs
  - Analyze start/finish JSON events; durations and counters.
  - Optional text progress bar (off in JSON mode); metrics sidecar implemented.

- [ ] Docs & Packaging
  - mdBook: add “Taxonomy & Domains” and “Scoring” pages; realistic examples.
  - Quickstart outputs updated; optional Docker with DIAMOND/HMMER via Pixi.

## Implementation Notes

- Taxonomy (Phase 2)
  - Source: `src/taxonomy.rs`; call after DIAMOND stats. Cache from FASTA or TSV; use `new_taxdump` for full lineage.
  - Outputs: JSONL taxonomy block + taxonomy_score; CSV taxonomy_score when enabled; adjustable scoring weight.

- HMMER/Pfam (Phase 2)
  - Source: `src/hmmer.rs` (present). Add executor for batch/threaded hmmscan, then compute `domains_score`.
  - Outputs: JSONL domains list and score; optional CSV columns later.

- DIAMOND
  - Chunked mode implemented; add `Auto` heuristics and per-chunk JSON logs.

- Manifest
  - `schema_version=1.0`; resolved config snapshot; consider adding classification thresholds explicitly.

- Docs
  - `book/src/analyze.md` updated; next add real JSONL/CSV samples from fixtures and a new “Taxonomy & Domains” page.

## Packaging & Docs

- [ ] mdBook pages aligned for prepare/analyze/metrics & scoring.
- [x] Update analyze docs to reflect new JSONL/CSV fields and manifest contents.
- [ ] Release artifacts (Linux/macOS/Windows), optional Docker image.
