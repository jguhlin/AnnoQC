# AnnoQC TODO

Converted from COMPREHENSIVE_PLAN.md. All items start unchecked.

## TODO Tracker

- [ ] Stabilize `prepare` flow (download, validation, clustering artifacts).
- [ ] Implement initial homology + intrinsic scoring with JSONL/CSV outputs.
- [ ] Expose scoring thresholds/weights via TOML config (`scoring.*`).
- [ ] Add reciprocal coverage delta + fusion/split heuristic to scorecards.
- [ ] Add conserved-region analysis (MAFFT conserved fractions & pairwise identity).
- [ ] Extend intrinsic metrics (low complexity windows, ORF checks, entropy stats).
- [ ] Integrate Bevy ECS scheduler for workload balancing and task orchestration.
- [x] Minimal ECS scheduler scaffolding to compute per-gene metrics from FASTA + DIAMOND TSV.
- [ ] Support batch-mode DIAMOND invocations to reduce process churn.
- [ ] Implement taxonomy/domain evidence pillar (optional hmmscan integration).
- [x] Build golden fixtures + CLI smoke tests in CI.
- [ ] Document metrics & scoring in mdBook (`book/`) with examples.
- [ ] Package releases (Linux/macOS/Windows) and optional Docker image.
- [x] Build offline UniProt taxonomy cache + lineage resolver.
- [ ] Integrate hmmscan + Pfam taxonomy heuristics into scoring model.

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
- [ ] Schema versioning for outputs; CSV headers include tool versions.

## Caching & Resume

- [ ] Cache downloads (ETag/Last-Modified) and support resume.
- [ ] Step checkpoints (`.done` files) for idempotent prepare.

## Observability

- [ ] Human-readable logs; optional JSON log format.
- [ ] Progress and basic counters (genes/sec, queue sizes).
- [x] Step-level timing and resume logging for prepare stages.

## Completed in this iteration

- DIAMOND pre-run caching before processing.
- ECS scheduler scaffold consuming FASTA + DIAMOND TSV.
- Enriched homology metrics in outputs.
- Prepare checkpoints (`makedb`, `linclust`, `cluster`) and `.done` markers.


## Next Up (prioritized)

- [x] Schema versioning (run.json) and include basic config weights/hash info.
- [x] JSON log mode and progress counters in analyze loop.
- [x] Reciprocal coverage delta, fusion/split heuristics, and initial scoring weights.
- [ ] Optional MAFFT conserved-region metrics for top-N hits integrated into CSV.

## Packaging & Docs

- [ ] mdBook pages aligned for prepare/analyze/metrics & scoring.
- [ ] Release artifacts (Linux/macOS/Windows), optional Docker image.
