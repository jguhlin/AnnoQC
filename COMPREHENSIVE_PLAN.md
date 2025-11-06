# Comprehensive Plan

AnnoQC is a Rust CLI that evaluates gene/protein annotations at scale by combining fast similarity search (DIAMOND) with an ECS scheduler (Bevy 0.17) to distribute heterogeneous workloads efficiently. It outputs reproducible QC metrics and summaries (CSV/JSON), aiming for high throughput and clear diagnostics on commodity hardware.

## TODO Tracker

- [x] Stabilize `prepare` flow (download, validation, clustering artifacts).
- [x] Implement initial homology + intrinsic scoring with JSONL/CSV outputs.
- [x] Expose scoring thresholds/weights via TOML config (`scoring.*`).
- [x] Add reciprocal coverage delta + fusion/split heuristic to scorecards.
- [x] Add conserved-region analysis (MAFFT conserved fractions & pairwise identity).
- [x] Extend intrinsic metrics (low complexity windows, ORF checks, entropy stats).
- [x] Integrate Bevy ECS scheduler for workload balancing and task orchestration.
- [x] Support batch-mode DIAMOND invocations to reduce process churn.
- [x] Implement taxonomy/domain evidence pillar (optional hmmscan integration).
- [x] Build golden fixtures + CLI smoke tests in CI.
- [x] Document metrics & scoring in mdBook (`book/`) with examples.
- [x] Package releases (Linux/macOS/Windows) and optional Docker image.
- [x] Build offline UniProt taxonomy cache + lineage resolver.
- [ ] Integrate hmmscan + Pfam taxonomy heuristics into scoring model.

## Goals

- Fast, scalable QC for gene annotations using DIAMOND for similarity and Rust for throughput.
- Balance heterogeneous workloads (short/long genes; few/many matches) via Bevy ECS scheduling.
- Produce a transparent, evidence-based quality score for each gene, with detailed metrics suitable for downstream pipelines and reports.

## Architecture Overview

- CLI (`prepare`, `analyze`) in `src/main.rs`.
- External tools: DIAMOND (`makedb`, `cluster`, `realign`, `recluster`, later `blastp`). Optional `hmmscan` for advanced analysis.
- FASTA handling: `needletail` for parsing/validation (gz supported).
- ECS runtime (Bevy 0.17): systems orchestrate per‑gene tasks; dynamic work batching keeps all cores busy.

## Data Flow

1) Prepare
   - Download SwissProt (`uniprot_sprot.fasta.gz`).
   - Validate FASTA (sample read with needletail).
   - DIAMOND: `makedb` → `cluster` (approx‑id 50, member‑cover 80, threads 16) → `realign` → `recluster` → artifacts: `uniprot_sprot.dmnd`, `clusters*`.
2) Analyze
   - Read input FASTA; queue genes into ECS/worker system.
   - For each gene, schedule tasks for homology search and intrinsic sequence analysis.
   - Parse DIAMOND outputs, compute intrinsic metrics, and aggregate evidence into a final scorecard.
   - Emit metrics to CSV/JSONL plus optional per-gene reports.

## Analyses & Scoring Model

AnnoQC evaluates gene models using a multi-pillar evidence-based approach. Each gene receives a final score (0-1) derived from weighted sub-scores from each evidence pillar. The output includes both the final score and the underlying data for full transparency.

* **Pillar 1: Homology Evidence**
  - **Hit Significance:** Bitscore, e-value, and **bitscore density** (`bitscore / alignment_length`) of top hits against SwissProt.
  - **Alignment Quality:** **Reciprocal coverage analysis** to detect potential gene fusions/splits. Penalization for high gap counts or internal stop codons.
  - **Conserved Regions:** Alignment of top hits (MAFFT or alternative) to flag inconsistencies in highly conserved areas.

* **Pillar 2: Intrinsic Sequence Quality**
  - **ORF Integrity:** Checks for canonical start/stop codons and absence of in-frame stop codons.
  - **Sequence Composition:** Low-complexity region detection (sdust-like), GC content.
  - **Genomic Context (when available):** Intron/exon length distributions and validation of canonical splice sites (e.g., GT-AG).

* **Pillar 3: Advanced Evidence (Future/Optional)**
  - **Protein Domain Architecture:** Identification of known protein domains using `hmmscan` (HMMER) against a database like Pfam. Checks for domain completeness and consistency with homologs.
  - **Taxonomic Congruence:** Analysis of the taxonomic lineage of top hits to flag potential contamination or horizontal gene transfer events.

## ECS Load Balancing Design

- Entities: `Gene`, `MatchTask`, `AnalysisTask`.
- Components: `Seq`, `GeneId`, `WorkSize`, `Status`, `HomologyResult`, `IntrinsicMetrics`, `FinalScorecard`.
- Resources: `Config`, `DiamondPaths`, `ThreadPool`, `Reporter`, `ScoringRubric`.
- Systems:
  - **Intake:** Chunk input genes and estimate work size.
  - **Schedule:** Push tasks to sized queues; prefer small tasks when queue is imbalanced.
  - **ExecuteHomology:** Spawn external DIAMOND jobs; parse streams incrementally.
  - **CalculateIntrinsics:** Compute ORF, complexity, and composition metrics from the sequence.
  - **ScoreAggregator:** Queries for entities with completed `HomologyResult` and `IntrinsicMetrics` components, applies a configurable scoring rubric, and adds a `FinalScorecard`.
  - **Report:** Write structured outputs from `FinalScorecard` components; progress/logging.

## CLI & Config

- Introduce optional `--config <file>.toml` to set thresholds, scoring weights, queues, and tool paths; CLI flags override config.
- Extend `analyze`:
  - `--fasta <path>` (input), `--db uniprot_sprot.dmnd`, `--threads N`, `--approx-id P`, `--member-cover P`, `--out <dir>`.
  - Optional `--mode {centroid, members, auto}` and `--top N`.
- Add `--log-level` via `RUST_LOG` and friendly defaults.

## Testing Strategy

- Unit tests for parsing, metrics, and report formatting (no network/DIAMOND).
- Fixture‑based tests with tiny FASTA/TSV to verify pipelines and scoring logic.
- Golden tests for summaries; property tests for metrics where applicable.
- Integration tests behind `--ignored` that require DIAMOND and local cache.

## Performance & Reliability

- Use streaming readers and bounded memory; avoid loading whole FASTA when possible.
- Expose thread/approx-id/member-cover as CLI options; sensible defaults.
- Retry transient DIAMOND errors; surface concise diagnostics.

## Reproducibility & Provenance

- Record inputs, SwissProt release date (or file checksum), DIAMOND version/flags, and CLI/config snapshot in a run manifest (`results/run.json`).
- Embed schema version in all outputs; include tool versions as headers in CSV.

## Caching & Resume

- Cache downloads with ETag/Last-Modified; skip if unchanged. Support resume for partial downloads.
- Use step checkpoints (e.g., `.done` files) so `prepare` is idempotent and resumable.

## Logging & Observability

- Human-readable logs by default; `--log-format json` for structured logs. Progress bars for long steps.
- Emit basic counters (genes processed/sec, queue sizes) at intervals.

## I/O Schema

- Outputs: `results/qc_report.jsonl`, `results/qc_summary.csv`, `results/errors.csv`.
- The primary JSONL output will contain a "scorecard" for each gene, including `final_score`, `score_components` (per-pillar scores), and detailed `evidence` blocks.
- Provide small example files and a schema note in the docs.

## Packaging & Release

- CI: lint (`fmt`, `clippy -D warnings`), unit tests, tiny fixture integration tests.
- Release: build artifacts for Linux/macOS (x86_64/aarch64) and Windows; optional Docker image with DIAMOND.

## Preflight Checks

- On startup, verify external tools (`diamond`, optional `mafft`, optional `hmmscan`) and print versions; allow `DIAMOND_BIN`/`HMMER_BIN` overrides.
- Validate write permissions and disk space in output directories.

## Cross-Platform Notes

- Avoid shell-specific features; use `std::process::Command` with explicit args.
- Normalize paths; test on Linux and macOS first, then Windows.

## Risks & Mitigations

- Large artifacts (multi-GB): keep paths configurable; never commit outputs.
- DIAMOND/HMMER version differences: capture command flags in logs; consider `--version` check.
- Bevy API drift: track in `BEVY_GUIDE.md` and pin to 0.17.x until stabilized.
- Network/CDN outages: retry with backoff; mirror URL support.

## Documentation

- Keep `book/` pages aligned: `prepare.md`, `analyze.md` examples, and a “Metrics & Scoring” page.
- Maintain `BEVY_GUIDE.md` with recurring errors and fixes.
- Add a quickstart with sample data and expected outputs for smoke testing.

## Roadmap Update (2025-11-06)

Recent completions
- DIAMOND: pre-run caching, retries; optional chunked mode.
- ECS scheduler scaffold with progress logs.
- Homology coverage metrics (coverage_delta/ratio) + fusion/split flag.
- Intrinsic metrics; MAFFT metrics (optional) in JSONL and CSV.
- Scoring: homology+intrinsic weighted; thresholds; CSV final_score/classification.
- Prepare checkpoints with `--resume` + JSON logs.
- run.json manifest: schema_version, tool versions, file hashes, config snapshot.
- HMMER/Pfam scaffolding and domtblout parser (gated).

Planned next
- Taxonomy pillar (Phase 2): lineage resolution and taxonomy_score (enabled via config/flag); expose in JSONL/CSV; add scoring weight.
- Pfam/HMMER (Phase 2): batch/threaded hmmscan; domains_score; extend domains block and optional CSV columns.
- DIAMOND Auto chunking + per-chunk JSON progress logs; better external-tool diagnostics.
- Observability: analyze start/finish JSON events; optional text progress; sidecar metrics.
- Docs: add “Taxonomy & Domains” and “Scoring” pages; realistic examples from fixtures.
