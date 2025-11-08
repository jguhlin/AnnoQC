# Next Feature: Orphan Domain Analysis

_Status: implemented — hmmscan now parses HMM coordinates, detects N/C-terminal orphan domains, emits `orphan_status`/`orphan_domain_score`, and wires the new pillar into scoring/docs._

This is the next priority feature to implement.

### Implementation Plan

**1. Objective**

To improve the accuracy of protein "completeness" assessment by detecting partial (i.e., "orphan") protein domains at the N- or C-terminus of a query sequence. The presence of a partial domain is a strong indicator of a truncated gene model.

**2. Background**

Currently, AnnoQC's HMMER integration checks for the *presence* of domains and calculates a `domains_arch_score` based on the set of domains found. It does not evaluate the *completeness* of the individual domain hits themselves. A protein that is missing its true start codon may be annotated with a partial N-terminal domain. This feature will explicitly detect and score this.

**3. Detailed Implementation Steps**

*   **A. Enhance HMMER Output Parsing:**
    *   The current `hmmscan` runs with `--domtblout`. This format contains the necessary information.
    *   The key columns to parse are `ali from`, `ali to` (coordinates on the alignment), `hmm from`, `hmm to` (coordinates on the HMM model), and `hmm length`.
    *   The `hmmer.rs:parse_domtblout` function and the `HmmscanHit` struct need to be extended to parse and store these additional fields.

*   **B. Define "Orphan" Domain Logic:**
    *   For each domain hit, calculate its completeness. A simple metric would be `(hmm_to - hmm_from + 1) / hmm_length`.
    *   Define a domain as an **N-terminal orphan** if:
        1.  It's the first domain in the protein.
        2.  The protein sequence alignment starts far from the beginning of the HMM model (e.g., `hmm from` > 20 amino acids, threshold is tunable).
        3.  The domain hit itself is not a small, repetitive domain that is expected to be partial.
    *   Define a domain as a **C-terminal orphan** if:
        1.  It's the last domain in the protein.
        2.  The protein sequence alignment ends far from the end of the HMM model (e.g., `hmm_length - hmm_to` > 20 amino acids).
    *   Create a new function in `scoring.rs` or `hmmer.rs` that takes an `HmmscanSummary` and returns an `OrphanDomainResult` (e.g., an enum: `None`, `NTerminalOrphan`, `CTerminalOrphan`, `Both`).

*   **C. Integrate into Scoring Model:**
    *   Create a new score component, `orphan_domain_score`. This could be a simple binary score (1.0 if no orphans, 0.0-0.5 if orphans are detected).
    *   Add a new weight, `weights.orphan`, to the `[scoring]` section of the config file to control its influence on the `final_score`.
    *   Update `scoring.rs` to incorporate this new score into the final weighted average.

*   **D. Update Output Files:**
    *   **`qc_report.jsonl`:**
        *   Add a new block, e.g., `"orphan_analysis": { "status": "NTerminalOrphan", "details": "PfamID:PF00069, completeness:0.6" }`.
        *   Add `orphan_domain_score` to the `score_components` object.
    *   **`qc_summary.csv`:**
        *   Add new columns: `orphan_status` (e.g., "None", "N_Orphan") and `orphan_domain_score`.

**4. Acceptance Criteria**

*   A new feature flag or config setting under `[hmmer]` or `[scoring]` enables/disables this analysis.
*   When enabled, the `final_score` for proteins with terminal domain fragments is penalized.
*   The JSONL and CSV outputs are updated with the new fields describing the orphan status and score.
*   Unit tests are added for the orphan detection logic using fixture `domtblout` data representing complete, N-truncated, and C-truncated domain hits.
*   Documentation in the `mdBook` is updated to explain the new analysis and its corresponding outputs.

---

# Next Feature: Advanced Taxonomic Analysis

### Implementation Plan

**1. Objective**

To move beyond single-hit taxonomic validation and instead analyze the entire taxonomic neighborhood of a query's homologs. This will enable robust detection of taxonomic outliers and potential sequence contamination.

**2. Background**

The current taxonomy feature provides a score based on whether the single top DIAMOND hit can be resolved to a taxon. This is a good first step, but it's sensitive to the top hit being an outlier itself. By analyzing the distribution of the top N hits, we can generate a much more stable "consensus" taxonomy and identify queries that don't belong with their neighbors, which is a strong signal of either contamination or interesting biological events like horizontal gene transfer.

**3. Detailed Implementation Steps**

*   **A. Broaden Taxonomic Data Collection:**
    *   In `main.rs`, after parsing the DIAMOND results, don't just look at the top hit for each query. Instead, for each query, collect the top N (e.g., N=20, configurable) `sseqid`s from its `DiamondHitRow` list.
    *   Modify the `TaxonomyResolver` lookup logic to accept a list of accessions and efficiently return a list of `TaxonSummary` objects. This may involve batching lookups for performance.

*   **B. Implement Lowest Common Ancestor (LCA) Algorithm:**
    *   In `taxonomy.rs`, create a new function `find_lca(taxids: &[u32]) -> Option<u32>`.
    *   This function will take a list of NCBI taxonomy IDs. For each ID, it will retrieve its lineage (parent chain up to the root) from the `TaxonomyResolver`.
    *   It will then find the deepest node that is common to all (or a quorum, e.g., >80%) of the lineages. This node is the LCA. The algorithm involves traversing the lineage paths.

*   **C. Develop Contamination & Congruence Scores:**
    *   **Consensus Taxon:** For each query, find the LCA of its top N hits. This is the "consensus taxon" for the query's homologs.
    *   **Congruence Score:** If the query itself has a known taxonomy (e.g., from the input file's metadata, if available), compare its lineage to the consensus taxon's lineage. The score is higher the closer they are. For example, sharing a Family is better than only sharing a Kingdom.
    *   **Contamination Score:** Calculate the fraction of the top N hits that fall outside the consensus taxon's phylum or class. A high score (e.g., >0.9) where most hits are from a completely different domain of life (e.g., a bacterial gene with all eukaryotic hits) is a strong flag for contamination.

*   **D. Integrate into Scoring and Output:**
    *   **Scoring:**
        *   Replace the current simple `taxonomy_score` with the more nuanced `congruence_score`.
        *   The `contamination_score` can be a separate, standalone metric used for flagging rather than for the main quality score.
    *   **`qc_report.jsonl`:**
        *   Modify the `taxonomy` block to include the new information:
            ```json
            "taxonomy": {
              "status": "enabled",
              "top_hit_taxid": 123,
              "consensus_taxid": 456,
              "consensus_name": "Some Consensus Name",
              "consensus_lineage": [...],
              "congruence_score": 0.85,
              "contamination_score": 0.05
            }
            ```
    *   **`qc_summary.csv`:**
        *   Replace `taxonomy_score` with `congruence_score`.
        *   Add new columns: `consensus_taxon`, `contamination_score`.

**4. Acceptance Criteria**

*   A new section in the config file (e.g., `[taxonomy]`) allows configuration of N (number of hits to consider) and thresholds for scoring.
*   The `TaxonomyResolver` is updated with an efficient LCA implementation.
*   The final JSONL and CSV reports contain the new consensus and contamination fields.
*   The overall `final_score` is influenced by the new `congruence_score`.
*   Unit tests are added for the LCA algorithm.
*   The mdBook documentation is updated to explain the new taxonomic analysis, how to interpret the scores, and how it can be used to detect contamination.

---
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

- [x] CLI with `prepare` and `analyze` subcommands.
- [x] Implemented `prepare` (makedb only) and `analyze`.
- [x] External tools wired: DIAMOND, optional MAFFT/HMMER.
- [x] DIAMOND wired: `--version`, `blastp` pre-run, `makedb` in `prepare`.
- [x] DIAMOND TSV parsed for top-hit stats (bitscore, evalue, coverage, density).
- [x] FASTA handling via `needletail` (gz supported).
- [x] FASTA parse for ids/lengths via needletail.
- [ ] `prepare`: makedb → cluster → realign → recluster artifacts.
- [x] Implement makedb + linclust + cluster with `.done` checkpoints; recluster placeholder.
- [x] Add `--resume` support and progress logs for all prepare steps (checkpointed).
- [x] `analyze`: read input → schedule tasks → parse DIAMOND → compute metrics → emit outputs.
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
- [x] Expose threads/approx-id/member-cover; sensible defaults.
- [x] Retry transient DIAMOND errors; concise diagnostics.

## Reproducibility

- [x] Run manifest (tool versions, config snapshot, checksums) at `results/run.json`.
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
