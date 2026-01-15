# Next Steps (Quality/Scoring Pipeline)

This file tracks known gaps in the current per-gene / per-protein QC scoring pipeline and related outputs.

## Scoring inputs computed but not scored (diagnostics-only)

- [ ] ORF/start/stop heuristics: fold `IntrinsicMetrics.{orf_start_score,internal_stop_count,terminal_stop,start_methionine,alt_start_pos}` into `compute_intrinsic_score` or a new pillar (source: `src/metrics.rs::compute_intrinsic`; scoring: `src/scoring.rs::compute_intrinsic_score`).
- [ ] Alignment QC beyond termini/divergence: expand/tune the `conserved_regions` pillar so it fully captures GV-style “missing/extra conserved blocks” (source metrics: `src/mafft.rs::compute_alignment_metrics`; scoring today: `src/scoring.rs::compute_conserved_regions_score` + `src/main.rs::build_scores_map`; warnings still emitted in `src/main.rs::render_gene_record`).
- [ ] Structural-variation tuning: refine `structvar_multiplier` mapping and/or add config overrides for thresholds and penalty strength (analysis: `src/structvar.rs::analyze`; penalty applied in `src/main.rs::build_scores_map`; multiplier mapping: `src/scoring.rs::compute_structvar_multiplier`).
- [ ] Domains: implement and weight an order/architecture penalty (currently `w_ord = 0.0` / `order_pen = 0.0`) (source: `src/hmmer.rs::domains_architecture_diagnostics`).
- [ ] Domains strength (`domains_score` based on top e-value) is diagnostic-only; decide whether to incorporate it into final score (computed in `src/main.rs` when building the JSONL `domains` block).

## Config sections that exist but are not wired into scoring

- [ ] Wire `[scoring.homology]` (e.g., `qcov_min, scov_min, density_min, evalue_max`) into either:
  - homology score computation (`src/scoring.rs::compute_homology_score`), and/or
  - panel selection thresholds (`src/consensus.rs::ConsensusConfig` / `select_panel_with_result`), and/or
  - diagnostics/warnings.
  (`config.example.toml` has this section but `src/main.rs::ScoringConfigOverride` does not define it.)
- [ ] Wire `[scoring.intrinsic]` thresholds (e.g., `ambiguous_max, max_homopolymer`) into scoring/warnings (same mismatch as above).
- [ ] Consider adding explicit config for termini/divergence scoring cutoffs (currently hard-coded in `src/mafft.rs` and `src/scoring.rs::compute_divergence_score`).

## Output/schema gaps (harder to audit than it should be)

- [ ] JSONL `score_components` is missing divergence (only `alignment.divergence_ratio` is emitted); add `divergence_score` to `score_components` for symmetry with CSV (JSON built in `src/main.rs::render_gene_record`).
- [ ] Classification is computed pre-plugin-penalty; decide whether to:
  - recompute classification after penalties, or
  - emit both `classification_base` and `classification_final`, or
  - explicitly annotate that classification is “pre-penalty” (classification today: `src/main.rs::build_scores_map`; penalties applied: `src/main.rs::render_gene_record`).
- [ ] Consider emitting a more explicit score breakdown in JSONL:
  - `raw_score` (weighted mean)
  - `base_score` (after calibration)
  - `final_score` (after plugins/rules)

## Implementation notes / touchpoints

- Final weighted score + presence rules live in `src/main.rs::build_scores_map` (special-case: missing homology/domains still count in denominator).
- Pillar score helpers live in `src/scoring.rs` (homology/intrinsic/taxonomy/subject_cov/divergence/genomic).
- Docs mapping is in `book/src/scoring.md` under “Code map (per-gene scores)”; update it as these TODOs get implemented.
