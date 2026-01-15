# Future TODOs (Planned Features / Pipeline Extensions)

This file tracks larger “next wave” ideas that require new inputs or new data products (e.g. RNA-seq), and therefore are not immediately actionable in the current codebase.

## RNA-seq evidence pillar (run-provided, per-gene variable)

Goal: if RNA-seq is provided for the *run* (genome-wide), then *each gene* should be mildly penalized when it lacks RNA-seq support, even if other pillars look good. If RNA-seq is not provided (proteins-only runs), there should be no penalty.

- [ ] Define the supported RNA-seq input artifact(s) and ingestion path:
  - Option A: quant table (Salmon/Kallisto) keyed by transcript/gene: TPM/NumReads, effective length, etc.
  - Option B: junction support + coverage summaries from BAM/CRAM (per intron, per exon, per transcript).
  - Option C: precomputed per-gene “expression support” JSON/TSV produced by an external pipeline.
- [ ] Add a new per-gene metrics struct (e.g. `RnaseqMetrics`) and compute it during `analyze`.
- [ ] Add a new scoring helper (e.g. `src/scoring.rs::compute_rnaseq_score`) mapping to `[0,1]`.
- [ ] Add a new weight key (e.g. `scoring.weights.rnaseq`) and presence semantics:
  - “Run has RNA-seq” (global): include `w_rnaseq` in denominator for **all genes**.
  - “Gene has RNA-seq support” (per-gene): contribute `w_rnaseq * rnaseq_score` to numerator; missing support should imply a low score (not “N/A”).
- [ ] Emit RNA-seq fields in outputs:
  - JSONL: `rnaseq` block + `score_components.rnaseq_score`
  - CSV: `rnaseq_score` + key diagnostics (TPM/support/junction counts)

## Split the “weighted mean” denominator (required vs optional)

Goal: make the scoring math and reporting clearer by explicitly separating:

- `expected_denom`: weights for pillars that are expected given the run inputs (e.g. genome+gff present ⇒ genomic expected; rnaseq present ⇒ rnaseq expected; hmmer enabled ⇒ domains expected).
- `optional_denom`: weights for pillars that can legitimately be “not applicable” per-gene even when the run has the input (e.g. termini requires enough aligned sequences).

This would generalize the current special-case denominator rule in `src/main.rs::build_scores_map` (today: missing homology/domains still count in the denominator) into a consistent “run-provided vs gene-provided” policy.

