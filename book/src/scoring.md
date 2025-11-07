% Scoring

AnnoQC combines three pillars into a final score in [0,1]:

- Homology: DIAMOND-based evidence, factoring hit count, bitscore density, and coverage agreement.
- Intrinsic: sequence quality signals like ambiguous fraction, low complexity, and homopolymer length.
- Taxonomy (optional): placeholder presence score (1.0 if top hit resolves to a taxid; else 0.0). Will expand to congruence checks later.

Weights and thresholds
- Configure in `config.example.toml` under `[scoring.weights]` and `[scoring.thresholds]`.
- Example defaults: `homology=0.6`, `intrinsic=0.4`, `taxonomy=0.0`.

Outputs
- JSONL: `score_components` block lists per-pillar scores; `final_score` is the weighted average.
- CSV: includes `final_score` and `classification` based on thresholds; `taxonomy_score` is present when taxonomy is enabled.

Current formulas (as implemented)

- Homology pillar h in [0,1]
  - Top-hit stats are parsed from DIAMOND (`bitscore`, `length`, `qcovhsp`, `scovhsp`, `pident`).
  - bitscore density d = bitscore / max(length, 1) / 5, clamped to [0,1].
  - coverage quality c = (qcovhsp + scovhsp) / 2, clamped to [0,1].
  - coverage penalty p = 1 - min(|qcovhsp - scovhsp|, 1).
  - hit count term n = min(hits/10, 1).
  - h = clamp(0.3·n + 0.3·d + 0.3·c + 0.1·p, 0, 1).

- Intrinsic pillar i in [0,1]
  - i = clamp(0.5·(1 - ambiguous_fraction) + 0.4·(1 - low_complexity_fraction) + 0.1·(1 - max_homopolymer/30), 0, 1).

- Taxonomy pillar t in [0,1]
  - Presence-only placeholder: t = 1.0 if top hit resolves to a taxid; else 0.0.
  - Future: congruence between top hits and expected lineage.

- Domains pillar (JSON only for now)
  - domains_score s_dom = clamp(-log10(top_evalue)/20, 0, 1), derived from hmmscan top Pfam hit.
  - Included in JSONL `domains.domains_score`; not yet part of `final_score`.

- Final score
  - final_score = (w_h·h + w_i·i + w_t·t) / max(w_h + w_i + w_t, ε).
  - Defaults: w_h=0.6, w_i=0.4, w_t=0.0 (taxonomy disabled by default).

Thresholds and classification
- High if final_score ≥ high; Medium if final_score ≥ medium; else Low.
- Defaults: high=0.8, medium=0.5.

Notes & caveats
- Coverage values use DIAMOND `qcovhsp/scovhsp` when available; otherwise qcov is estimated from qstart/qend and query length, and scov may be 0.0 if subject length isn’t available.
- `domains_score` is experimental and excluded from `final_score` pending weighting research.
- Future work: add taxonomy congruence metrics and optional inclusion of `domains_score` in the weighted final score via a new weight.
