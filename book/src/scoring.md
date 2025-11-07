% Scoring

AnnoQC combines pillar scores in [0,1] into a weighted average. Pillars:

- Homology (h): DIAMOND-based evidence (hit count, bitscore density, coverage agreement).
- Intrinsic (i): sequence quality (ambiguous, low-complexity, homopolymer; ORF heuristics).
- Taxonomy (t, optional): presence-only for now (1.0 if taxid resolved).
- Domains Architecture (d, optional): Pfam clan-collapsed architecture agreement vs homolog panel.
- Length Consistency (l, optional): query length vs homolog panel median using robust z-scores.

Weights and thresholds
- Configure in `config.example.toml` under `[scoring.weights]` and `[scoring.thresholds]`.
- Example defaults: `homology=0.6`, `intrinsic=0.4`, others `0.0`.

Outputs
- JSONL: `score_components` lists per-pillar scores; `final_score` is the weighted average.
- CSV (`--csv-verbose`): includes homology_score, intrinsic_score, taxonomy_score, domains_score, domains_arch_score, length_score, length_ratio, length_class, start_concordance, start_class.

Formulas (implemented)

- Homology h
  - Inputs: DIAMOND `bitscore`, `length`, `qcovhsp`, `scovhsp`, `pident`.
  - bitscore density d = clamp(bitscore / max(length, 1) / 5, 0, 1).
  - coverage quality c = (qcovhsp + scovhsp)/2; penalty p = 1 − min(|qcovhsp − scovhsp|, 1).
  - hits term n = min(hits/10, 1).
  - h = clamp(0.3·n + 0.3·d + 0.3·c + 0.1·p, 0, 1).

- Intrinsic i
  - i = clamp(0.5·(1 − ambiguous_fraction) + 0.4·(1 − low_complexity_fraction) + 0.1·(1 − max_homopolymer/30), 0, 1).

- Taxonomy t
  - Presence-only placeholder: t = 1.0 if top hit maps to a taxid; else 0.0.

- Domains architecture d (Pfam clan-collapsed)
  - Collapse accessions→clans; for each clan, compute frequency across homolog panel.
  - Core if freq≥0.7; Accessory if freq≥0.3.
  - d = clamp(0.6·recall_core + 0.3·precision_acc − 0.1·extras_pen, 0, 1).

- Length consistency l
  - Let Lq be query length; Ls subject lengths from the consensus panel.
  - median m = median(Ls), robust spread MADn = 1.4826·MAD(Ls).
  - z = (Lq − m) / max(MADn, 1.0).
  - l = exp(−|z|/2). Also report length_ratio = Lq/m and length_class:
    - LikelyNTruncated if ratio < 0.8; LikelyNExtended if ratio > 1.2; else InRange.

- Final score
  - final_score = (w_h·h + w_i·i + w_t·t + w_d·d + w_l·l) / max(w_h + w_i + w_t + w_d + w_l, ε)

Thresholds and classification
- High if final_score ≥ high; Medium if final_score ≥ medium; else Low. Defaults: high=0.8, medium=0.5.

Worked example (a9_head200)

Config weights: `homology=0.5, intrinsic=0.3, domains=0.1, length=0.1, taxonomy=0.0`.

```
final_score ≈ 0.86
homology ≈ 0.78, intrinsic ≈ 0.99, domains_arch_score ≈ 0.65, length_score ≈ 0.99
length_ratio ≈ 1.00 (InRange), start_concordance ≈ 1.0 (LikelyComplete)
```

Notes & caveats
- Coverage uses DIAMOND `qcovhsp/scovhsp` when available; we request `qlen/slen` to enable robust length metrics.
