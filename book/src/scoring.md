% Scoring

AnnoQC combines pillar scores in [0,1] into a weighted average. Pillars:

- Homology (h): DIAMOND-based evidence (hit count, bitscore density, coverage agreement).
- Intrinsic (i): sequence quality (ambiguous, low-complexity, homopolymer; ORF heuristics).
- Taxonomy (t, optional): consensus-based lineage congruence (LCA across top hits, penalized by contamination/outlier rate).
- Domains Architecture (d, optional): Pfam clan-collapsed architecture agreement vs homolog panel.
- Length Consistency (l, optional): query length vs homolog panel median using robust z-scores.
- Orphan Domain Integrity (o, optional): detects N- or C-terminal Pfam domains that align far from the model edges, a strong hint of truncated gene models.

Weights and thresholds
- Configure in `config.example.toml` under `[scoring.weights]` and `[scoring.thresholds]`.
- Example defaults: `homology=0.6`, `intrinsic=0.4`, others `0.0`.

Outputs
- JSONL: `score_components` lists per-pillar scores; `final_score` is the weighted average.
- CSV (`--csv-verbose`): includes homology_score, intrinsic_score, taxonomy_score, domains_score, domains_arch_score, orphan_domain_score, length_score, length_ratio, length_class, start_concordance, start_class (compact mode keeps `domains_score`, `domains_arch_score`, `orphan_domain_score`, and `orphan_status`).

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
  - Inputs: top `[taxonomy].top_hits` DIAMOND subjects (default 20) after redundancy trimming.
  - Require at least `[taxonomy].min_consensus` resolved hits (default 5). Collect their lineage IDs and find the deepest node supported by ≥ `[taxonomy].min_support` fraction (default 0.75).
  - Let `support_fraction = support_hits / considered_hits` and `depth_norm = depth / max_depth`. The congruence score is `t = clamp(support_fraction · depth_norm, 0, 1)`. Contamination is surfaced separately as `1 − support_fraction`.
  - JSONL/CSV now emit `congruence_score`, `contamination_score`, support counts, and a consensus taxon label. The weighted scoring pipeline uses the congruence score whenever `[scoring.weights].taxonomy > 0`.

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

- Orphan domain integrity o
  - Enabled when `hmmscan` runs (`--hmmscan-bin …`) and `[hmmer].orphan_analysis = true` (override with `--disable-orphan-analysis`).
  - Collapse Pfam hits (clan-aware), sort domains by alignment start, and inspect the first and last domains only.
  - Flag the N-terminal domain if it begins >20 HMM residues away from the model start (and the HMM length ≥60 aa). Flag the C-terminal domain if it ends >20 residues shy of the model end.
  - o = 1.0 when no domain is flagged, 0.5 when only one terminus looks truncated, and 0.0 when both N and C appear truncated. JSON/CSV expose `orphan_status` along with the Pfam accession(s) that triggered the warning.

Start-concordance (gated)
- Alignment-based assessment of whether the query begins at the consensus N-terminus of its homologs.
- Requires a homolog panel of at least `min_hits` (default 5). If panel_size < `min_hits`, MAFFT is skipped and `start_concordance`/`start_class` are omitted.
- Computation (when enabled):
  - Align query + panel with MAFFT, using `mafft_threads_per_job` threads per alignment (default 4 or 8 based on `--threads`) while pooling up to `mafft_max_jobs` concurrent MAFFT invocations.
  - Find first non-gap column per sequence; use the panel median as consensus start.
  - start_concordance = clamp(1 − |query_start − consensus_start|/30, 0, 1).
  - start_class: LikelyComplete (|Δ|≤3), LikelyNTruncated (Δ>3), LikelyNExtended (Δ<−3).

- Final score
  - final_score = (w_h·h + w_i·i + w_t·t + w_d·d + w_l·l + w_o·o) / max(w_h + w_i + w_t + w_d + w_l + w_o, ε)

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
