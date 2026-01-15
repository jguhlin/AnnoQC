# Scoring architecture (per-gene / per-protein)

This document describes how `annoqc analyze` turns per-gene evidence into pillar scores and a final score, including optional gates.

For the exact per-pillar formulas, see `book/src/scoring.md` (the “Code map (per-gene scores)” section points to the implementing Rust functions).

## High-level data flow

```mermaid
flowchart TD
    A[Inputs] --> B[DIAMOND blastp TSV]
    A --> C[FASTA sequences]
    A --> D[Optional: reference_fasta]
    A --> E[Optional: hmmscan + pfam_db]
    A --> F[Optional: gff + genome]
    A --> G[Optional: taxonomy resolver inputs]
    A --> H[Optional: plugins / rules]

    B --> I[ECS intake: per-gene DIAMOND stats]
    C --> J[Intrinsic metrics per gene]
    B --> K[Consensus panel selection (aggregated HSPs)]
    K --> L[Length consistency vs panel]

    D --> M{Alignment enabled?}
    M -->|yes| N[MAFFT/SPOA alignment metrics]

    E --> O{HMMER enabled?}
    O -->|yes| P[Query hmmscan summaries]
    O -->|yes| Q[Reference hmmscan summaries (panel)]
    P --> R[Domains architecture score]
    P --> RS[Domains strength score]
    P --> S{Orphan analysis enabled?}
    S -->|yes| T[Orphan domain score]

    G --> U{Taxonomy enabled?}
    U -->|yes| V[Taxonomy evidence + congruence]

    F --> W{Genomic enabled?}
    W -->|yes| X[Genomic splice metrics]

    I --> Y[Homology raw score]
    L --> Z[Homology adjustment (length)]
    N --> AA[Homology adjustment (divergence)]
    Y --> AB[Adjusted homology score]
    Z --> AB
    AA --> AB

    AB --> AC[Weighted mean → raw_score]
    J --> AC
    V --> AC
    R --> AC
    RS --> AC
    T --> AC
    L --> AC
    I --> AC
    N --> AC
    X --> AC

    AC --> AD{Calibration enabled + enough samples?}
    AD -->|yes| AE[base_score (calibrated)]
    AD -->|no| AE[base_score = raw_score]

    AE --> AF[Classification (base_score)]
    H --> AG[Plugin/rule penalties]
    AE --> AH[final_score = clamp(base_score - penalties)]
    AG --> AH
```

## Pipeline stages (where things are computed)

### 1) DIAMOND stats (homology inputs)

- DIAMOND TSV is parsed into:
  - per-gene top-hit summary (`DiamondHitStats`): `src/diamond.rs::parse_tsv_stats`
  - per-gene grouped hit rows for panel selection: `src/diamond.rs::parse_tsv_grouped`
- The `annoqc` scheduler uses Bevy ECS for ingestion/parallelism: `src/ecs.rs::run_scheduler` (invoked from `src/main.rs`).

### 2) Intrinsic metrics (intrinsic pillar input)

- FASTA is read and intrinsic metrics are computed per gene: `src/main.rs::compute_intrinsic_for_ids`
- Per-sequence metrics: `src/metrics.rs::compute_intrinsic`
- The intrinsic *pillar score* currently uses only ambiguity/low-complexity/homopolymer: `src/scoring.rs::compute_intrinsic_score`.

### 3) Consensus panel selection + aggregated HSPs (shared foundation)

- For each gene, DIAMOND hits are aggregated by subject (union-of-HSP spans) to compute `qcov/scov/len_ratio` used for panel selection: `src/consensus.rs::aggregate_hits_by_subject` + `src/consensus.rs::union_span_len`.
- The panel selection phases and dynamic length window live in: `src/consensus.rs::select_panel_with_result`.
- Length consistency is computed against subject lengths from the **selected panel only**:
  - `src/length.rs::compute_length_consistency`
  - `src/main.rs` builds `len_map` during the consensus loop.

### 4) Optional heavy pipelines (alignment + HMMER)

These are executed through the “heavy pipeline” ECS runner which caps concurrency and streams results:
- `src/ecs.rs::run_heavy_pipelines` (called from `src/main.rs`)

**Alignment (MAFFT/SPOA)**
- Gate: requires `reference_fasta` + enough panel hits to form jobs (and a configured aligner backend).
- Alignment metrics (termini, divergence ratio, gap runs, etc.) are computed in: `src/mafft.rs` (see `compute_alignment_metrics`).

**HMMER (hmmscan)**
- Gate: requires `--hmmscan-bin` and a Pfam HMM database (`pfam_db` from args/config).
- Query hmmscan summaries are produced per gene; reference hmmscan summaries are produced for panel sequences.
- Optional Pfam clan collapsing: `src/hmmer.rs::load_pfam_clans` + `src/hmmer.rs::collapse_by_clan`.

### 5) Optional taxonomy

- Gate: taxonomy is enabled only when a `TaxonomyResolver` can be built and the feature isn’t explicitly disabled.
- Evidence + consensus: `src/taxonomy.rs::TaxonomyResolver::summarize_panel`
- Pillar score is the clamped `congruence_score`: `src/scoring.rs::compute_taxonomy_score`.

### 6) Optional genomic context

- Gate: requires `--gff` + `--genome`.
- Per-gene splice-site counts: `src/genomic.rs::analyze_gff_context`
- Genomic pillar score (with optional thresholds): `src/scoring.rs::compute_genomic_score_with_cfg`.

### 7) Pillar scoring, weighting, calibration, penalties

**Per-pillar scoring helpers**
- `src/scoring.rs`: homology raw, intrinsic, taxonomy, subject coverage, divergence score mapping, genomic score.
- `src/hmmer.rs`: domains architecture diagnostics score + orphan analysis.
- `src/length.rs`: length score.
- `src/mafft.rs`: termini + divergence_ratio source.

**Homology adjustments**
- Adjusted homology `h = h_raw × adj_div × adj_len` is applied in `src/main.rs::adjust_homology_score`.
  - Note: divergence/length can be used both as standalone pillars *and* as multipliers on homology when present.

**Weighted mean + presence rules**
- The weighted score + classification are computed in `src/main.rs::build_scores_map`.
- Presence (“is this pillar available for this gene?”) is computed per pillar; missing pillars are generally excluded from the denominator.
- Special-case: **homology**, **domains architecture**, and **domains strength** still count in the denominator when missing (penalizes “no evidence” genes).

**Calibration**
- Gate: `calibration.mode != Off` AND the run has at least `min_samples` entries and `min_unique` unique raw scores.
- Implementations: `src/main.rs::{apply_percentile_calibration, apply_isotonic_calibration, calibration_has_min_samples}`.

**Plugins / rules**
- Gate: only if any Extism plugins and/or Rhai rules are configured.
- Plugins run via `src/plugins.rs::run_plugin`; rules via `src/rhai_rules.rs::RhaiRuntime::run`.
- The summed penalty is applied per gene in `src/main.rs::render_gene_record`:
  - `final_score = clamp(base_score - plugin_penalty)`
  - Classification is currently based on the base score (pre-penalty).

## Optional gates (summary)

| Feature | Gate (must be true) | Outputs / pillars affected |
| --- | --- | --- |
| Taxonomy | resolver builds successfully + not disabled | taxonomy pillar; taxonomy warnings/labels |
| HMMER | `hmmscan_bin` + `pfam_db` | domains architecture pillar; domains strength pillar; orphan (if enabled) |
| Orphan analysis | HMMER enabled + `hmmer.orphan_analysis=true` + not `--disable-orphan-analysis` | orphan pillar |
| Alignment | `reference_fasta` present + enough panel hits to form jobs + aligner configured | conserved_regions pillar; termini pillar; divergence pillar; homology adjustment; many diagnostics |
| Genomic | `--gff` + `--genome` | genomic pillar |
| Calibration | mode != Off + enough samples/unique scores | final score remap (base score) |
| Plugins / rules | any configured | final score penalty (post-calibration) |

## Known gaps

Some metrics are computed and exported but are not yet first-class scoring pillars (or some config sections are not wired). See `NEXT_STEPS_TODO.md` for the current checklist and code touchpoints.
