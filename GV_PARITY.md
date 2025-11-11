# GeneValidator Parity — TODO Checklist

Actionable items to reach (or surpass) GeneValidator parity. Boxes indicate status.

Core parity tasks
- [x] Advanced structural variation detection (fusion/split/duplication)
  - Group DIAMOND HSPs per subject; classify tiling/overlap patterns; emit explicit JSONL warnings and CSV columns.
  - Status: completed — per-subject coverage/orientation summaries feed JSON/CSV, fusion pairs must agree in orientation/order, and MAFFT warnings surface downstream via the shared warnings column.
  - Inputs: BLAST6 with HSP coordinates; thresholds for minimal span and gap.
- [ ] Nucleotide input + ORF validation
  - Add nt mode: FASTA(nt) → ORF calling → AA; run same pipeline.
  - Validate ORF vs homolog panel: frameshifts, early stops, retained introns.
  - Status: paused (we added a simple placeholder to enable nt mode; robust ORF finder and validations will follow later).
- [x] Detailed alignment-based warnings
  - From MAFFT, detect missing exons (large query gaps) and retained introns (large insertions); flag motif disruptions.
  - Status: implemented — the MAFFT block now records `missing_exon_run`/`retained_intron_run`, and warnings (`MissingExonPossible`, `RetainedIntronPossible`) flow through JSON/CSV.

Implemented / partial
- [x] DIAMOND coverage delta as fusion/split hint (exposes `coverage_delta`, fusion flag via threshold).
- [x] Consensus homolog selection (min_hits=5, max_panel=20; phased filters; redundancy reduction).
- [x] Domains architecture scoring (Pfam clan-collapsed) with diagnostics; CSV/JSONL outputs.
- [x] Length distribution pillar (median/MAD z; ratio; class) weighted via `[scoring.weights].length`.
- [x] Start-concordance via MAFFT, gated by panel size ≥ `min_hits`.
- [x] CSV verbose per-pillar; JSONL `score_components`.
- [x] Default NoData classification when no DIAMOND hits and no Pfam domains.
- [x] Dump matched reference sequences (`--dump-matches-best/--dump-matches-gene`).

Planned enhancements (GV-inspired)
- [x] Improved panel selection heuristics (length/score distributions; maintain ≥5 informative homologs via coverage/identity phases, length-ratio windows, and redundancy trimming; emits panel_debug.csv diagnostics).
- [x] Taxonomy pillar 2.0 (consensus hits + lineage LCA, congruence/contamination scores, CSV/JSONL outputs, weighted in scoring).
- [ ] Per-step observability (ECS inflight/queue traces; phase timing; optional TTY progress bar when not JSON).
