% Taxonomy & Domains

This page outlines how AnnoQC enriches results with taxonomy lineage and protein domains.

- Taxonomy inputs
  - Cache (tab-separated `accession\ttaxid`) via `--taxonomy-cache` or `[taxonomy].cache_path`.
  - Reference FASTA with UniProt headers to build a cache on-the-fly when the cache is missing.
  - NCBI taxdump (`nodes.dmp`, `names.dmp`) via `--taxonomy-taxdump-dir` or `[taxonomy].taxdump_dir` to reconstruct scientific names and lineage.

- Resolution
  - For each gene we take up to `[taxonomy].top_hits` (default 20) DIAMOND subjects, canonicalize their accessions, deduplicate, and look them up via the cache.
  - The resulting taxids feed a Lowest Common Ancestor pass that demands at least `[taxonomy].min_consensus` resolved hits (default 5) and a quorum of `[taxonomy].min_support` (default 0.75).
  - The deepest node that satisfies the quorum becomes the consensus taxon; we record support counts, support fraction, and depth to approximate how specific the consensus is.
  - If no node meets the strict quorum but the panel still contains hits, we fall back to a coarse consensus at `[taxonomy].coarse_rank_index` (default index 1 ⇒ superkingdom). This is useful for “bacterial vs eukaryote” checks when species-level agreement is sparse.
  - Transcripts that share a root ID like `geneX.t1/geneX.t2` borrow consensus evidence across isoforms (`detail = Borrowed`) so alternative transcripts inherit the lineage call of their best-supported sibling.
  - Even when the consensus cannot be established (too few hits, or no resolver) we still emit a `detail` string so downstream dashboards can distinguish "NoHits" from "NoResolver" cases.

- JSONL output
  - `taxonomy.status`: `enabled` or `disabled` plus a `detail` reason (`NoHits`, `InsufficientHits`, `Consensus`, `NoResolver`).
  - When a consensus exists the block includes `top_hit` (taxid/name/lineage), `consensus_taxid`, `consensus_name`, `consensus_lineage`, `consensus_depth`, `consensus_rank`, `support_hits`, `considered_hits`, `support_fraction`, and both `congruence_score` + `contamination_score`.
  - The taxonomy pillar in `score_components` now mirrors the congruence score so weighting `[scoring.weights].taxonomy` actually reflects lineage agreement.

- CSV output
  - Compact mode now includes `taxonomy_score` (the congruence score), `taxonomy_contamination`, `taxonomy_support`, `taxonomy_considered`, `consensus_taxon`, and `taxonomy_status` (detail string).
  - `--csv-verbose` keeps the per-pillar `taxonomy_score` column near the front and appends the same contamination/support/consensus columns near the end for filtering.

- Domains (hmmscan)
  - When `--hmmscan-bin` and `pfam_db` are provided, JSONL includes a `domains` block with `hits_count`, `top_accession`, `top_evalue`, and per-hit details (`target_name`, `accession`, `evalue`, `score`, `bias`).
  - `domains_score` = clamp(-log10(top_evalue)/20, 0, 1) (fast proxy for domain strength).
  - Clan collapse (recommended): provide `pfam_clans` TSV to collapse Pfam accessions to clans for cross-gene architecture comparisons.
  - Architecture score (`domains_arch_score`): compare the query’s clan set to the consensus homolog panel.
    - Compute clan frequency across homologs; Core if ≥0.7; Accessory if ≥0.3.
    - Score = clamp(0.6·recall_core + 0.3·precision_acc − 0.1·extras_pen, 0, 1).
  - CSV (`--csv-verbose`) exposes `domains_score` and `domains_arch_score`.
  - Flags/config:
    - `--hmmer-top-n` or `[hmmer].top_n` to cap per-gene hits in JSONL (default 5).
    - `--hmmer-threads` or `[hmmer].threads` to override thread count for hmmscan.
