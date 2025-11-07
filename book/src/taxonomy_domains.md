% Taxonomy & Domains

This page outlines how AnnoQC enriches results with taxonomy lineage and protein domains.

- Taxonomy inputs
  - Cache (tab-separated `accession\ttaxid`) via `--taxonomy-cache` or `[taxonomy].cache_path`.
  - Reference FASTA with UniProt headers to build a cache on-the-fly when the cache is missing.
  - NCBI taxdump (`nodes.dmp`, `names.dmp`) via `--taxonomy-taxdump-dir` or `[taxonomy].taxdump_dir` to reconstruct scientific names and lineage.

- Resolution
  - For each gene, the DIAMOND top hit is canonicalized (e.g., `sp|P12345|...` → `P12345`) and looked up in the cache.
  - If a taxid is found, lineage is reconstructed using taxdump; otherwise only the taxid is surfaced if available.

- JSONL output
  - `taxonomy.status`: `enabled` or `disabled`.
  - When enabled and resolved, emits `taxonomy.taxid`, `taxonomy.name`, and `taxonomy.lineage` (root → leaf).
  - Score components include a placeholder `taxonomy` score (1.0 if resolved, 0.0 otherwise) weighted by `[scoring.weights].taxonomy`.

- CSV output
  - Adds `taxonomy_score` and `taxonomy_status` columns; currently presence/absence based.

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
