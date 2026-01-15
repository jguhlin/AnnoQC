# Analyze

AnnoQC’s `analyze` subcommand combines DIAMOND homology searches, MAFFT conserved-region alignments, intrinsic sequence checks, and optional taxonomy resolution to produce a per-gene scorecard. The JSONL output exposes detailed evidence, while the CSV mirrors the headline metrics for downstream dashboards.

## Prerequisites

- DIAMOND installed and the SwissProt database prepared via `cargo run -- prepare`.
- (Optional) MAFFT installed when you want conserved-region analysis (`--mafft-bin`).
- A reference FASTA that contains the DIAMOND subjects (defaults to `uniprot_sprot.fasta.gz`).

## Example Run

1. Copy and adjust `config.example.toml`:

```toml
fasta = "examples/demo.faa"
db = "uniprot_sprot.dmnd"
out = "results"
top = 5
diamond_bin = "diamond"
mafft_bin = "/usr/bin/mafft"
reference_fasta = "uniprot_sprot.fasta.gz"
[scoring.weights]
homology = 0.6
intrinsic = 0.4
```

2. Run the analyzer:

```bash
cargo run -- analyze --config config.example.toml
```

This produces `results/qc_report.jsonl`, `results/qc_summary.csv`, and `results/run.json` (manifest with schema_version, tool versions, file hashes, and a config snapshot). The JSONL and CSV files now begin with a small metadata preamble so downstream tooling can quickly detect schema/tool changes: the JSONL's first line is a `{"type":"metadata",...}` record, and the CSV's first line is a comment (prefixed with `#`) that names the schema version, the AnnoQC version, and the detected DIAMOND/MAFFT/HMMER versions.

If you only want specific outputs, pass `--report-format jsonl|csv|parquet|all` (default: `all`). JSONL/CSV can be appended with `--resume`, which skips genes already rendered in the existing output files.

To print the resolved scoring rubric without running DIAMOND, use `--dry-run`:

```bash
cargo run -- analyze --config config.example.toml --dry-run
```

To inspect a single gene after a run, use:

```bash
cargo run -- explain Gene00042 --out results
```

The explain output prints the per-pillar scores, key homology/taxonomy highlights, warnings, and any plugin penalties.

## JSONL Scorecard Snapshot

Each JSON line records the final score, per-pillar components, evidence blocks, and warnings. Fields include homology coverage metrics and a fusion/split flag, intrinsic metrics, and (optionally) an alignment block when MAFFT is enabled. The `panel_provenance` block reports how many consensus panel IDs came from SwissProt vs reference proteomes vs cluster backfill (IDs added from `clusters.recluster` when the panel is scarce), so you can immediately see when a gene’s downstream optional pillars relied on backfill/rescue.

```json
{
  "gene_id": "Gene00042",
  "taxonomy": {
    "status": "enabled",
    "detail": "Consensus",
    "top_hit": { "taxid": 562, "name": "Escherichia coli", "lineage": ["root", "Bacteria", "Gammaproteobacteria", "Escherichia coli"] },
    "consensus_taxid": 561,
    "consensus_name": "Escherichia",
    "consensus_lineage": ["root", "Bacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia"],
    "support_hits": 17,
    "considered_hits": 20,
    "support_fraction": 0.85,
    "congruence_score": 0.83,
    "contamination_score": 0.17
  },
  "final_score": 0.87,
  "score_components": { "homology": 0.92, "intrinsic": 0.81, "taxonomy": 0.83 },
  "homology": {
    "hits_count": 12,
    "top_hit": "sp|P12345|REF_HUMAN",
    "top_bitscore": 212.6,
    "top_evalue": 3.2e-68,
    "top_qcov": 0.94,
    "top_scov": 0.91,
    "bitscore_density": 3.21,
    "coverage_ratio": 1.03,
    "coverage_delta": 0.03,
    "fusion_split_flag": false
  },
  "intrinsic": {
    "ambiguous_fraction": 0.0,
    "max_homopolymer": 5,
    "low_complexity_fraction": 0.08,
    "low_complexity_windows": 2,
    "orf_has_start": true,
    "orf_has_stop": true,
    "orf_internal_stop": false
  },
  "alignment": {
    "mafft_enabled": true,
    "strategy_used": "Add",
    "conserved_fraction": 0.89,
    "pairwise_identity": 0.86,
    "sequences_aligned": 6,
    "top_hits_considered": 5,
    "query_gap_fraction": 0.04,
    "gap_run_count": 1,
    "max_gap_run": 4,
    "motif_mismatch_fraction": 0.11
  },
  "warnings": []
}
```

Rows with potential issues advertise warnings such as missing start/stop codons, internal stops, excessive low-complexity content, or suspicious fusion/split coverage. Alignment-aware signals (`MissingExonPossible`, `RetainedIntronPossible`) and structural hints (`FusionPossible`, `SplitPossible`) are now emitted whenever the MAFFT panel or DIAMOND tiling heuristics detect the corresponding anomalies.

The `structvar` block enriches fusion/split diagnostics with a `subjects[]` array listing per-subject query and subject coverage, dominant strand, and order conflicts so you can quickly inspect which references are tiling the query when a warning fires.

## CSV Columns

`qc_summary.csv` contains headline homology numbers, coverage metrics, final_score, optional classification, and (when MAFFT is enabled) alignment summary metrics. The file starts with a metadata comment line before the header row, so make sure your CSV parser tolerates commented lines (anything beginning with `#`). Columns `panel_swissprot`, `panel_refprot`, and `panel_cluster` track how many homologs from each source made it into the final panel, mirroring the JSON `panel_provenance` block. A `warnings` column eases quick filtering:

```
gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,domains_score,domains_arch_score,orphan_domain_score,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_score,taxonomy_contamination,taxonomy_support,taxonomy_considered,consensus_taxon,taxonomy_status,warnings
Gene00042,12,sp|P12345|REF_HUMAN,212.600,3.20e-68,0.940,0.910,3.210,0.030,1.030,0,0.870,High,1,0.890,0.860,6,0.040,1,4,0.980,0.640,1.000,,,,,,,LikelyComplete,0.8300,0.1700,17,20,"Escherichia (561)",Consensus,
```

The CSV is suitable for spreadsheets or dashboards, while the JSONL is richer for downstream pipelines.

## Useful Options

- HMMER/Pfam
  - `--hmmscan-bin /path/to/hmmscan` and `pfam_db` in config enable domain summaries.
  - `--hmmer-top-n N` or `[hmmer].top_n` (default 5) caps per-gene domain hits in JSONL.
  - `--hmmer-threads N` or `[hmmer].threads` overrides threads used for hmmscan; otherwise uses global `--threads`.
  - Disable the orphan-domain integrity pillar with `--disable-orphan-analysis` (or set `[hmmer].orphan_analysis = false`).

- Taxonomy
  - `--enable-taxonomy` flips the pillar on (or set `[taxonomy].enabled = true`).
  - `--taxonomy-top-hits N` overrides `[taxonomy].top_hits` (default 20) when selecting subjects for the LCA.
  - `--taxonomy-min-consensus N` enforces how many resolved hits are required before we attempt a full-depth consensus (default 5).
  - `--taxonomy-min-support F` overrides `[taxonomy].min_support` (default 0.75) for the quorum fraction used in the LCA walk.
  - `--taxonomy-coarse-rank-index K` tells the fallback to evaluate lineage level `K` (default 1 ⇒ superkingdom) when strict consensus fails.
  - `--taxonomy-coarse-min-support F` loosens or tightens the fallback quorum (default 0.6).

- Alignment warnings
  - `--alignment-missing-exon N` (or `alignment_missing_exon = N` in your config) sets the run length that triggers `MissingExonPossible` (default 30 columns where only the query is gapped).
  - `--alignment-retained-intron N` / `alignment_retained_intron = N` controls the retained-intron warning (default 30 columns where the query is filled but ≥70% of homologs are gapped).

- Score calibration
  - `--calibration-mode Percentile|Isotonic` (or `[calibration].mode = "Percentile|Isotonic"`) remaps the weighted final scores within the current run. The default `Off` mode emits the raw weighted scores.
  - Calibration is skipped if the run has too few genes or unique scores; override with `[calibration].min_samples` / `[calibration].min_unique`.

- Alignment backend
  - `--aligner spoa` uses the built-in SIMD Partial Order Alignment engine. This is significantly faster than spawning external MAFFT processes and is recommended for high-throughput runs.
  - `--aligner mafft` (default) uses the external `mafft` binary. This requires `mafft` to be in your PATH or specified via `--mafft-bin`.

- MAFFT pooling
  - `--mafft-threads-per-job N` or `mafft_threads_per_job = N` chooses how many MAFFT worker threads each alignment receives (default heuristic: 4 or 8 depending on `--threads`).
  - `--mafft-max-jobs M` or `mafft_max_jobs = M` caps concurrent MAFFT alignments; the default is `floor((threads−1)/mafft_threads_per_job)` so the scheduler always leaves at least one CPU for coordination.
  - `--render-max-jobs M` (or `render_max_jobs = M`) controls how many streaming render tasks write JSONL/CSV concurrently; the default matches `--threads`, but you can lower it if disk bandwidth is saturated or raise it when emit becomes the bottleneck.
  - Alignment metrics in JSON/CSV now include both `start_concordance/start_class` and `end_concordance/end_class`. These are only evaluated when ≥5 homologs are aligned; otherwise the fields are blank/`InsufficientPanel`.

- StructVar diagnostics
  - The JSON `structvar` block includes a `subjects[]` array with per-subject coverage, dominant strand, and ordering information; CSV warnings echo the major structvar states (`FusionPossible`, `SplitPossible`).

- DIAMOND mode
  - `--diamond-mode auto|single|batch` (default `auto`).
  - Auto uses Single up to `--diamond-auto-threshold` or `[diamond].auto_threshold` (default 200,000 queries), otherwise Batch.
  - In Batch and `--log-format json`, per-chunk progress events are logged.

- High-score exports
  - Use `--export-high` (or `export_high = true`) to dump every gene classified as `High` into `<out>/high_scoring.faa`. Override the destination with `--export-high-path path` or `export_high_path = "..."` in your config.

- Reports & resume
  - `--report-format jsonl|csv|parquet|all` controls which outputs are emitted (default `all`).
  - `--resume` appends to existing JSONL/CSV and skips already-rendered genes; parquet append is not supported (parquet output is skipped with a warning when resuming).

- Explain
  - `annoqc explain <gene_id> --out <results_dir>` prints a per-gene breakdown from `qc_report.jsonl`.

- Rules (Rhai)
  - `--rhai path/to/rule.rhai` loads a Rhai script and runs `fn analyze(input)` for each gene.
  - Return either a map with `name/score/penalty/metadata` or a numeric score. Rules run alongside WASM plugins and their penalties are subtracted from `final_score`.

## Release packaging & Docker images

- Run `scripts/package_release.sh` to build a release tarball under `dist/`. The script auto-detects the host triple, but you can cross-compile via `--target x86_64-unknown-linux-gnu` (or any Rust target that your toolchain supports).
- Add `--docker` (or `--docker-tag registry/annoqc:0.1.1`) to produce a container image from the repo `Dockerfile` while still emitting the tarball. This keeps release artifacts consistent across environments.
- All packages bundle the `annoqc` binary plus `LICENSE`, `config.example.toml`, and `config.weights_demo.toml`, so users can unpack-and-run immediately.

## Consensus Panels & Diagnostics

Before MAFFT, domain comparisons, and length scoring kick in, AnnoQC builds a consensus panel for each query from its ranked DIAMOND hits. The selector enforces:

- Minimum informative hits (`[consensus].min_hits`, default 5) and a maximum panel size (`max_panel`, default 20).
- Phase-based filters on query/subject coverage, e-value, and identity. The primary phase honors `[consensus].filt_qcov`, `filt_scov`, `filt_evalue`, and `filt_pident`; later phases relax toward GeneValidator-like heuristics when the panel would otherwise be too small.
- Length-ratio windows anchored to the query. `[consensus].len_ratio_tolerance` (default ±30%) defines the tight window; later phases widen to 0.5–1.5× and 0.4–2.0× if more homology evidence is needed.
- Redundancy trimming. Hits with DIAMOND `pident` ≥ `[consensus].redundancy_pident` (default 90%) are capped by `max_high_identity` (default 3) so near-identical isoforms do not crowd out diverse evidence.
- Taxonomy-aware diversity. `[consensus].diversity_rank_index` chooses which lineage level to treat as a "clade" (default index 3 ≈ class), and `[consensus].diversity_rank_cap` limits how many homologs from the same clade may enter the panel (default 3). When hard caps would leave the panel below `min_hits`, we fall back and admit the deferred hits, and the effect is summarized in `panel_debug.csv` along with the number of clades represented.

Each run now emits `panel_debug.csv` in the results directory. Columns capture how many hits passed the filters, which phase produced the selected panel, the observed length-ratio span, median % identity, and how many high-identity hits were dropped. Use this file to tune cutoffs on new datasets (e.g., loosen coverage filters for fragmented assemblies or tighten identity caps for clonal panels).

## Packaging Releases

Ship ready-to-run binaries with the helper script:

```
pixi run bash scripts/package_release.sh
```

The script builds `target/release/AnnoQC`, stages it alongside `LICENSE` and the sample configs, and emits a compressed tarball plus `.sha256` checksum under `dist/` (e.g., `annoqc-0.1.1-x86_64-unknown-linux-gnu.tar.gz`). Run it on each target platform you need binaries for; the archive name encodes both the version (from `Cargo.toml`) and the local Rust target triple.
- RefProt fallback
  - `--refprot-trigger-k` (or `[refprot].trigger_k`) controls how many SwissProt members must be selected before reference-proteome hits are considered (default = `[consensus].min_hits`).
  - Tighten filters with `--refprot-min-qcov`, `--refprot-min-scov`, `--refprot-max-evalue`, `--refprot-min-pident`, and `--refprot-max-hits` (all mirrored under `[refprot]`). These run after the refprot DIAMOND pass so you can keep the rescue set small without editing the DB.
  - Limit how many sequences from any single proteome survive panel selection with `--refprot-proteome-cap` (or `[refprot].proteome_cap`), complementing the lineage-based diversity caps.
