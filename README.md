# AnnoQC

[![CI](https://github.com/jguhlin/AnnoQC/workflows/CI/badge.svg)](https://github.com/jguhlin/AnnoQC/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

AnnoQC is a Rust CLI for protein and gene-model quality control. The current codebase is a multi-stage analysis pipeline built around DIAMOND homology search, consensus-panel construction, optional heavy evidence stages (MAFFT or SPOA, HMMER/Pfam), contextual evidence (taxonomy, genomic context, RNA-seq, structural-variation heuristics), and streaming report generation.

## Status

- Implemented and code-backed today: `prepare`, `analyze`, `explain`, `taxonomy-cache`, `refprot-index`, and `taxonomy-count`.
- The analysis pipeline is substantially broader than the older README/docs implied. It includes config/profile overrides, tool preflight, DIAMOND auto/batch/single modes, consensus/refprot rescue, streaming JSONL/CSV/Parquet output, run manifests, debug sidecars, WASM plugins, and Rhai rules.
- The main gaps are no longer "missing pipeline" items. They are mostly scoring/config refinement and output semantics:
  - ECS intake is still lightweight bookkeeping rather than rich per-gene computation.
  - Parquet output cannot append on `--resume`.
  - Some scoring knobs and sample config sections are only partially wired.
  - A few classification/reporting behaviors are still heuristic.

See [ARCHITECTURE.md](ARCHITECTURE.md), [NEXT_STEPS_TODO.md](NEXT_STEPS_TODO.md), and [FUTURE_TODO.md](FUTURE_TODO.md) for the deeper design and gap ledger.

## Current CLI

- `prepare`: build DIAMOND and clustering artifacts from an existing SwissProt FASTA, and refresh/build the Aves refprot sidecar when the UniProt README + taxonomy inputs are present.
- `analyze`: run the full QC pipeline.
- `explain <gene_id>`: read `qc_report.jsonl` from a prior run and print a per-gene breakdown.
- `taxonomy-cache`: build an accession-to-taxid cache from a FASTA.
- `refprot-index`: inspect/build a reference-proteome index from the UniProt README + taxdump.
- `taxonomy-count`: summarize lineage coverage in the reference-proteome catalog.

## Recommended Workflow

### 1. Create the tool environment

Pixi is the most accurate local workflow because it places DIAMOND, HMMER, and MAFFT on `PATH`.

```bash
pixi run cargo build
```

### 2. Fetch reference inputs

`prepare` does not fetch SwissProt by itself. The current repo workflow is to seed reference data first:

```bash
./scripts/fetch_reference_data.sh
```

That script downloads and stages:

- `share/uniprot/uniprot_sprot.fasta.gz`
- `share/taxonomy/new_taxdump`
- `share/pfam/Pfam-A.hmm` and related Pfam files
- `share/uniprot/reference_proteomes/README`

### 3. Build derived reference artifacts

```bash
pixi run cargo run -- prepare \
  --fasta share/uniprot/uniprot_sprot.fasta.gz \
  --db-out share/uniprot/uniprot_sprot.dmnd \
  --resume
```

This produces or refreshes:

- `share/uniprot/uniprot_sprot.dmnd`
- `clusters`, `clusters.realign`, `clusters.recluster`
- `share/refprot/aves/aves_refprot.fasta.gz`
- `share/refprot/aves/aves_refprot.dmnd`

`prepare` is checkpointed with `.done` markers and validates taxonomy metadata in rebuilt DIAMOND databases.

### 4. Run analysis

The easiest starting point is to copy `config.example.toml` and update the paths for your environment.

Minimal direct CLI example:

```bash
pixi run cargo run -- analyze \
  --fasta tests/fixtures/input.faa \
  --db share/uniprot/uniprot_sprot.dmnd \
  --reference-fasta share/uniprot/uniprot_sprot.fasta.gz \
  --out results
```

Config-driven example:

```bash
pixi run cargo run -- analyze --config config.example.toml
```

Useful flags:

- `--report-format jsonl|csv|parquet|all`
- `--resume`
- `--dry-run`
- `--aligner mafft|spoa`
- `--plugin path/to/plugin.wasm`
- `--rhai path/to/rule.rhai`
- `--gff annotations.gff3 --genome genome.fa`
- `--rnaseq-file expression.tsv`

### 5. Inspect a gene

```bash
pixi run cargo run -- explain Gene00042 --out results
```

## Implemented Feature Set

### Core analysis

- DIAMOND preflight and DB validation before analysis starts.
- DIAMOND execution modes: `auto`, `single`, and `batch`.
- Optional nucleotide input translation before protein QC.
- Intrinsic metrics including ambiguity, homopolymers, low-complexity windows, and ORF diagnostics.
- Consensus panel construction from DIAMOND HSP aggregation.
- Panel provenance tracking across `swissprot`, `refprot`, and `cluster` sources.
- Length expectation summaries derived from the selected panel.
- Structural-variation heuristics from DIAMOND tiling patterns.

### Optional heavy evidence

- Multiple-sequence alignment via external MAFFT or in-process SPOA.
- Alignment-derived metrics including conserved fraction, pairwise identity, divergence, termini concordance, missing/extra conserved block warnings, and exon/intron-like run heuristics.
- HMMER `hmmscan` support against Pfam.
- Domain-architecture scoring and `domains_arch_debug.csv`.
- Orphan-domain analysis for N/C-terminal partial domain fragments.

### Optional contextual evidence

- Taxonomy consensus and contamination/congruence summaries from DIAMOND hits.
- Taxonomy cache generation and reference-proteome lineage utilities.
- Genomic context scoring from GFF plus genome FASTA, including splice-site and intron statistics.
- RNA-seq expression-file parsing and per-gene support metrics.
- RefProt rescue panels and Aves reference-proteome sidecar generation.

### Scoring and extensibility

- Weighted multi-pillar scoring with threshold-based classification.
- Calibration modes: `off`, `percentile`, `isotonic`.
- Scoring profiles and TOML overrides.
- WASM plugins via Extism.
- Rhai scripting rules as a lighter-weight extension path.
- Optional export of high-confidence sequences to FASTA.

### Rendering and observability

- Streaming render pipeline instead of building one giant in-memory result map.
- Output formats: JSONL, CSV, Parquet.
- Resume support for JSONL and CSV.
- Run manifest with tool versions, checksums, config snapshot, plugins, and rules.
- Run summary, run metrics, slowest-gene report, and sidecar debug artifacts.

## Support Matrix

This table is meant to answer three practical questions quickly: is the feature implemented, is it exercised by automation, and what extra data or binaries does it need?

| Feature | Release status | Coverage | Requirements / notes |
| --- | --- | --- | --- |
| `prepare` artifact build | Supported | Some automated coverage | Needs DIAMOND plus a pre-fetched SwissProt FASTA; taxdump/refprot inputs unlock richer outputs |
| Core `analyze` pipeline | Supported | Smoke tested | Smoke test uses a DIAMOND stub and verifies the main JSONL/CSV path |
| DIAMOND modes (`auto`/`single`/`batch`) | Supported | Partial | Core behavior is covered; mode-specific end-to-end coverage is still lighter than ideal |
| JSONL output | Supported | Smoke tested | Standard output path |
| CSV output | Supported | Smoke tested | Standard output path |
| Parquet output | Supported | Manual/external validation recommended | Implemented in the render pipeline; parquet append on `--resume` is not supported |
| `--dry-run` scoring rubric | Supported | Tested | No external reference data needed beyond normal CLI startup |
| `explain` subcommand | Supported | Manual/external validation recommended | Reads an existing `qc_report.jsonl` and prints a per-gene breakdown |
| SPOA alignment backend | Supported | Smoke tested plus unit-tested helpers | In-process backend; does not require MAFFT |
| MAFFT alignment backend | Supported | Unit-tested helpers | Requires MAFFT plus a matching reference FASTA; explicit end-to-end MAFFT smoke coverage is still limited |
| HMMER / Pfam domain analysis | Supported | Unit-tested parsing and scoring | Requires `hmmscan` plus a Pfam database |
| Orphan-domain analysis | Supported | Unit-tested | Runs as part of the HMMER/Pfam stage |
| Consensus panel selection / length scoring | Supported | Unit-tested | Core part of the main analysis pipeline |
| Taxonomy consensus / contamination | Supported | Unit-tested with lighter end-to-end coverage | Requires taxonomy metadata, cache, or taxdump-backed resolution |
| Genomic context scoring | Supported | Unit-tested | Requires GFF plus genome FASTA, ideally with `.fai` indexing |
| RNA-seq parsing / support metrics | Supported | Unit-tested parsing/normalization | Requires an expression input file; scoring semantics are still an active refinement area |
| Structural-variation heuristics | Supported | Unit-tested | Derived from DIAMOND HSP tiling patterns |
| WASM plugins | Supported | Unit-tested | Requires plugin `.wasm` modules |
| Rhai rules | Supported | Manual/external validation recommended | Requires `.rhai` rule files; implemented and wired into scoring/rendering |
| Resume for JSONL/CSV | Supported | Manual/external validation recommended | Reuses existing outputs by loading already-rendered gene IDs |
| `taxonomy-cache` / `refprot-index` / `taxonomy-count` | Supported | Manual/external validation recommended | Utility commands that depend on FASTA and/or UniProt/taxdump inputs |

## Output Files

Primary outputs in `results/`:

- `qc_report.jsonl`: rich per-gene JSON records.
- `qc_summary.csv`: per-gene tabular summary.
- `qc_summary.parquet`: Parquet summary when enabled.
- `run.json`: manifest, resolved settings, tool versions, checksums, plugin/rule metadata.
- `run_summary.json`: run-level counts, feature flags, throughput.
- `run_metrics.json`: step durations and aggregate totals.
- `slowest_genes.json`: slowest alignment/HMMER jobs.

Common sidecars/debug outputs:

- `panel_debug.csv`
- `panel_agg_debug.csv`
- `panel_sources.csv`
- `domains_arch_debug.csv`
- `structvar_summary.json`
- `refprot_selected.txt`

## Testing

Current test coverage is broader than the old README implied:

- Integration smoke tests run `analyze` with a compiled DIAMOND stub and verify JSONL/CSV output plus the SPOA path.
- `analyze --dry-run` is tested.
- Unit tests exist across scoring, taxonomy, genomic context, consensus, HMMER, MAFFT, plugins, RNA-seq, ORF logic, and related helpers.

Recommended local commands:

```bash
pixi run cargo test
pixi run cargo fmt --all
pixi run cargo clippy -- -D warnings
```

## Project Maturity

If "done" means "usable software release for other groups", the project is fairly close. The pipeline and architecture are already in place; the remaining work is mostly scoring-policy cleanup, output/schema consistency, full test validation, and deciding which features are supported versus still experimental.

Roughly, the current state looks like this:

- software/tooling maturity: about 75-85% of the way to a solid release
- paper readiness: about 40-60% of the way to a defensible manuscript

The difference matters. The software is close because the core pipeline already exists. A paper is farther away because it still needs a proper validation package:

- benchmark datasets with curated truth or strong proxies
- baseline comparisons against relevant tools
- ablation studies across evidence pillars
- threshold/calibration evaluation
- runtime and scaling measurements
- case studies and error analysis

Practical estimate:

- "done" as a tool: likely a few focused weeks of cleanup and validation
- "done" as a paper: more like a dedicated benchmarking and writing phase after the tool is stabilized

## Known Gaps / TODO

The most credible current TODOs are in [NEXT_STEPS_TODO.md](NEXT_STEPS_TODO.md). The main code-verified gaps are:

- ORF/start/stop diagnostics are computed but not fully folded into the final score.
- Conserved-block and alignment scoring still need further tuning for GV-style missing/extra block behavior.
- Structural-variation penalty tuning is still evolving.
- Domain order/architecture penalties are not fully weighted in final scoring.
- Some sample config sections are ahead of implementation, especially parts of `[scoring.homology]` and `[scoring.intrinsic]`.
- Classification/reporting semantics are still being refined around plugin penalties and the `X` bucket.
- Parquet append on `--resume` is not supported.
- The ECS intake stage is still a lightweight scheduler placeholder.

The larger roadmap in [FUTURE_TODO.md](FUTURE_TODO.md) still calls out:

- stronger RNA-seq scoring semantics at run scope
- a clearer denominator model for "expected" versus "optional" evidence pillars

## Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md)
- [INSTALLATION.md](INSTALLATION.md)
- [book/](book/)
- [BEVY_GUIDE.md](BEVY_GUIDE.md)

## License

MIT. See [LICENSE](LICENSE).
