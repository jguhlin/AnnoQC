# Plugins & Rules

AnnoQC supports two extension paths:

- WASM plugins (Extism) loaded with `--plugin`.
- Rhai rules loaded with `--rhai`.

Both receive a per-gene JSON payload and can return a score and/or penalty. Penalties are subtracted from `final_score`.

## WASM Plugins (Extism)

### Contract

- Export a function named `analyze`.
- Input: JSON string representing the gene input snapshot.
- Output: JSON with fields:
  - `name` (string, required)
  - `score` (number, optional)
  - `penalty` (number, optional)
  - `metadata` (object, optional)

Example output:

```json
{"name":"CysteineCheck","score":0.9,"metadata":{"cysteines":4}}
```

### Input shape (abbreviated)

```json
{
  "gene_id": "Gene00042",
  "sequence": "MKT...",
  "homology": {
    "hits_count": 12,
    "top_hit": "sp|P12345|...",
    "top_bitscore": 212.6,
    "top_evalue": "3.2e-68",
    "top_qcov": 0.94,
    "top_scov": 0.91,
    "bitscore_density": 3.21,
    "coverage_delta": 0.03,
    "coverage_ratio": 1.03
  },
  "intrinsic": {
    "ambiguous_fraction": 0.0,
    "max_homopolymer": 5,
    "low_complexity_fraction": 0.08,
    "low_complexity_windows": 2,
    "orf_start_score": 1.0
  },
  "taxonomy": {
    "detail": "Consensus",
    "congruence_score": 0.83,
    "contamination_score": 0.17,
    "support_fraction": 0.85,
    "support": 17,
    "considered": 20,
    "consensus_rank": "family",
    "consensus_taxid": 561,
    "consensus_name": "Escherichia"
  },
  "panel": {"swissprot": 6, "refprot": 2, "cluster": 0},
  "genomic": {
    "introns_total": 4,
    "splice_canonical": 3,
    "splice_noncanonical": 1,
    "splice_weird": 0,
    "intron_len_min": 72,
    "intron_len_max": 410,
    "intron_len_avg": 221.5
  }
}
```

Fields may be `null` if the corresponding pillar is disabled or has no data.

### CLI usage

```
annoqc analyze --fasta input.faa --db uniprot_sprot.dmnd \
  --plugin path/to/plugin.wasm
```

You can pass `--plugin` multiple times to load multiple modules.

### Output integration

Plugin results appear in:

- JSONL: `plugins.total_penalty` and `plugins.results[]`.
- CSV/Parquet: `plugin_penalty`, `plugin_names`, `plugin_scores`, `plugin_penalties` columns.
- `run.json`: `plugins[]` manifest with file hashes.

Failures are logged at `debug` and do not stop the run.

### Optional metadata sidecar

If `path/to/plugin.wasm.json` exists and contains `name` or `version`, those values are recorded in `run.json` for traceability:

```json
{"name":"CysteineCheck","version":"1.2.0"}
```

## Rhai Rules

Rhai scripts are lightweight rules that share the same input shape as WASM plugins.

### CLI usage

```
annoqc analyze --fasta input.faa --db uniprot_sprot.dmnd \
  --rhai rules/penalize_short.rhai
```

Rhai rules return either:

- A map with `name/score/penalty/metadata`, or
- A numeric score.

Multiple rules may be provided; penalties are summed.

## Scoring impact

Plugin penalties are subtracted from the calibrated score:

```
final_score = clamp(final_score_raw - total_penalty, 0.0, 1.0)
```

Use `annoqc explain <gene_id>` to inspect per-gene plugin results.
