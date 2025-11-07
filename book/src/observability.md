% Observability

Analyze emits structured progress and a compact metrics sidecar to help monitor runs:

- Step events (JSON mode)
  - `{"event":"step_start","name":"diamond"}`
  - `{"event":"step_finish","name":"diamond","seconds":"12.34"}`
  - Emitted for: `diamond`, `ecs`, `intrinsic`, `mafft`, `hmmer`, `scoring`, `emit_outputs`.

- Chunk events (DIAMOND batch)
  - `{"event":"diamond_chunk","chunk_index":0,"queries":8,"processed":8,"seconds":"0.50","rate":"16.00"}`

- Sidecar: `run_metrics.json`
  - Schema: `{ "schema_version": "1.0", "steps": [{"name":"ecs","seconds":1.23}, ...], "totals": { ... } }`
  - Totals include: `total_genes`, `mafft_alignments`, `hmmer_genes`, `hmmer_hits_total`, `diamond_mode`.

Text mode uses simple `step_start:` / `step_finish:` log lines; chunk lines are suppressed.

