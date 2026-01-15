# Quickstart

This is the fastest path to a working run on a small dataset.

## 1) Prepare reference DB

```
pixi run cargo run -- prepare
```

This downloads SwissProt and builds `uniprot_sprot.dmnd` plus clustering artifacts. Outputs are large; keep them out of git.

## 2) Analyze a FASTA

```
pixi run cargo run -- analyze --config config.example.toml
```

Outputs land in `results/`:

- `qc_report.jsonl` (rich per‑gene scorecards)
- `qc_summary.csv` (headline fields)
- `run.json` (manifest + config snapshot)

## 3) Optional: dry‑run scoring rubric

```
pixi run cargo run -- analyze --config config.example.toml --dry-run
```

This prints the resolved weights, thresholds, and calibration mode without running DIAMOND.

## 4) Optional: plugins and rules

```
pixi run cargo run -- analyze --config config.example.toml \
  --plugin plugins/my_plugin.wasm \
  --rhai rules/my_rule.rhai
```

## Notes

- Use `--report-format jsonl|csv|parquet|all` to control outputs.
- Use `--resume` to continue rendering into an existing output directory.
- If you have MAFFT and HMMER installed, add `mafft_bin` and `hmmscan_bin` to the config to enable those pillars.
