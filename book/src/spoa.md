# SPOA Parity & Benchmarks

This page tracks parity checks and performance comparisons between SPOA (in-process MSA) and MAFFT (external MSA).

## Parity checklist

When SPOA is enabled, the pipeline applies a lightweight parity check before trusting the alignment. If any check fails, SPOA is treated as failed and the run falls back to MAFFT with a warning.

Checklist:

- Alignment has at least one row and a non-zero length.
- All aligned rows share the same length.
- The non-gap residue count for each row matches the original sequence length.
- If any of the above fails, log a warning and fall back to MAFFT.

Operational behavior:

- Warning line (example): `SPOA failed for <gene_id>: parity check failed: ...; falling back to MAFFT`.
- `strategy_used` becomes `mafft_fallback` in JSONL/CSV for those genes.

## Benchmarks

Use the benchmark script to compare MAFFT vs SPOA on representative datasets. The script runs the same input twice and records `run_metrics.json` for each run.

### How to run

```
./scripts/bench_aligners.sh --fasta a9_head200.faa --db uniprot_sprot.dmnd --out-root /tmp/annoqc_bench
```

### Results

Fill in the table after running the script. Keep runs on the same machine and config.

| Dataset | Threads | Aligner | Total (s) | MSA step (s) | Notes |
| --- | --- | --- | --- | --- | --- |
| a9_head200 | 8 | mafft | 34.854 | 0.000 | mafft_fast=true |
| a9_head200 | 8 | spoa | 34.546 | 0.000 | mafft_fast=true |
| a9_head1000 | 8 | mafft | 37.876 | 0.000 | mafft_fast=true |
| a9_head1000 | 8 | spoa | 37.073 | 0.000 | mafft_fast=true |
