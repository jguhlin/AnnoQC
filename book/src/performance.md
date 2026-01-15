% Performance & Tuning

This page collects the practical knobs we watch when scaling AnnoQC from toy subsets to proteome-wide runs.

## DIAMOND

- **Modes**: `--diamond-mode single` runs a single blastp pass; `auto` switches to chunked mode when the query FASTA exceeds `[diamond].auto_threshold` (default 200k sequences). Chunked mode reuses the same DB and logs per-chunk JSON progress, which is helpful on clusters.
- **Sensitivity vs throughput**: the default config requests `--sensitive`, `--max-target-seqs 50`, and disables motif masking. When you need faster screening, drop to `--mid-sensitive` or lower `max-target-seqs`. When you need more homologs (e.g., sparse proteomes), you can bump `max-target-seqs` to 100 and enable `--very-sensitive`; expect MAFFT/HMMER load to grow accordingly.
- **Outfmt tokens**: we pass each field token separately (`--outfmt 6 qseqid … staxids slineages`). If you modify the list, keep `qlen/slen/staxids` so that coverage and taxonomy pillars continue to function.
- **Retries**: both `prepare` and `analyze` wrap DIAMOND invocations in retry loops (three attempts) and log stderr to `diamond.log`. Use `RUST_LOG=debug` to see command lines when debugging failed runs.

## Reference proteomes

- Run `cargo run -- prepare --resume` periodically so the Aves sidecar database stays in sync. The downloader now caches ETags / Last-Modified headers, splits transfers across `--refprot-workers` (default 4), and retries each proteome up to `--download-retries` times with exponential backoff.
- During `analyze`, refprot hits are only merged when SwissProt contributes fewer than `refprot_trigger_k` panel members. Adjust `[consensus].refprot_proteome_cap` if a single proteome floods the panel; `[consensus].diversity_rank_cap` limits per-class diversity so that bird panels spread across multiple clades automatically.

## MAFFT & HMMER pipelines

- `--mafft-fast` (or `mafft_fast=true` in TOML) enforces FFT-NS-1/300 settings, which keeps alignments short without sacrificing termini concordance signals.
- `mafft_threads_per_job` × `mafft_max_jobs` controls concurrency. The default auto-allocates `(threads-1)` cores across MAFFT workers while leaving one CPU for orchestrating ECS. For example, with `--threads 16` we typically run `mafft_threads_per_job = 2`, `mafft_max_jobs = 6`.
- `render_max_jobs` caps the number of JSON/CSV render tasks in flight. If disk I/O becomes the bottleneck, lower this value; if emit lags behind MAFFT/HMMER, increase it.
- HMMER obeys `[hmmer].threads` (defaults to the global `--threads`). Raising `--hmmer-top-n` dramatically increases domtblout size; keep it ≤200 unless you need exhaustive domain dumps for debugging.

## ECS scheduling & observability

- The scheduler sorts gene entities by `length + 100*hit_count` and then feeds them into async pools, so long, hit-rich genes start early. Use `RUST_LOG=info` to watch `step_start`, `mafft progress`, and `hmmscan progress` JSON lines. They fire every 5% or 60s (whichever comes first) to avoid spamming multi-hour runs.
- `panel_debug.csv`, `panel_agg_debug.csv`, and `panel_sources.csv` provide ground truth for how many homologs were scheduled, filtered, or backfilled per gene. Inspect these when MAFFT fails to trigger as often as expected.

## Quick recipes

- **200-seq sanity**: `pixi run cargo run -- --config config.weights_demo.toml analyze --fasta a9_head200.faa --db uniprot_sprot.dmnd --out results_head200 --threads 12 --log-format json --mafft-fast --hmmscan-bin hmmscan --csv-verbose` (≈2 min on a 16-core workstation).
- **Full proteome with fast MAFFT**: `RUST_LOG=info pixi run cargo run -- --config config.weights_demo.toml analyze --fasta a9.faa --db uniprot_sprot.dmnd --out results_full --threads 32 --mafft-fast --mafft-threads-per-job 2 --mafft-max-jobs 12 --render-max-jobs 16 --hmmscan-bin hmmscan --csv-verbose` (≈15 min after DIAMOND).
- **Debugging slow downloads**: rerun `cargo run -- prepare --download-retries 5 --refprot-workers 2` with `RUST_LOG=debug` to see each proteome’s URL and retry backoff.

Adjust these templates to match your hardware; the defaults aim to keep MAFFT/HMMER busy while leaving enough headroom for streaming emit and file I/O.
