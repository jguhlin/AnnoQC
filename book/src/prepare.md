# Prepare

`cargo run -- prepare` builds the SwissProt DIAMOND database and (optionally) the Aves reference-proteome sidecar used for rescue panels. The command is checkpointed, so re-running with `--resume` skips steps whose `.done` marker already exists.

## Download caching & resume

Every HTTP artifact (`*.fasta.gz`, Pfam files, the UniProt reference-proteomes README, etc.) now drops two helper files next to the download:

- `filename.part` holds in-progress bytes so that interrupted transfers can resume via HTTP `Range`/`If-Range`.
- `filename.httpmeta` caches the most recent `ETag`/`Last-Modified` headers. When you rerun `prepare`, we issue conditional requests (`If-None-Match` / `If-Modified-Since`) and skip transfers that the server reports as unmodified. This keeps the SwissProt + Pfam payloads warm even if the step is repeated daily.

Partial downloads are automatically resumed, and empty/corrupt files trigger a retry/backoff loop before the run fails. You can tune the retry behavior with `--download-retries N` (default 3 exponential backoff attempts).

## Reference proteomes (Aves)

When `scripts/fetch_reference_data.sh` has seeded the UniProt reference-proteomes README + taxdump, `prepare` automatically selects every Aves proteome listed in the README, downloads the corresponding FASTA, concatenates them into `share/refprot/aves/aves_refprot.fasta.gz`, and builds a taxonomy-aware DIAMOND database (`aves_refprot.dmnd`).

- `--refprot-workers N` controls how many concurrent download workers run (default 4). Each worker logs per-proteome progress so you can see which accessions were refreshed.
- `--download-retries N` sets the retry count per proteome (default 3). Failures back off exponentially (2s, 4s, 8s, …) before giving up.

The downloader only rebuilds the concatenated FASTA/DIAMOND database when new proteomes arrive or when the aggregate files are missing. Existing databases are auto-verified with `diamond dbinfo`; if taxonomy metadata is absent we delete the `.done` marker and rebuild with `--taxonmap/--taxonnodes/--taxonnames` so downstream analyses can safely emit `staxids`/`slineages`.

SwissProt makedb receives the same auto-repair treatment: if `uniprot_sprot.dmnd` exists but lacks taxonomy entries, `prepare` deletes `uniprot_sprot.dmnd.done` and rebuilds the DB in-place with the cached accession→taxid map plus the downloaded taxdump.

## Taxonomy counts CLI

Need a quick answer to "how many reference proteomes sit under Aves?" Run the new helper:

```bash
cargo run -- taxonomy-count \
  --name Aves \
  --rank order \
  --top 15
```

By default the subcommand scans `share/uniprot/reference_proteomes/README` and the cached taxdump; you can override both paths via `--readme` and `--taxdump-dir`. Provide either `--name` or `--taxid` to pick the target lineage, and optionally `--rank <level>` to group descendants by a particular rank (e.g., `order`, `family`). The tool reports the total number of proteomes at/under the target lineage plus the top-N descendants at the requested rank so you can gauge coverage before kicking off a large run.
