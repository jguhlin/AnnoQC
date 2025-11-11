  RUST_LOG=info pixi run cargo run -- \
    --config config.weights_demo.toml \
    analyze \
    --fasta a9.faa \
    --db uniprot_sprot.dmnd \
    --out results_a9_full_mafftfast \
    --threads 12 \
    --log-format json \
    --diamond-mode single \
    --mafft-bin mafft \
    --reference-fasta share/uniprot/uniprot_sprot.fasta.gz \
    --hmmscan-bin hmmscan \
    --csv-verbose \
    --mafft-fast \
    --render-max-jobs 16
