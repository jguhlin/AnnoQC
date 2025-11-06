# External Reference Data Checklist

Use this list to provision all third-party datasets required for offline taxonomy scoring and domain-based evidence. Download into `share/` (or another path referenced by `--config`). Refresh on release updates or whenever UniProt/Pfam issue revisions.

## SwissProt Reference
- **File**: `uniprot_sprot.fasta.gz`
- **Source**: `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
- **Why**: Primary DIAMOND subject database; headers include `OX=` (NCBI TaxID) and lineage tokens.
- **Post-download**: Generate DIAMOND DB (`diamond makedb`), retain the FASTA for taxonomy parsing, and record release version/checksum in the run manifest.

## NCBI Taxonomy Tree
- **Archive**: `new_taxdump.tar.gz`
- **Source**: `https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz`
- **Contents**: `names.dmp`, `nodes.dmp`, `merged.dmp`, etc.
- **Usage**: Build in-memory TaxID → lineage/rank lookup when evaluating taxonomy consensus.
- **Notes**: Extract into `share/taxonomy/` and keep the accompanying `taxdump_readme.txt` with release dates.

## Pfam HMM Library
- **Models**: `Pfam-A.hmm.gz`
- **Metadata**: `Pfam-A.hmm.dat.gz`, `Pfam-A.clans.tsv.gz`
- **Seeds (optional)**: `Pfam-A.seed.gz`
- **Source**: `http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/`
- **Post-download**: Decompress `Pfam-A.hmm.gz`, run `hmmpress Pfam-A.hmm` to create `.h3{f,i,m,p}` files required by `hmmscan`. Unpack `.dat`/`.clans` for domain annotations and clade heuristics.

## Optional Supplements
- **SwissProt ID mapping** (`idmapping.dat.gz`) if cross-referencing UniProt accessions with RefSeq/Ensembl.
- **UniProt release notes** for provenance in manifests.
- **GTDB/PhyEco marker sets** if additional phylogenetic markers are desired.

## Verification & Storage
1. Store downloads under `share/` with subfolders: `share/uniprot/`, `share/taxonomy/`, `share/pfam/`.
2. Record SHA256 checksums in `share/CHECKSUMS.txt` to ensure reproducibility.
3. Update `config.example.toml` or your project config with paths:
   - `reference_fasta = "share/uniprot/uniprot_sprot.fasta.gz"`
   - `taxonomy.profile_db = "share/pfam/Pfam-A.hmm"` (after decompression).
   - Additional keys (planned): `taxonomy.taxdump_dir`, `taxonomy.cache_tsv`.
4. Keep raw archives alongside processed outputs so regeneration is repeatable.

All downloads should be performed manually or via `scripts/fetch_reference_data.sh` (see below); avoid hammering remote servers by caching artifacts in artifact storage or shared buckets.
