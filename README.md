# AnnoQC

[![CI](https://github.com/jguhlin/AnnoQC/workflows/CI/badge.svg)](https://github.com/jguhlin/AnnoQC/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**AnnoQC** is a production-grade tool for protein sequence quality control. It combines homology search, intrinsic sequence metrics, conserved-region alignment, taxonomy contamination detection, and genomic context to identify problematic gene models in annotated genomes.

## Features

- **Homology-based QC**: DIAMOND searches against SwissProt with intelligent panel selection
- **Intrinsic quality**: Ambiguity, low complexity, and homopolymer detection
- **Conserved regions**: MAFFT/SPOA alignment-based completeness assessment
- **Taxonomy**: Contamination detection and consensus building
- **Genomic context**: Intron/exon structure validation
- **Extensible**: WASM plugins (Extism) and Rhai scripting support
- **Fast**: Parallel processing using Bevy ECS
- **Rich output**: JSONL, CSV, and Parquet with detailed scorecards

## Quick Start

### Installation

#### Option 1: Pixi (reproducible environment, recommended)

```bash
pixi install annoqc
pixi run annoqc --help
```

#### Option 2: Precompiled binary

Download from [GitHub Releases](https://github.com/jguhlin/AnnoQC/releases):

```bash
wget https://github.com/jguhlin/AnnoQC/releases/latest/download/annoqc-linux-x86_64.tar.gz
tar -xzf annoqc-linux-x86_64.tar.gz
chmod +x annoqc
```

See [INSTALLATION.md](INSTALLATION.md) for more installation options including Docker and Cargo.

### First Run

1. **Download reference data** (~1GB):
   ```bash
   annoqc prepare
   ```

2. **Configure** (copy and edit `config.example.toml`):
   ```toml
   [general]
   fasta = "my_proteins.faa"
   db = "uniprot_sprot.dmnd"
   out = "results"

   [general.resources]
   threads = 16
   ```

3. **Analyze**:
   ```bash
   annoqc analyze --config config.toml
   ```

4. **View results**:
   - `results/qc_report.jsonl`: Detailed per-gene scorecards
   - `results/qc_summary.csv`: Headline metrics
   - `results/run.json`: Run manifest and config snapshot

## Documentation

Full documentation: [https://jguhlin.github.io/AnnoQC](https://jguhlin.github.io/AnnoQC)

- [Quickstart Guide](https://jguhlin.github.io/AnnoQC/quickstart.html)
- [Installation Guide](https://jguhlin.github.io/AnnoQC/installation.html)
- [Scoring Methods](https://jguhlin.github.io/AnnoQC/scoring.html)
- [Plugin Development](https://jguhlin.github.io/AnnoQC/plugins.html)
- [Genomic Context](https://jguhlin.github.io/AnnoQC/genomic_context.html)
- [Performance Tuning](https://jguhlin.github.io/AnnoQC/performance.html)

## Example Output

```json
{
  "gene_id": "Gene00042",
  "final_score": 0.87,
  "quality_category": "high",
  "homology_score": 0.92,
  "intrinsic_score": 0.95,
  "taxonomy_score": 0.78,
  "domains_score": 0.85,
  "warnings": [
    "Low DIAMOND hit count (3 hits)"
  ],
  "panel_provenance": {
    "swissprot": 8,
    "refprot": 2,
    "cluster": 0
  }
}
```

## Performance

| Dataset | Genes | Threads | Runtime | Memory |
|---------|-------|---------|---------|--------|
| Small   | 100   | 4       | ~2min   | 2GB    |
| Medium  | 10K   | 16      | ~30min  | 8GB    |
| Large   | 100K  | 32      | ~4hr    | 16GB   |

*Benchmarks on DIAMOND homology search, SwissProt database, Ubuntu 22.04*

See [performance.md](https://jguhlin.github.io/AnnoQC/performance.html) for detailed tuning guidance.

## Citation

If you use AnnoQC in your research, please cite:

```
Guhlin J. AnnoQC: Production-grade protein sequence quality control.
GitHub repository: https://github.com/jguhlin/AnnoQC
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

Built with [Bevy ECS](https://bevyengine.org/), [DIAMOND](https://github.com/bbuchfink/diamond), [MAFFT](https://mafft.cbrc.jp/alignment/software/), and [HMMER](http://hmmer.org/).
