# Future Roadmap: AnnoQC v2.0

This document outlines the strategic technical expansion of AnnoQC to move beyond "Protein FASTA Quality" into "Genomic Integrity" and "Extensible Architecture."

## 1. The Rule Engine (Scriptable Pillar)
*Goal: Allow users to define custom logic without recompiling.*

We will embed **Rhai**, a tiny scripting language for Rust, to allow users to define "soft rules" in their `config.toml`.

### Usage
```toml
[scripts]
# Define a custom penalty rule
my_check = """
    if intrinsic.gc_content < 0.2 {
        return 0.5; // Penalty
    }
    if homology.top_hit.contains("Uncharacterized") {
        return 0.2;
    }
    return 0.0;
"""
```

### Implementation
1.  **Dependency:** Add `rhai = "1.10"` to `Cargo.toml`.
2.  **Context:** Expose the `GeneMetrics`, `IntrinsicMetrics`, and `DiamondHitStats` structs to the Rhai engine.
3.  **Execution:** In the `scoring` system (ECS), compile the user scripts once and execute them for every gene.
4.  **Output:** The return value (0.0 to 1.0) is added to a new `custom_penalty` field in the final score.

## 2. WASM Plugin System
*Goal: Allow community-developed, high-performance modules (e.g., SignalP, TMHMM wrappers).*

For complex logic that requires heavy computation or proprietary algorithms, we will support WebAssembly (WASM) plugins via **Wasmer** or **Wasmtime**.

### Architecture
*   **Host (AnnoQC):** Provides an API (`annoqc_get_sequence(id)`, `annoqc_emit_score(key, val)`).
*   **Guest (Plugin):** A `.wasm` file that exports `analyze_gene()`.

### Workflow
1.  User provides `--plugin my_signal_peptide.wasm`.
2.  AnnoQC loads the WASM module.
3.  For each gene, AnnoQC calls the plugin.
4.  The plugin processes the sequence and returns a JSON object with scores/metadata.
5.  AnnoQC merges this into the `custom` block of the JSONL/Parquet output.

## 3. Genomic Context (GFF3 Integration)
*Goal: Validate the "Gene Model" structure, not just the protein sequence.*

This is the biggest scientific leap. It moves AnnoQC from "Protein QC" to "Annotation QC."

### The "Prepare-GFF" Command
New subcommand: `annoqc prepare-gff --gff annotation.gff3 --genome genome.fna`

### Features to Extract
1.  **Splice Site Validity:**
    *   Iterate all CDS features for a gene.
    *   Extract the dinucleotides at the intron boundaries from the genome.
    *   **Metric:** `% Canonical Splice Sites (GT-AG, GC-AG, AT-AC)`.
    *   **Penalty:** Penalize non-canonical sites heavily (indicates frameshift or assembly error).
2.  **Intron Statistics:**
    *   Mean intron length.
    *   Micro-introns (<10bp, usually errors).
3.  **Start/Stop Context:**
    *   Does the CDS start with ATG/GTG/TTG in the genome?
    *   Does it end with a valid stop codon?
4.  **Overlap/Antisense:**
    *   Detect "Shadow Genes" (ORFs on the opposite strand of a real gene).

### Integration
*   The `prepare-gff` command produces a `genomic_context.jsonl` (or Parquet) file indexed by `gene_id`.
*   The `analyze` command accepts `--genomic-context` as an input.
*   ECS System: A new `GenomicPillar` reads this data and adds it to the ScoreCard.

## 4. Golden Set Auto-Tuning
*Goal: Remove the guesswork from weighting.*

1.  **Input:** User provides `golden.faa` (high-confidence genes from their organism).
2.  **Training:** AnnoQC runs the full pipeline on `golden.faa`.
3.  **Distribution:** It builds a statistical profile (Mean/Stdev) of Intrinsic scores, Homology density, etc.
4.  **Application:**
    *   Sets `thresholds.high` to the 10th percentile of the Golden Set.
    *   Adjusts weights to maximize the score of the Golden Set.
    *   Saves this as a `[profile]` for the actual run.

## 5. Summary of Priorities

| Feature | Effort | Impact | Status |
| :--- | :--- | :--- | :--- |
| **GFF3/Genomic** | High | **Very High** | Planned |
| **Rule Engine (Rhai)** | Low | Medium | Planned |
| **WASM Plugins** | High | High (Community) | Future |
| **Auto-Tuning** | Medium | Medium | Future |

**Recommendation:** Start with **GFF3 Integration**. It addresses the most common critique of protein-based QC ("You can't see splice site errors").