# Remaining Work: RNA-seq Support Implementation

## Status

**Plan Approved:** RNA-seq Support Implementation Plan saved at `/home/guhjo98p/.claude/plans/linear-splashing-wombat.md`

**Progress:** Plan approved, implementation not yet started.

**Conversation Summary:**
- Completed: SPOA alignment integration (✅)
- Completed: ORF translation multi-frame support (✅)
- In Progress: RNA-seq support (plan approved, ready to implement)

## Overview

Add RNA-seq as a new optional evidence pillar to AnnoQC. When RNA-seq data is provided for a run, genes without RNA-seq support will receive a penalty (low RNA-seq score), while genes with strong expression evidence will receive a boost.

**Key Design Principle**: RNA-seq is "run-provided" - if present in the run configuration, it contributes to the denominator for ALL genes, but only contributes to the numerator for genes with actual expression support.

## Implementation Approach

Support two tabular input formats:
- **Option A**: Quant tables from Salmon/Kallisto (TPM/NumReads)
- **Option C**: Precomputed per-gene expression support JSON/TSV

## Critical Files to Modify

| File | Purpose | Status |
|------|---------|--------|
| `src/rnaseq.rs` | NEW: RNA-seq metrics struct and parsing logic | Not started |
| `src/scoring.rs` | Add `compute_rnaseq_score()` function | Not started |
| `src/main.rs` | CLI args, config, scoring integration, output | Not started |
| `src/explain.rs` | Add RNA-seq display to explain command | Not started |
| `config.example.toml` | RNA-seq configuration examples | Not started |

## Implementation Tasks

### Task 1: Create RNA-seq Module (`src/rnaseq.rs`)

**Status:** Not started

Create new file with:
- `RnaseqMetrics` struct (tpm, num_reads, effective_length, has_support, expression_score)
- `parse_quant_table()` function for Salmon/Kallisto TSV format
- `parse_expression_file()` function for JSON/TSV support files
- `parse_expression_json()` and `parse_expression_tsv()` helpers
- `normalize_tpm()` function for log-scale TPM normalization [0, 1]
- Unit tests for TPM normalization

**Code template:** See plan file at `/home/guhjo98p/.claude/plans/linear-splashing-wombat.md` lines 28-213

### Task 2: Add Scoring Function (`src/scoring.rs`)

**Status:** Not started

Add function:
```rust
pub fn compute_rnaseq_score(rnaseq: Option<&RnaseqMetrics>) -> f64 {
    let Some(rnaseq) = rnaseq else {
        return 0.0;
    };
    rnaseq.expression_score.clamp(0.0, 1.0)
}
```

### Task 3: Add CLI Arguments (`src/main.rs`)

**Status:** Not started

Add to `AnalyzeArgs` struct (around line 500-520):
- `rnaseq_file: Option<String>` - RNA-seq data file path
- `rnaseq_min_tpm: f64` - Minimum TPM threshold (default 1.0)

Add to `FileConfig` struct (around line 600):
- `rnaseq: Option<RnaseqConfig>`

Add `RnaseqConfig` struct (around line 650):
```rust
#[derive(Deserialize, Default, Clone)]
struct RnaseqConfig {
    enabled: Option<bool>,
    file: Option<String>,
    min_tpm: Option<f64>,
}
```

### Task 4: Update Scoring System (`src/main.rs`)

**Status:** Not started

**A. Add to `WeightSet` struct** (around line 4300):
```rust
struct WeightSet {
    // ... existing fields ...
    rnaseq: f64,  // NEW
}
```

**B. Update `scoring_weights()`** (around line 4304):
- Set default `rnaseq: 0.0`
- Apply user config from `scoring.weights.get("rnaseq")`

**C. Update `WeightSet::sum()` method**:
- Add `self.rnaseq` to sum calculation

### Task 5: Add PillarPresence Field

**Status:** Not started

Update `PillarPresence` struct (around line 4400):
```rust
struct PillarPresence {
    // ... existing fields ...
    has_rnaseq: bool,  // NEW
}
```

### Task 6: Integrate into `build_scores_map()`

**Status:** Not started

In scoring loop (around line 4700):
```rust
let (rnaseq_score, rnaseq_present) = if rnaseq_enabled {
    let evidence = rnaseq_map.as_ref().and_then(|m| m.get(&m.gene_id));
    let score = compute_rnaseq_score(evidence);
    let present = rnaseq_map.as_ref().map_or(false, |m| m.contains_key(&m.gene_id));
    (score, present)
} else {
    (0.0, false)
};

// In numerator calculation:
if presence.rnaseq {
    numerator += weights.rnaseq * rnaseq_score;
}
```

### Task 7: Update ScoreCard Struct

**Status:** Not started

Update `ScoreCard` struct (around line 200):
```rust
pub struct ScoreCard {
    // ... existing fields ...
    pub rnaseq_score: f64,
    pub rnaseq_tpm: String,
    pub rnaseq_num_reads: String,
}
```

### Task 8: Update Output Formats

**Status:** Not started

**JSONL** (around line 5500):
```rust
"rnaseq": {
    "tpm": rnaseq_metrics.as_ref().and_then(|r| r.tpm),
    "num_reads": rnaseq_metrics.as_ref().and_then(|r| r.num_reads),
    "expression_score": rnaseq_metrics.as_ref().map(|r| r.expression_score).unwrap_or(0.0),
},
"score_components": {
    "rnaseq": rnaseq_score,
    // ... other scores ...
}
```

**CSV** (around line 5860, 5920):
- Add `rnaseq_score` column
- Add `rnaseq_tpm` column
- Add `rnaseq_num_reads` column

### Task 9: Update AnalysisContext

**Status:** Not started

Update `AnalysisContext` struct (around line 150):
```rust
struct AnalysisContext {
    // ... existing fields ...
    rnaseq_enabled: bool,
    rnaseq_map: Arc<HashMap<String, RnaseqMetrics>>,
}
```

### Task 10: Parse RNA-seq Input

**Status:** Not started

Add loading logic in main analyze function (around line 1200-1400):
```rust
let (rnaseq_enabled, rnaseq_map) = if let Some(rnaseq_path) = args.rnaseq_file
    .or_else(|| file_cfg.rnaseq.as_ref().and_then(|r| r.file.clone()))
{
    log::info!("loading RNA-seq data from: {}", rnaseq_path);
    let map = if rnaseq_path.ends_with(".tsv") || rnaseq_path.ends_with(".txt") {
        rnaseq::parse_quant_table(&rnaseq_path)
            .or_else(|_| rnaseq::parse_expression_file(&rnaseq_path))
    } else {
        rnaseq::parse_expression_file(&rnaseq_path)
    };
    match map {
        Ok(m) => {
            log::info!("loaded RNA-seq data for {} genes", m.len());
            (true, Arc::new(m))
        }
        Err(e) => {
            log::warn!("failed to load RNA-seq data: {}", e);
            (false, Arc::new(HashMap::new()))
        }
    }
} else {
    (false, Arc::new(HashMap::new()))
};
```

### Task 11: Update Explain Command

**Status:** Not started

Add to `src/explain.rs`:
```rust
if let Some(rnaseq) = val.get("rnaseq").and_then(|v| v.as_object()) {
    if let Some(score) = rnaseq.get("expression_score").and_then(|v| v.as_f64()) {
        println!("RNA-seq expression score: {:.4}", score);
    }
    if let Some(tpm) = rnaseq.get("tpm").and_then(|v| v.as_f64()) {
        println!("RNA-seq TPM: {:.2}", tpm);
    }
}
```

## Verification Steps

1. **Unit tests:** `cargo test rnaseq`
   - Test TPM normalization
   - Test quant table parsing
   - Test JSON/TSV parsing
   - Test scoring function

2. **Integration test:**
   ```bash
   annoqc analyze --fasta proteins.faa --db dmnd.dmnd --rnaseq-file quant.tsv --out results/
   ```
   - Verify RNA-seq scores appear in output
   - Check that genes without support get lower scores
   - Verify CSV/JSONL output has new columns

3. **Score calculation verification:**
   - Genes with TPM >= 100 should score 1.0
   - Genes with TPM < 1 should score 0.0
   - Genes without RNA-seq data should score 0.0

4. **Backward compatibility:**
   - Runs without `--rnaseq-file` should work identically to before
   - Existing tests should still pass

## Important Notes

### Module Declaration
Don't forget to add to `src/main.rs`:
```rust
mod rnaseq;
use rnaseq::{RnaseqMetrics, parse_quant_table, parse_expression_file};
```

### Architectural Pattern
RNA-seq follows the same pattern as other optional pillars (taxonomy, genomic, domains):
- Run-provided: contributes to denominator for ALL genes when enabled
- Gene-provided: contributes to numerator only for genes with actual data
- Missing data = 0.0 score (not "N/A")

### TPM Normalization
- TPM < 1.0 → score 0.0 (no expression support)
- TPM >= 100.0 → score 1.0 (max expression)
- Uses log scale for smooth interpolation

## Example Input Files

**Salmon/Kallisto Quant Table (TSV):**
```
Name	Length	EffectiveLength	TPM	NumReads
gene1	1500	1450	50.5	500
gene2	2000	1950	0.5	5
gene3	1800	1750	120.0	1200
```

**Expression Support JSON:**
```json
{
  "gene1": {"tpm": 50.5, "num_reads": 500, "support_score": 0.7},
  "gene2": {"tpm": 0.5, "num_reads": 5, "support_score": 0.1},
  "gene3": {"tpm": 120.0, "num_reads": 1200, "support_score": 0.95}
}
```

**Expression Support TSV:**
```
gene_id	tpm	num_reads	support_score
gene1	50.5	500	0.7
gene2	0.5	5	0.1
gene3	120.0	1200	0.95
```

## Configuration Example

Add to `config.example.toml`:
```toml
[rnaseq]
enabled = true
file = "path/to/quant.tsv"
min_tpm = 1.0

[scoring.weights]
homology = 0.5
intrinsic = 0.25
taxonomy = 0.0
genomic = 0.0
domains = 0.0
structvar = 0.0
conserved = 0.0
rnaseq = 0.25  # New weight
```

## Quick Start for New Conversation

When starting a new conversation, provide this context:

```
I'm working on AnnoQC RNA-seq support. See REMAINING_WORK.md for details.
The implementation plan is at /home/guhjo98p/.claude/plans/linear-splashing-wombat.md

Please start with Task 1: Create src/rnaseq.rs module.
```

## Recent Completed Work

For reference, recently completed features:
- SPOA alignment integration (✅ Complete)
- ORF translation multi-frame support (✅ Complete)
- Both have tests passing

These are separate from RNA-seq work.
