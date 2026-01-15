# Extism Plugin System Plan

**Goal:** Enable users to extend AnnoQC with custom logic written in Python, JavaScript, Rust, or Go, without recompiling the main binary.

**Why Extism?** It abstracts away the complexity of WASM memory management and provides robust SDKs (PDKs) for multiple languages, making it the most user-friendly choice for bioinformatics customization.

## 1. Architecture

*   **Host (AnnoQC):**
    *   Loads the `.wasm` plugin file provided by the user.
    *   Serializes the current gene's data (Metrics, Sequence, Context) into a JSON string.
    *   Calls the plugin's exported function (default: `analyze`).
    *   Deserializes the JSON response and merges it into the `ScoreCard`.

*   **Guest (Plugin):**
    *   Reads the input JSON.
    *   Performs arbitrary logic (e.g., regex search, GC content calc, specific motif check).
    *   Returns a JSON object with optional scores, penalties, and metadata.

## 2. Data Contract (JSON)

### Input (Host -> Guest)
```json
{
  "gene_id": "Gene123",
  "sequence": "MKT...",
  "metrics": {
    "homology": { "bitscore": 100.0, "qcov": 0.9 },
    "intrinsic": { "ambiguous": 0.0 }
  },
  "genomic": {
    "introns": 2,
    "splice_canonical": 2
  }
}
```

### Output (Guest -> Host)
```json
{
  "status": "Success",  // or "Skip", "Error"
  "name": "MyCustomCheck",
  "score": 0.8,         // Optional: Direct score contribution
  "penalty": 0.0,       // Optional: Subtracted from final score
  "metadata": {         // Arbitrary data to log
    "signal_peptide_detected": true
  }
}
```

## 3. Implementation Plan

### Phase A: Dependencies & Setup
1.  Add `extism = "1.0"` to `Cargo.toml`.
2.  Add `serde_json` (already present, verify features).

### Phase B: Module `src/plugins.rs`
1.  Define the `PluginManager` struct.
2.  Implement `load_plugin(path: &str) -> Plugin`.
3.  Implement `run_plugin(plugin, data: &GeneData) -> PluginResult`.

### Phase C: Integration
1.  **CLI:** Add `--plugin <path>` argument to `AnalyzeArgs` (supports multiple).
2.  **Config:** Support `[plugins]` section in `config.toml` to load by default.
3.  **ECS:**
    *   Create a `PluginSystem` that runs after metrics are calculated but before scoring.
    *   Or, integrate into the `scoring` loop directly for simplicity.
4.  **Output:** Add a `custom_scores` map to the `ScoreCard` and columns to the Parquet/JSONL output.

### Phase D: "Batteries Included" (Examples)
1.  Create `examples/plugins/python_gc_content/` with a Python script and instructions on how to compile it to WASM (using `extism-py` or `componentize-py` if applicable, though Extism usually has its own tooling).
    *   *Note:* Python-to-WASM currently often requires a simplified runtime or specific Extism PDK usage. JS/Rust are often easier "first examples".

## 4. Todo List

- [ ] **Dependency:** Add `extism` crate.
- [ ] **Structs:** Define `PluginInput` and `PluginOutput` structs (derive Serialize/Deserialize).
- [ ] **Logic:** Implement `src/plugins.rs`.
- [ ] **CLI:** Add `--plugin` flag.
- [ ] **Wiring:** Call plugins inside the `render_gene_record` or a dedicated system.
- [ ] **Output:** Extend `ScoreCard` to hold plugin results.
- [ ] **Documentation:** Add `PLUGINS.md` explaining how to write a plugin.

## 5. Example Python Plugin Logic (Conceptual)

```python
import extism
import json

def analyze():
    input_str = extism.input_str()
    data = json.loads(input_str)
    
    seq = data.get("sequence", "")
    if "Cys" in seq: # Pseudo-code
        result = {
            "name": "CysteineCheck",
            "score": 1.0,
            "metadata": { "has_cysteine": True }
        }
    else:
        result = {
            "name": "CysteineCheck",
            "penalty": 0.1,
            "metadata": { "has_cysteine": False }
        }
        
    extism.output_str(json.dumps(result))
```
