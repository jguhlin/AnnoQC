use extism::{Manifest, Plugin, Wasm};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Serialize, Clone)]
pub struct PluginInput {
    pub gene_id: String,
    pub sequence: String,
    pub homology: Option<PluginHomology>,
    pub intrinsic: PluginIntrinsic,
    pub taxonomy: Option<PluginTaxonomy>,
    pub panel: Option<PluginPanel>,
    pub genomic: Option<PluginGenomic>,
}

#[derive(Serialize, Clone)]
pub struct PluginHomology {
    pub hits_count: usize,
    pub top_hit: Option<String>,
    pub top_bitscore: f64,
    pub top_evalue: String,
    pub top_qcov: f64,
    pub top_scov: f64,
    pub bitscore_density: f64,
    pub coverage_delta: f64,
    pub coverage_ratio: f64,
}

#[derive(Serialize, Clone, Default)]
pub struct PluginIntrinsic {
    pub ambiguous_fraction: f64,
    pub max_homopolymer: usize,
    pub low_complexity_fraction: f64,
    pub low_complexity_windows: usize,
    pub orf_start_score: f64,
}

#[derive(Serialize, Clone)]
pub struct PluginTaxonomy {
    pub detail: String,
    pub congruence_score: f64,
    pub contamination_score: f64,
    pub support_fraction: f64,
    pub support: usize,
    pub considered: usize,
    pub consensus_rank: Option<String>,
    pub consensus_taxid: Option<u32>,
    pub consensus_name: Option<String>,
}

#[derive(Serialize, Clone)]
pub struct PluginPanel {
    pub swissprot: usize,
    pub refprot: usize,
    pub cluster: usize,
}

#[derive(Serialize, Clone)]
pub struct PluginGenomic {
    pub introns_total: usize,
    pub splice_canonical: usize,
    pub splice_noncanonical: usize,
    pub splice_weird: usize,
    pub intron_len_min: usize,
    pub intron_len_max: usize,
    pub intron_len_avg: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PluginOutput {
    pub name: String,
    pub score: Option<f64>,
    pub penalty: Option<f64>,
    pub metadata: Option<serde_json::Value>,
}

#[derive(Clone)]
pub struct PluginDefinition {
    pub name: String,
    manifest: Arc<Manifest>,
}

impl PluginDefinition {
    pub fn load(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let bytes = std::fs::read(path)?;
        let name = std::path::Path::new(path)
            .file_stem()
            .ok_or_else(|| format!("plugin path '{}' has no file stem", path))?
            .to_string_lossy()
            .to_string();
        let manifest = Manifest::new([Wasm::data(bytes)]);
        Ok(Self {
            name,
            manifest: Arc::new(manifest),
        })
    }
}

const MAX_PLUGIN_OUTPUT_BYTES: usize = 1024 * 1024;

/// Runs the plugin's `analyze` entrypoint using `PluginInput` and returns `PluginOutput`.
/// Output JSON must include `name` and may include `score`, `penalty`, and `metadata`.
pub fn run_plugin(
    def: &PluginDefinition,
    input: &PluginInput,
) -> Result<PluginOutput, Box<dyn std::error::Error>> {
    let mut plugin = Plugin::new(def.manifest.as_ref(), [], true)?;

    let input_json = serde_json::to_string(input)?;
    let output_bytes = plugin.call::<&str, &str>("analyze", &input_json)?;

    if output_bytes.len() > MAX_PLUGIN_OUTPUT_BYTES {
        return Err(format!(
            "plugin {} output exceeded {} bytes",
            def.name, MAX_PLUGIN_OUTPUT_BYTES
        )
        .into());
    }

    let output: PluginOutput = serde_json::from_str(&output_bytes)
        .map_err(|err| format!("plugin {} produced invalid JSON: {}", def.name, err))?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn mock_wasm_plugin_returns_output() {
        let output = r#"{"name":"MockPlugin","score":0.75}"#;
        let bytes = output.as_bytes();
        let output_len = bytes.len();
        let mut stores = String::new();
        for (idx, b) in bytes.iter().enumerate() {
            let offset = idx as i64;
            stores.push_str(&format!(
                "    (call $store_u8 (i64.add (local.get $ptr) (i64.const {})) (i32.const {}))\n",
                offset, b
            ));
        }
        let wat = format!(
            r#"(module
                (import "extism:host/env" "alloc" (func $alloc (param i64) (result i64)))
                (import "extism:host/env" "store_u8" (func $store_u8 (param i64 i32)))
                (import "extism:host/env" "output_set" (func $output_set (param i64 i64)))
                (memory (export "memory") 1)
                (func (export "analyze") (result i32)
                    (local $ptr i64)
                    (local.set $ptr (call $alloc (i64.const {})))
{}
                    (call $output_set (local.get $ptr) (i64.const {}))
                    (i32.const 0)
                )
            )"#,
            output_len, stores, output_len
        );
        let wasm = wat::parse_str(&wat).expect("compile wat");
        let manifest = Manifest::new([Wasm::data(wasm)]);
        let def = PluginDefinition {
            name: "mock".to_string(),
            manifest: Arc::new(manifest),
        };
        let input = PluginInput {
            gene_id: "gene1".to_string(),
            sequence: "ACGTACGT".to_string(),
            homology: None,
            intrinsic: PluginIntrinsic::default(),
            taxonomy: None,
            panel: None,
            genomic: None,
        };
        let out = run_plugin(&def, &input).expect("run plugin");
        assert_eq!(out.name, "MockPlugin");
        assert_eq!(out.score, Some(0.75));
    }
}
