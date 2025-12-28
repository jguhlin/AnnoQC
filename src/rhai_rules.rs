use std::path::Path;
use std::sync::Arc;

use rhai::serde::{from_dynamic, to_dynamic};
use rhai::{Dynamic, Engine, Scope, AST};

use crate::plugins::{PluginInput, PluginOutput};

/// A compiled Rhai rule with its source path metadata.
#[derive(Clone)]
pub struct RhaiRule {
    pub name: String,
    #[allow(dead_code)]
    pub path: String,
    pub ast: AST,
}

/// A shared Rhai runtime that compiles and evaluates a set of rules.
#[derive(Clone)]
pub struct RhaiRuntime {
    engine: Arc<Engine>,
    rules: Vec<RhaiRule>,
}

impl RhaiRuntime {
    /// Load and compile Rhai scripts from the provided paths.
    ///
    /// Returns an error if any file cannot be read or compiled.
    pub fn load(paths: &[String]) -> Result<Self, String> {
        let mut engine = Engine::new();
        engine.set_max_expr_depths(64, 64);
        engine.set_max_call_levels(32);
        let engine = Arc::new(engine);
        let mut rules = Vec::new();
        for path in paths {
            let script =
                std::fs::read_to_string(path).map_err(|e| format!("read {}: {}", path, e))?;
            let ast = engine
                .compile(&script)
                .map_err(|e| format!("compile {}: {}", path, e))?;
            let name = Path::new(path)
                .file_stem()
                .or_else(|| Path::new(path).file_name())
                .map(|name| name.to_string_lossy().to_string())
                .unwrap_or_else(|| path.clone());
            rules.push(RhaiRule {
                name,
                path: path.clone(),
                ast,
            });
        }
        Ok(Self { engine, rules })
    }

    /// Return the compiled rules in load order.
    pub fn rules(&self) -> &[RhaiRule] {
        &self.rules
    }

    /// Evaluate all loaded rules for a single plugin input.
    pub fn run(&self, input: &PluginInput) -> Vec<PluginOutput> {
        if self.rules.is_empty() {
            return Vec::new();
        }
        let input_dyn = match to_dynamic(input) {
            Ok(v) => v,
            Err(e) => {
                log::debug!("rhai input serialize failed: {}", e);
                return Vec::new();
            }
        };
        let mut out = Vec::with_capacity(self.rules.len());
        let mut scope = Scope::new();
        for rule in &self.rules {
            scope.clear();
            scope.push_dynamic("input", input_dyn.clone());
            let res = self.engine.call_fn::<Dynamic>(
                &mut scope,
                &rule.ast,
                "analyze",
                (input_dyn.clone(),),
            );
            let res = match res {
                Ok(v) => v,
                Err(e) => {
                    log::debug!("rhai rule {} failed: {}", rule.name, e);
                    continue;
                }
            };
            let mut parsed = if let Ok(val) = from_dynamic::<PluginOutput>(&res) {
                val
            } else if res.is::<f64>() {
                PluginOutput {
                    name: rule.name.clone(),
                    score: Some(res.cast::<f64>()),
                    penalty: None,
                    metadata: None,
                }
            } else if res.is::<i64>() {
                PluginOutput {
                    name: rule.name.clone(),
                    score: Some(res.cast::<i64>() as f64),
                    penalty: None,
                    metadata: None,
                }
            } else {
                PluginOutput {
                    name: rule.name.clone(),
                    score: None,
                    penalty: None,
                    metadata: None,
                }
            };
            if parsed.name.is_empty() {
                parsed.name = rule.name.clone();
            }
            out.push(parsed);
        }
        out
    }
}
