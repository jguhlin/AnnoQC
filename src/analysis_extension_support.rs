use crate::plugins::PluginDefinition;
use crate::rhai_rules;

pub struct LoadedExtensions {
    pub plugins: Vec<PluginDefinition>,
    pub rhai_runtime: Option<rhai_rules::RhaiRuntime>,
    pub plugin_count: usize,
    pub rule_count: usize,
}

pub struct FeatureSummary<'a> {
    pub profile_name: &'a str,
    pub has_genomic_input: bool,
    pub nucleotide_input: bool,
    pub taxonomy_enabled: bool,
    pub genomic_enabled: bool,
    pub plugins_enabled: bool,
    pub rules_enabled: bool,
    pub hmmer_enabled: bool,
    pub orphan_enabled: bool,
    pub alignment_enabled: bool,
}

pub fn load_extensions(plugin_paths: &[String], rhai_paths: &[String]) -> LoadedExtensions {
    let mut plugins = Vec::new();
    for path in plugin_paths {
        match PluginDefinition::load(path) {
            Ok(plugin) => {
                log::info!("loaded plugin: {}", plugin.name);
                plugins.push(plugin);
            }
            Err(e) => {
                log::warn!("failed to load plugin {}: {}", path, e);
            }
        }
    }

    let mut rhai_runtime = None;
    if !rhai_paths.is_empty() {
        match rhai_rules::RhaiRuntime::load(rhai_paths) {
            Ok(runtime) => {
                let names: Vec<_> = runtime.rules().iter().map(|r| r.name.clone()).collect();
                log::info!("loaded rhai rules: {}", names.join(", "));
                rhai_runtime = Some(runtime);
            }
            Err(e) => {
                log::warn!("failed to load rhai rules: {}", e);
            }
        }
    }

    LoadedExtensions {
        plugin_count: plugins.len(),
        rule_count: rhai_paths.len(),
        plugins,
        rhai_runtime,
    }
}

pub fn build_features_string(summary: &FeatureSummary<'_>) -> String {
    let mut features_list = vec!["Homology", "Intrinsic"];
    if summary.taxonomy_enabled {
        features_list.push("Taxonomy");
    }
    if summary.genomic_enabled {
        features_list.push("Genomic");
    }
    if summary.plugins_enabled {
        features_list.push("Plugins");
    }
    if summary.rules_enabled {
        features_list.push("Rules");
    }
    if summary.hmmer_enabled {
        features_list.push("Domains");
    }
    if summary.orphan_enabled {
        features_list.push("Orphan");
    }
    if summary.alignment_enabled {
        features_list.push("Alignment");
        features_list.push("Divergence");
    }

    let input_label = if summary.has_genomic_input {
        "GFF+Genome"
    } else if summary.nucleotide_input {
        "Nucleotide FASTA"
    } else {
        "Protein FASTA"
    };

    format!(
        "Input: {}, Profile: {}, Features: [{}]",
        input_label,
        summary.profile_name,
        features_list.join("+")
    )
}
