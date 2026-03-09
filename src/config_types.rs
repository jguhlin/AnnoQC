use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::{AlignmentStrategy, CalibrationMode, DiamondMode, Mode, ReportFormat};

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(dead_code)]
pub(crate) struct FileConfig {
    pub(crate) fasta: Option<String>,
    pub(crate) db: Option<String>,
    pub(crate) threads: Option<usize>,
    pub(crate) approx_id: Option<u32>,
    pub(crate) member_cover: Option<u32>,
    pub(crate) out: Option<String>,
    pub(crate) report_format: Option<ReportFormat>,
    pub(crate) mode: Option<Mode>,
    pub(crate) top: Option<usize>,
    pub(crate) diamond_bin: Option<String>,
    pub(crate) hmmscan_bin: Option<String>,
    pub(crate) scoring: Option<ScoringConfigOverride>,
    pub(crate) reference_fasta: Option<String>,
    pub(crate) mafft_bin: Option<String>,
    pub(crate) alignment_top_hits: Option<usize>,
    pub(crate) alignment_strategy: Option<AlignmentStrategy>,
    pub(crate) conserved_identity_min: Option<f64>,
    pub(crate) alignment_missing_exon: Option<usize>,
    pub(crate) alignment_retained_intron: Option<usize>,
    pub(crate) mafft_threads_per_job: Option<usize>,
    pub(crate) mafft_max_jobs: Option<usize>,
    pub(crate) mafft_fast: Option<bool>,
    pub(crate) mafft_backend: Option<String>,
    pub(crate) render_max_jobs: Option<usize>,
    pub(crate) batch_size: Option<usize>,
    pub(crate) taxonomy: Option<TaxonomyConfigOverride>,
    pub(crate) taxonomy_cache: Option<String>,
    pub(crate) taxonomy_taxdump_dir: Option<String>,
    pub(crate) pfam_metadata: Option<String>,
    pub(crate) pfam_clans: Option<String>,
    pub(crate) pfam_db: Option<String>,
    pub(crate) diamond_mode: Option<DiamondMode>,
    pub(crate) hmmer: Option<HmmerConfigOverride>,
    pub(crate) diamond: Option<DiamondConfigOverride>,
    pub(crate) consensus: Option<ConsensusConfigOverride>,
    pub(crate) refprot: Option<RefProtConfigOverride>,
    pub(crate) structvar: Option<StructVarConfigOverride>,
    pub(crate) export_high: Option<bool>,
    pub(crate) export_high_path: Option<String>,
    pub(crate) calibration: Option<CalibrationConfigOverride>,
    pub(crate) rhai: Option<Vec<String>>,
    pub(crate) rnaseq: Option<RnaseqConfig>,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub(crate) struct RnaseqConfig {
    pub(crate) enabled: Option<bool>,
    pub(crate) file: Option<String>,
    pub(crate) min_tpm: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct TaxonomyConfigOverride {
    pub(crate) enabled: Option<bool>,
    pub(crate) min_support: Option<f64>,
    pub(crate) top_hits: Option<usize>,
    pub(crate) min_consensus: Option<usize>,
    pub(crate) coarse_rank_index: Option<usize>,
    pub(crate) coarse_min_support: Option<f64>,
    pub(crate) expected_domain: Option<String>,
    pub(crate) warn_non_target_min_frac: Option<f64>,
    pub(crate) warn_non_target_min_hits: Option<usize>,
    pub(crate) warn_non_target_strong_frac: Option<f64>,
    pub(crate) warn_non_target_strong_hits: Option<usize>,
    pub(crate) warn_genus_min_frac: Option<f64>,
    pub(crate) warn_genus_min_hits: Option<usize>,
    pub(crate) low_coverage_frac: Option<f64>,
    pub(crate) profile_db: Option<String>,
    #[serde(alias = "cache_path")]
    pub(crate) cache_path: Option<String>,
    #[serde(alias = "taxdump_dir")]
    pub(crate) taxdump_dir: Option<String>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct HmmerConfigOverride {
    pub(crate) top_n: Option<usize>,
    pub(crate) threads: Option<usize>,
    pub(crate) ievalue: Option<f64>,
    pub(crate) ref_ievalue: Option<f64>,
    pub(crate) orphan_analysis: Option<bool>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct DiamondConfigOverride {
    pub(crate) auto_threshold: Option<usize>,
    pub(crate) max_hsps: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct StructVarConfigOverride {
    pub(crate) min_hsp_len: Option<usize>,
    pub(crate) min_hsp_frac: Option<f64>,
    pub(crate) fusion_min_gap: Option<usize>,
    pub(crate) dup_max_gap: Option<usize>,
    pub(crate) split_delta: Option<f64>,
    pub(crate) min_subject_cov: Option<f64>,
    pub(crate) orient_majority: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct ConsensusConfigOverride {
    pub(crate) min_hits: Option<usize>,
    pub(crate) max_panel: Option<usize>,
    pub(crate) filt_qcov: Option<f64>,
    pub(crate) filt_scov: Option<f64>,
    pub(crate) filt_evalue: Option<f64>,
    pub(crate) filt_pident: Option<f64>,
    pub(crate) redundancy_pident: Option<f64>,
    pub(crate) max_high_identity: Option<usize>,
    pub(crate) len_ratio_tolerance: Option<f64>,
    pub(crate) backfill_enabled: Option<bool>,
    pub(crate) backfill_min_primary_hits: Option<usize>,
    pub(crate) backfill_max_added: Option<usize>,
    pub(crate) refprot_proteome_cap: Option<usize>,
    pub(crate) diversity_rank_index: Option<usize>,
    pub(crate) diversity_rank_cap: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct CalibrationConfigOverride {
    pub(crate) mode: Option<CalibrationMode>,
    pub(crate) min_samples: Option<usize>,
    pub(crate) min_unique: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct RefProtConfigOverride {
    pub(crate) enabled: Option<bool>,
    pub(crate) base_dir: Option<String>,
    pub(crate) readme_path: Option<String>,
    pub(crate) taxon_scope_rank: Option<String>,
    pub(crate) max_scopes: Option<usize>,
    pub(crate) max_proteomes: Option<usize>,
    pub(crate) trigger_k: Option<usize>,
    pub(crate) min_qcov: Option<f64>,
    pub(crate) min_scov: Option<f64>,
    pub(crate) max_evalue: Option<f64>,
    pub(crate) min_pident: Option<f64>,
    pub(crate) max_hits: Option<usize>,
    pub(crate) proteome_cap: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct ScoringConfigOverride {
    #[serde(default)]
    pub(crate) weights: HashMap<String, f64>,
    pub(crate) thresholds: Option<ScoringThresholds>,
    pub(crate) caps: Option<ScoringCapsConfigOverride>,
    pub(crate) genomic: Option<ScoringGenomicConfigOverride>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct ScoringThresholds {
    pub(crate) high: Option<f64>,
    pub(crate) medium: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct ScoringGenomicConfigOverride {
    pub(crate) min_canonical: Option<f64>,
    pub(crate) max_noncanonical: Option<f64>,
    pub(crate) max_weird: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub(crate) struct ScoringCapsConfigOverride {
    pub(crate) structvar_fusion_max: Option<f64>,
    pub(crate) structvar_split_max: Option<f64>,
    pub(crate) structvar_dup_max: Option<f64>,
}
