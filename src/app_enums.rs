use clap::ValueEnum;
use serde::{Deserialize, Serialize};

use crate::ecs::RenderOutputConfig;
use crate::mafft::AlignerBackend;

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize, Default)]
pub(crate) enum CalibrationMode {
    #[default]
    Off,
    Percentile,
    Isotonic,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
pub(crate) enum Mode {
    Centroid,
    Members,
    Auto,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
pub(crate) enum AlignmentStrategy {
    Auto,
    Add,
    AddFragments,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
pub(crate) enum AlignerCliBackend {
    #[serde(alias = "mafft")]
    Mafft,
    #[serde(alias = "spoa")]
    Spoa,
}

impl From<AlignerCliBackend> for AlignerBackend {
    fn from(v: AlignerCliBackend) -> Self {
        match v {
            AlignerCliBackend::Mafft => AlignerBackend::Mafft,
            AlignerCliBackend::Spoa => AlignerBackend::Spoa,
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
pub(crate) enum DiamondMode {
    Auto,
    Batch,
    Single,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
pub(crate) enum LogFormat {
    Text,
    Json,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum ReportFormat {
    Jsonl,
    Csv,
    Parquet,
    All,
}

impl ReportFormat {
    pub(crate) fn output_config(self, resume: bool) -> RenderOutputConfig {
        match self {
            ReportFormat::Jsonl => RenderOutputConfig {
                jsonl: true,
                csv: false,
                parquet: false,
                resume,
            },
            ReportFormat::Csv => RenderOutputConfig {
                jsonl: false,
                csv: true,
                parquet: false,
                resume,
            },
            ReportFormat::Parquet => RenderOutputConfig {
                jsonl: false,
                csv: false,
                parquet: true,
                resume,
            },
            ReportFormat::All => RenderOutputConfig {
                jsonl: true,
                csv: true,
                parquet: true,
                resume,
            },
        }
    }
}
