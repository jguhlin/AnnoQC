use std::fs;

use crate::config_types::FileConfig;
use crate::{
    resolve_calibration_settings, resolve_effective_config, AnalyzeArgs, CalibrationSettings,
    EffectiveConfig, ReportFormat,
};

pub(crate) struct AnalyzeBootstrap {
    pub(crate) file_cfg: FileConfig,
    pub(crate) cfg: EffectiveConfig,
    pub(crate) calibration: CalibrationSettings,
    pub(crate) report_format: ReportFormat,
    pub(crate) rhai_paths: Vec<String>,
}

pub(crate) fn bootstrap_analyze(
    cli_config_path: Option<&str>,
    args: &AnalyzeArgs,
) -> Result<AnalyzeBootstrap, Box<dyn std::error::Error>> {
    let config_path = args.config.as_deref().or(cli_config_path);
    let mut file_cfg = if let Some(path) = config_path {
        let text = fs::read_to_string(path)?;
        toml::from_str::<FileConfig>(&text)?
    } else {
        FileConfig::default()
    };

    let mut scoring = file_cfg.scoring.unwrap_or_default();
    args.profile.apply_to(&mut scoring);
    file_cfg.scoring = Some(scoring);

    let cfg = resolve_effective_config(&file_cfg, args)?;
    let calibration = resolve_calibration_settings(args, &file_cfg);
    let report_format = args
        .report_format
        .or(file_cfg.report_format)
        .unwrap_or(ReportFormat::All);
    let mut rhai_paths = file_cfg.rhai.clone().unwrap_or_default();
    if !args.rhai.is_empty() {
        rhai_paths.extend(args.rhai.clone());
    }

    Ok(AnalyzeBootstrap {
        file_cfg,
        cfg,
        calibration,
        report_format,
        rhai_paths,
    })
}
