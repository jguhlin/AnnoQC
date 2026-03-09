mod build;
mod refprot;
pub(crate) mod support;

use crate::*;

pub(crate) fn run_prepare(p: PrepareArgs) -> Result<(), Box<dyn std::error::Error>> {
    let diamond_bin = p
        .diamond_bin
        .clone()
        .unwrap_or_else(|| "diamond".to_string());
    let prep_json = matches!(p.log_format, LogFormat::Json);
    build::run_primary_prepare_steps(&p, &diamond_bin, prep_json)?;
    refprot::run_aves_refprot_step(&p, &diamond_bin, prep_json)?;
    Ok(())
}
