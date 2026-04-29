pub use matrix::{MatrixHeader, MatrixRow};
pub mod core;
pub mod matrix;
mod reference_point;
mod run;
mod scale_regions;
pub mod zones;

pub(crate) use run::run_pipeline;

use crate::config::ModeConfig;

pub fn execute(config: crate::config::Config) -> anyhow::Result<()> {
    let general = &config.general;
    let io = &config.io;
    let gtf = &config.gtf;

    if !general.quiet {
        eprintln!(
            "Running computeMatrix (Rust) in {} mode.",
            describe_mode(&config.mode)
        );
        if general.verbose {
            eprintln!(
                "Processing {} region file(s) against {} score file(s).",
                io.regions.len(),
                io.scores.len()
            );
        }
    }

    match &config.mode {
        ModeConfig::ScaleRegions(options) => {
            scale_regions::run(general, io, gtf, options)?;
        }
        ModeConfig::ReferencePoint(options) => {
            reference_point::run(general, io, gtf, options)?;
        }
    }

    Ok(())
}

fn describe_mode(mode: &ModeConfig) -> &'static str {
    match mode {
        ModeConfig::ScaleRegions(_) => "scale-regions",
        ModeConfig::ReferencePoint(_) => "reference-point",
    }
}
