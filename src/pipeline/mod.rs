pub use matrix::{MatrixData, MatrixHeader, MatrixRow};
pub mod matrix;
mod reference_point;

use crate::config::{Config, ModeConfig};
use crate::io::writers;

pub fn execute(config: Config) -> anyhow::Result<()> {
    let general = &config.general;
    let io = &config.io;

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
        ModeConfig::ScaleRegions(_options) => {
            // TODO: implement scale-regions pipeline.
        }
        ModeConfig::ReferencePoint(options) => {
            let matrix = reference_point::run(general, io, options)?;
            writers::write_outputs(&matrix, io)?;
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
