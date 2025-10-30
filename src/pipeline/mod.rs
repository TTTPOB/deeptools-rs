pub use matrix::{MatrixData, MatrixHeader, MatrixRow};
pub mod core;
pub mod matrix;
mod reference_point;
mod scale_regions;
pub mod zones;

use crate::config::{Config, IoOptions, ModeConfig};
use crate::io::writers;

pub enum RunOutcome {
    Matrix(MatrixData),
    Streamed,
}

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
        ModeConfig::ScaleRegions(options) => match scale_regions::run(general, io, options)? {
            RunOutcome::Matrix(matrix) => spawn_writer_thread(matrix, io)?,
            RunOutcome::Streamed => {}
        },
        ModeConfig::ReferencePoint(options) => match reference_point::run(general, io, options)? {
            RunOutcome::Matrix(matrix) => spawn_writer_thread(matrix, io)?,
            RunOutcome::Streamed => {}
        },
    }

    Ok(())
}

fn spawn_writer_thread(matrix: MatrixData, io: &IoOptions) -> anyhow::Result<()> {
    let io_clone = io.clone();
    let handle = std::thread::Builder::new()
        .name("matrix-writer".into())
        .spawn(move || writers::write_outputs(matrix, &io_clone))
        .map_err(|err| anyhow::anyhow!("Failed to spawn matrix writer thread: {err}"))?;

    match handle.join() {
        Ok(result) => result,
        Err(err) => std::panic::resume_unwind(err),
    }
}

fn describe_mode(mode: &ModeConfig) -> &'static str {
    match mode {
        ModeConfig::ScaleRegions(_) => "scale-regions",
        ModeConfig::ReferencePoint(_) => "reference-point",
    }
}
