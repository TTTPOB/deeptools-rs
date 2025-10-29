use crate::config::{Config, ModeConfig};

pub fn execute(config: Config) -> anyhow::Result<()> {
    if !config.general.quiet {
        eprintln!(
            "Running computeMatrix (Rust) in {} mode.",
            describe_mode(&config.mode)
        );
        if config.general.verbose {
            eprintln!(
                "Processing {} region file(s) against {} score file(s).",
                config.io.regions.len(),
                config.io.scores.len()
            );
        }
    }

    match config.mode {
        ModeConfig::ScaleRegions(_) => {
            // TODO: implement scale-regions pipeline.
        }
        ModeConfig::ReferencePoint(_) => {
            // TODO: implement reference-point pipeline.
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
