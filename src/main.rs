use clap::Parser;

use compute_matrix_rs::cli::Cli;
use compute_matrix_rs::pipeline;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    let config = cli.build_config()?;
    pipeline::execute(config)
}
