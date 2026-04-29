use std::{env, process};

use clap::Parser;
use mimalloc::MiMalloc;

use compute_matrix_rs::cli::Cli;
use compute_matrix_rs::pipeline;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> anyhow::Result<()> {
    guard_multi_letter_short_flags();

    let cli = Cli::parse();
    let config = cli.build_config()?;
    pipeline::execute(config)
}

fn guard_multi_letter_short_flags() {
    for arg in env::args_os().skip(1) {
        if arg == "-bs" || arg == "-bl" {
            eprintln!(
                "Error: multi-letter short flags '-bs'/'-bl' are not supported; use --bs/--binSize and --bl/--blackListFileName instead."
            );
            process::exit(2);
        }
    }
}
