mod diff;
mod header;
mod parse;
mod values;

use std::path::PathBuf;
use std::process;

use clap::{Parser, Subcommand};

use diff::{full_diff, print_full_diff, print_value_result};
use header::compare_headers;
use parse::load_matrix;
use values::compare_values;

/// Exit codes:
///   0 — match
///   1 — mismatch
///   2 — error (I/O or parse failure)
const EXIT_MATCH: i32 = 0;
const EXIT_MISMATCH: i32 = 1;
const EXIT_ERROR: i32 = 2;

#[derive(Parser)]
#[command(
    name = "compare_matrix",
    about = "Dev tool: compare two computeMatrix output files",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Compare only the JSON header lines.
    Header {
        left: PathBuf,
        right: PathBuf,
        /// Header keys to ignore during comparison.
        #[arg(long = "ignore", value_name = "KEY")]
        ignore: Vec<String>,
    },
    /// Compare only the numeric values (BED fields are used for alignment but not checked).
    Values {
        left: PathBuf,
        right: PathBuf,
        /// Absolute tolerance for floating-point comparison [default: 1e-6].
        #[arg(long = "tolerance", default_value = "1e-6")]
        tolerance: f64,
        /// Maximum number of value diffs to print [default: 10].
        #[arg(long = "max-diffs", default_value = "10")]
        max_diffs: usize,
    },
    /// Full diff: header + BED field ordering + numeric values.
    Diff {
        left: PathBuf,
        right: PathBuf,
        /// Absolute tolerance for floating-point comparison [default: 1e-6].
        #[arg(long = "tolerance", default_value = "1e-6")]
        tolerance: f64,
        /// Header keys to ignore during comparison.
        #[arg(long = "ignore", value_name = "KEY")]
        ignore: Vec<String>,
        /// Maximum number of value diffs to print [default: 10].
        #[arg(long = "max-diffs", default_value = "10")]
        max_diffs: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    let exit_code = match run(cli) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("error: {:#}", e);
            EXIT_ERROR
        }
    };

    process::exit(exit_code);
}

fn run(cli: Cli) -> anyhow::Result<i32> {
    match cli.command {
        Command::Header {
            left,
            right,
            ignore,
        } => {
            let left_matrix = load_matrix(&left)?;
            let right_matrix = load_matrix(&right)?;

            let diffs =
                compare_headers(&left_matrix.header_json, &right_matrix.header_json, &ignore);

            if diffs.is_empty() {
                println!(
                    "MATCH — headers are identical (ignoring {} key(s))",
                    ignore.len()
                );
                Ok(EXIT_MATCH)
            } else {
                println!("MISMATCH — {} header difference(s):", diffs.len());
                for d in &diffs {
                    println!("  key {:?}:", d.key);
                    println!("    left:  {}", d.left);
                    println!("    right: {}", d.right);
                }
                Ok(EXIT_MISMATCH)
            }
        }

        Command::Values {
            left,
            right,
            tolerance,
            max_diffs,
        } => {
            let left_matrix = load_matrix(&left)?;
            let right_matrix = load_matrix(&right)?;

            let result =
                compare_values(&left_matrix.rows, &right_matrix.rows, tolerance, max_diffs);

            print_value_result(&result);

            if result.matches {
                Ok(EXIT_MATCH)
            } else {
                Ok(EXIT_MISMATCH)
            }
        }

        Command::Diff {
            left,
            right,
            tolerance,
            ignore,
            max_diffs,
        } => {
            let left_matrix = load_matrix(&left)?;
            let right_matrix = load_matrix(&right)?;

            let diff = full_diff(
                &left_matrix.header_json,
                &right_matrix.header_json,
                &left_matrix.rows,
                &right_matrix.rows,
                tolerance,
                &ignore,
                max_diffs,
            );

            print_full_diff(
                &diff,
                &left.display().to_string(),
                &right.display().to_string(),
            );

            if diff.matches {
                Ok(EXIT_MATCH)
            } else {
                Ok(EXIT_MISMATCH)
            }
        }
    }
}
