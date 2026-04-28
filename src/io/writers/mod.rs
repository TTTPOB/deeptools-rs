mod auxiliary;
mod formatting;
mod matrix_gz;

pub use matrix_gz::{
    StreamingMatrixWriter, build_padded_header_payload, ensure_streaming_header_capacity,
};

use anyhow::Result;

use crate::config::{IoOptions, SortRegions};
use crate::pipeline::matrix::MatrixData;

pub fn write_outputs(mut matrix: MatrixData, io: &IoOptions) -> Result<()> {
    if should_use_streaming(&matrix, io) {
        matrix_gz::write_matrix_gz_streaming(&io.matrix_output, &mut matrix)?;
        return Ok(());
    }

    matrix_gz::write_matrix_gz(&io.matrix_output, &matrix)?;

    if let Some(path) = &io.matrix_values_output {
        auxiliary::write_matrix_values(path, &matrix)?;
    }

    if let Some(path) = &io.sorted_regions_output {
        auxiliary::write_sorted_regions(path, &matrix)?;
    }

    Ok(())
}

fn should_use_streaming(matrix: &MatrixData, io: &IoOptions) -> bool {
    let sort_regions_str = &matrix.header.sort_regions;
    let sort_ok = sort_regions_str == "keep" || sort_regions_str == "no";
    should_use_streaming_for_plan(
        matrix.rows.len(),
        matrix.sample_count,
        matrix.bin_count,
        if sort_ok { SortRegions::Keep } else { SortRegions::Descend },
        io,
    )
}

pub fn should_use_streaming_for_plan(
    row_count: usize,
    sample_count: usize,
    bin_count: usize,
    sort_regions: SortRegions,
    io: &IoOptions,
) -> bool {
    if io.matrix_values_output.is_some() || io.sorted_regions_output.is_some() {
        return false;
    }

    if !matches!(sort_regions, SortRegions::Keep | SortRegions::No) {
        return false;
    }

    if row_count == 0 {
        return false;
    }

    let cell_count = row_count
        .saturating_mul(sample_count)
        .saturating_mul(bin_count);

    cell_count >= matrix_gz::STREAMING_CELL_THRESHOLD
}
