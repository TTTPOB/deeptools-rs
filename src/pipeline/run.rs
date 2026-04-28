use std::sync::Arc;

use anyhow::Result;

use crate::config::{GeneralOptions, GtfOptions, IoOptions};
use crate::io::writers;
use crate::pipeline::core::{
    self, FileCollector, InMemoryCollector, PipelineMode, RegionTask,
};
use crate::pipeline::matrix::MatrixHeader;

use super::RunOutcome;

pub fn run_pipeline<M>(
    mode: M,
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
) -> Result<RunOutcome>
where
    M: PipelineMode + Clone + Send + 'static,
    M::Metadata: Clone + 'static,
{
    let metadata = Arc::new(mode.validate(general)?);

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = sample_labels.len();

    let groups = core::load_groups(&io.regions, gtf)?;
    let group_labels: Vec<String> = groups.iter().map(|g| g.label.clone()).collect();
    let group_capacity: Vec<usize> = groups.iter().map(|g| g.records.len()).collect();

    let mut tasks = Vec::new();
    for (group_index, group) in groups.into_iter().enumerate() {
        for record in group.records {
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: Arc::new(record),
            });
        }
    }

    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let total_bins = mode.total_bins(metadata.as_ref());
    let row_count = tasks.len();
    let should_stream = writers::should_use_streaming_for_plan(
        row_count,
        sample_count,
        total_bins,
        general.sort_regions,
        io,
    );

    let sample_paths = Arc::new(io.scores.clone());

    if should_stream {
        let header_estimate = mode.build_header(
            general,
            metadata.as_ref(),
            &sample_labels,
            &group_labels,
            &group_capacity,
            thread_count,
            sample_count,
        );
        writers::ensure_streaming_header_capacity(&header_estimate)?;

        let writer = writers::StreamingMatrixWriter::start(&io.matrix_output)?;
        let collector = FileCollector::new(writer);
        let header_builder = {
            let general = general.clone();
            let sample_labels = sample_labels.clone();
            let group_labels = group_labels.clone();
            let metadata = Arc::clone(&metadata);
            let mode = mode.clone();
            move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
                Ok(mode.build_header(
                    &general,
                    metadata.as_ref(),
                    &sample_labels,
                    &group_labels,
                    &group_counts,
                    thread_count,
                    sample_count,
                ))
            }
        };

        core::execute_mode(
            tasks,
            general,
            Arc::clone(&sample_paths),
            collector,
            thread_count,
            &mode,
            Arc::clone(&metadata),
            header_builder,
            group_labels.len(),
        )?;
        return Ok(RunOutcome::Streamed);
    }

    let collector = InMemoryCollector::with_capacity(row_count, sample_count, total_bins);
    let header_builder = {
        let general = general.clone();
        let sample_labels = sample_labels.clone();
        let group_labels = group_labels.clone();
        let metadata = Arc::clone(&metadata);
        let mode = mode.clone();
        move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
            Ok(mode.build_header(
                &general,
                metadata.as_ref(),
                &sample_labels,
                &group_labels,
                &group_counts,
                thread_count,
                sample_count,
            ))
        }
    };

    let mut matrix = core::execute_mode(
        tasks,
        general,
        sample_paths,
        collector,
        thread_count,
        &mode,
        metadata,
        header_builder,
        group_labels.len(),
    )?;

    let sort_sample_indices =
        core::normalize_sort_sample_indices(general.sort_using_samples.as_ref(), sample_count)?;

    matrix.sort_groups(
        general.sort_regions,
        general.sort_using,
        sort_sample_indices.as_deref(),
    )?;

    matrix.prune_zero_rows();

    Ok(RunOutcome::Matrix(matrix))
}
