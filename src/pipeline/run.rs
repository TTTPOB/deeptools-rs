use std::sync::Arc;

use anyhow::Result;

use crate::config::{GeneralOptions, GtfOptions, IoOptions, SortRegions};
use crate::io::writers;
use crate::pipeline::core::{self, FileCollector, PipelineMode, RegionTask};
use crate::pipeline::matrix::MatrixHeader;

pub fn run_pipeline<M>(
    mode: M,
    general: &GeneralOptions,
    io: &IoOptions,
    gtf: &GtfOptions,
) -> Result<()>
where
    M: PipelineMode + Clone + Send + 'static,
    M::Metadata: Clone + 'static,
{
    let metadata = Arc::new(mode.validate(general)?);

    let sample_labels = core::derive_sample_labels(&io.scores, general)?;
    let sample_count = sample_labels.len();

    let groups = core::load_groups(&io.regions, gtf)?;
    let mut group_labels: Vec<String> = groups.iter().map(|g| g.label.clone()).collect();
    let mut group_capacity: Vec<usize> = groups.iter().map(|g| g.records.len()).collect();

    // ── Load blacklist & chromosome sizes ───────────────────────────────
    // Keep deepTools compatibility: Python subtracts blacklist intervals
    // from mapReduce genome chunks before region dispatch
    // (deeptools/mapReduce.py:87-104,239-263), then computes signal without
    // a blacklist mask (deeptools/heatmapper.py:531-538). This is not the
    // cleanest design, but output parity depends on it.
    let blacklist: Option<std::collections::HashMap<String, Vec<(u32, u32)>>> =
        if let Some(ref bl_path) = general.blacklist {
            Some(core::regions::load_blacklist(bl_path)?)
        } else {
            None
        };

    let chrom_sizes: Option<std::collections::HashMap<String, u32>> = if blacklist.is_some() {
        Some(core::regions::load_chrom_sizes(&io.scores)?)
    } else {
        None
    };

    let allowed_intervals: Option<std::collections::HashMap<String, Vec<(u32, u32)>>> =
        if let (Some(bl), Some(cs)) = (&blacklist, &chrom_sizes) {
            Some(core::regions::precompute_allowed_intervals(bl, cs))
        } else {
            None
        };

    // ── Generate tasks (with blacklist filtering) ───────────────────────
    let mut tasks = Vec::new();
    for (group_index, group) in groups.into_iter().enumerate() {
        for record in group.records {
            if let Some(ref bl) = blacklist {
                if !core::regions::record_passes_blacklist(
                    &record,
                    bl,
                    allowed_intervals.as_ref().unwrap(),
                    chrom_sizes.as_ref().unwrap(),
                ) {
                    continue;
                }
            }
            let index = tasks.len();
            tasks.push(RegionTask {
                index,
                group_index,
                record: Arc::new(record),
            });
        }
    }

    // ── Validate and remap groups after blacklist filtering ────────────
    if blacklist.is_some() {
        if tasks.is_empty() {
            anyhow::bail!(
                "No regions remain after blacklist filtering. \
                 All {} regions were removed by the blacklist.",
                group_capacity.iter().sum::<usize>()
            );
        }

        let mut post_filter_counts = vec![0usize; group_labels.len()];
        for task in &tasks {
            post_filter_counts[task.group_index] += 1;
        }

        let has_empty_group = post_filter_counts.iter().any(|&c| c == 0);

        if has_empty_group {
            if matches!(general.sort_regions, SortRegions::Keep) {
                // Python errors on empty groups only under --sortRegions keep
                // (computeMatrixOperations.py:729, via sortMatrix).
                let empty = post_filter_counts
                    .iter()
                    .enumerate()
                    .find(|(_, c)| **c == 0)
                    .unwrap()
                    .0;
                anyhow::bail!(
                    "No regions remain in group '{}' after blacklist filtering.",
                    group_labels[empty]
                );
            }

            // For no/ascend/descend, Python drops empty groups from the
            // header entirely (they don't appear in group_labels or
            // group_boundaries). Remap task group_index values and
            // rebuild group_labels/group_capacity to match.
            let mut old_to_new = vec![0usize; group_labels.len()];
            let mut new_labels = Vec::new();
            let mut new_capacity = Vec::new();
            let mut new_idx = 0usize;
            for (old_idx, &count) in post_filter_counts.iter().enumerate() {
                if count > 0 {
                    old_to_new[old_idx] = new_idx;
                    new_labels.push(group_labels[old_idx].clone());
                    new_capacity.push(count);
                    new_idx += 1;
                }
            }
            for task in &mut tasks {
                task.group_index = old_to_new[task.group_index];
            }
            group_labels = new_labels;
            group_capacity = new_capacity;
        }
    }

    let thread_count = std::cmp::max(1, general.number_of_processors.resolve() as usize);
    let task_count = tasks.len();

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
    let collector = FileCollector::new(
        writer,
        &group_labels,
        &header_estimate,
        io.sorted_regions_output.as_deref(),
        io.matrix_values_output.as_deref(),
    )?;

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

    let sample_paths = Arc::new(io.scores.clone());

    core::execute_mode(
        tasks,
        general,
        sample_paths,
        collector,
        thread_count,
        &mode,
        metadata,
        header_builder,
        group_labels.len(),
        task_count,
        sample_count,
    )?;

    Ok(())
}
