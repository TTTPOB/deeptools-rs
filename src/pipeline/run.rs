use std::sync::Arc;

use anyhow::Result;

use crate::config::{EmptyGroupPolicy, GeneralOptions, GtfOptions, IoOptions};
use crate::io::writers;
use crate::pipeline::core::{self, FileCollector, PipelineMode, RegionTask};
use crate::pipeline::matrix::MatrixHeader;

fn apply_pre_execution_empty_group_policy(
    tasks: &mut [RegionTask],
    group_labels: &mut Vec<String>,
    group_capacity: &mut Vec<usize>,
    sort_regions: crate::config::SortRegions,
    reason: &str,
) -> Result<()> {
    if tasks.is_empty() {
        anyhow::bail!("No regions remain after {reason} filtering.");
    }

    let mut post_filter_counts = vec![0usize; group_labels.len()];
    for task in tasks.iter() {
        post_filter_counts[task.group_index] += 1;
    }

    if !post_filter_counts.iter().any(|&count| count == 0) {
        *group_capacity = post_filter_counts;
        return Ok(());
    }

    match EmptyGroupPolicy::for_sort_regions(sort_regions) {
        EmptyGroupPolicy::Error => {
            let empty = post_filter_counts
                .iter()
                .enumerate()
                .find(|(_, count)| **count == 0)
                .unwrap()
                .0;
            anyhow::bail!(
                "No regions remain in group '{}' after {reason} filtering.",
                group_labels[empty]
            );
        }
        EmptyGroupPolicy::Drop => {
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
            for task in tasks.iter_mut() {
                task.group_index = old_to_new[task.group_index];
            }
            *group_labels = new_labels;
            *group_capacity = new_capacity;
        }
        EmptyGroupPolicy::PreserveWithZeroCount => {
            unreachable!(
                "pre-execution filtering cannot use PreserveWithZeroCount because groups are not yet fixed"
            );
        }
    }

    Ok(())
}

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

    let mut groups = core::load_groups(&io.regions, gtf)?;
    let score_chrom_aliases = core::regions::load_score_chrom_aliases(&io.scores)?;
    let common_score_chroms = core::regions::load_common_score_chroms(&io.scores)?;
    core::regions::remap_group_chroms_to_scores(&mut groups, &score_chrom_aliases);
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

    // ── Generate tasks (with pre-dispatch filters) ──────────────────────
    let mut tasks = Vec::new();
    for (group_index, group) in groups.into_iter().enumerate() {
        for record in group.records {
            if !core::regions::record_chrom_in_common_scores(&record, &common_score_chroms) {
                continue;
            }
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

    apply_pre_execution_empty_group_policy(
        &mut tasks,
        &mut group_labels,
        &mut group_capacity,
        general.sort_regions,
        "pre-dispatch",
    )?;

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

    // Runtime row filtering (--skipZeros, --minThreshold, --maxThreshold)
    // may produce zero-count groups; group structure is already fixed at this
    // point so the only valid policy is PreserveWithZeroCount.
    let runtime_empty_group_policy = EmptyGroupPolicy::PreserveWithZeroCount;

    let header_builder = {
        let general = general.clone();
        let sample_labels = sample_labels.clone();
        let group_labels = group_labels.clone();
        let metadata = Arc::clone(&metadata);
        let mode = mode.clone();
        move |group_counts: Vec<usize>| -> Result<MatrixHeader> {
            if group_counts.iter().all(|count| *count == 0) {
                anyhow::bail!("No regions remain after runtime row filtering.");
            }

            // Explicitly validate against the runtime policy — zero-count
            // groups are expected and allowed when filtering is active.
            runtime_empty_group_policy
                .validate(&group_counts, &group_labels)
                .map_err(|msg| anyhow::anyhow!("{msg}"))?;

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
