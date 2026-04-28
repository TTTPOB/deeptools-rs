use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::config::{GeneralOptions, GtfOptions};
use crate::io::{BedReadError, BedRecord, load_gtf_records};

pub struct Group {
    pub label: String,
    pub records: Vec<BedRecord>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RegionFormat {
    Bed,
    Gtf,
}

#[derive(Clone)]
pub struct RegionTask {
    pub index: usize,
    pub group_index: usize,
    pub record: Arc<BedRecord>,
}

pub fn load_groups(paths: &[PathBuf], gtf: &GtfOptions) -> Result<Vec<Group>> {
    let mut groups = Vec::new();
    let mut seen_labels = HashSet::new();
    // When there's only one file, Python uses "genes" as the default label
    let use_default_genes_label = paths.len() == 1;
    for path in paths {
        match infer_region_format(path) {
            RegionFormat::Bed => {
                let mut file_groups =
                    parse_grouped_bed(path, use_default_genes_label, &mut seen_labels)
                        .map_err(anyhow::Error::new)
                        .with_context(|| {
                            format!("Failed to parse regions file '{}'", path.display())
                        })?;
                groups.append(&mut file_groups);
            }
            RegionFormat::Gtf => {
                let mut file_groups =
                    parse_grouped_gtf(path, gtf, use_default_genes_label, &mut seen_labels)
                        .with_context(|| {
                            format!("Failed to parse regions file '{}'", path.display())
                        })?;
                groups.append(&mut file_groups);
            }
        }
    }

    Ok(groups)
}

pub fn derive_sample_labels(paths: &[PathBuf], general: &GeneralOptions) -> Result<Vec<String>> {
    if let Some(labels) = &general.samples_label {
        if labels.len() != paths.len() {
            bail!(
                "--samplesLabel expects {} entries but {} were provided",
                paths.len(),
                labels.len()
            );
        }
        return Ok(labels.clone());
    }

    Ok(paths
        .iter()
        .map(|path| label_from_path(path, general.smart_labels))
        .collect())
}

pub fn normalize_sort_sample_indices(
    raw: Option<&Vec<usize>>,
    sample_count: usize,
) -> Result<Option<Vec<usize>>> {
    let Some(raw_indices) = raw else {
        return Ok(None);
    };

    if raw_indices.is_empty() {
        return Ok(None);
    }

    let mut normalized = Vec::with_capacity(raw_indices.len());
    for &value in raw_indices {
        if value == 0 || value > sample_count {
            bail!(
                "The value {} for --sortUsingSamples is not valid. Only values from 1 to {} are allowed.",
                value,
                sample_count
            );
        }
        normalized.push(value - 1);
    }

    Ok(Some(normalized))
}

fn parse_grouped_bed(
    path: &Path,
    use_default_genes_label: bool,
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>, BedReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let default_label = if use_default_genes_label {
        "genes".to_string()
    } else {
        bed_file_label(path)
    };
    let mut groups = Vec::new();
    let mut current_records = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            finalize_group(
                trimmed.strip_prefix('#').unwrap_or("").trim(),
                &default_label,
                &mut current_records,
                &mut groups,
                seen_labels,
            );
            continue;
        }

        match BedRecord::parse(trimmed) {
            Ok(record) => current_records.push(record),
            Err(message) => {
                return Err(BedReadError::Parse {
                    line_number: line_number + 1,
                    message,
                    line,
                });
            }
        }
    }

    if !current_records.is_empty() {
        finalize_group(
            "",
            &default_label,
            &mut current_records,
            &mut groups,
            seen_labels,
        );
    } else if groups.is_empty() {
        let label = next_unique_label("", &default_label, seen_labels);
        groups.push(Group {
            label,
            records: Vec::new(),
        });
    }

    Ok(groups)
}

fn parse_grouped_gtf(
    path: &Path,
    options: &GtfOptions,
    use_default_genes_label: bool,
    seen_labels: &mut HashSet<String>,
) -> Result<Vec<Group>> {
    let default_label = if use_default_genes_label {
        "genes".to_string()
    } else {
        bed_file_label(path)
    };
    let mut groups = Vec::new();

    let records = load_gtf_records(path, options)?;
    let label = next_unique_label("", &default_label, seen_labels);
    groups.push(Group { label, records });

    Ok(groups)
}

fn finalize_group(
    raw_label: &str,
    default_label: &str,
    current_records: &mut Vec<BedRecord>,
    groups: &mut Vec<Group>,
    seen_labels: &mut HashSet<String>,
) {
    if current_records.is_empty() {
        return;
    }

    let label = next_unique_label(raw_label, default_label, seen_labels);
    let records = std::mem::take(current_records);
    groups.push(Group { label, records });
}

fn next_unique_label(
    raw_label: &str,
    default_label: &str,
    seen_labels: &mut HashSet<String>,
) -> String {
    let candidate = if raw_label.trim().is_empty() {
        default_label.to_string()
    } else {
        raw_label.trim().to_string()
    };

    if seen_labels.insert(candidate.clone()) {
        return candidate;
    }

    let mut suffix = 1;
    loop {
        let proposal = format!("{}_{}", candidate, suffix);
        if seen_labels.insert(proposal.clone()) {
            return proposal;
        }
        suffix += 1;
    }
}

fn label_from_path(path: &Path, use_stem: bool) -> String {
    if use_stem {
        if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
            return stem.to_string();
        }
    }

    path.file_name()
        .and_then(|name| name.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| path.display().to_string())
}

fn bed_file_label(path: &Path) -> String {
    label_from_path(path, false)
}

fn infer_region_format(path: &Path) -> RegionFormat {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".gtf")
        || lower.ends_with(".gtf.gz")
        || lower.ends_with(".gff")
        || lower.ends_with(".gff.gz")
        || lower.ends_with(".gff3")
        || lower.ends_with(".gff3.gz")
    {
        RegionFormat::Gtf
    } else {
        RegionFormat::Bed
    }
}
