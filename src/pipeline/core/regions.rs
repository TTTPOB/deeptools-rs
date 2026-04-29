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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::path::Path;

    // --- normalize_sort_sample_indices ---

    #[test]
    fn normalize_none_returns_ok_none() {
        let result = normalize_sort_sample_indices(None, 5).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn normalize_empty_vec_returns_ok_none() {
        let indices: Vec<usize> = vec![];
        let result = normalize_sort_sample_indices(Some(&indices), 5).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn normalize_valid_indices_converts_to_zero_based() {
        let indices = vec![1, 3, 5];
        let result = normalize_sort_sample_indices(Some(&indices), 5).unwrap();
        assert_eq!(result, Some(vec![0, 2, 4]));
    }

    #[test]
    fn normalize_index_zero_returns_error() {
        let indices = vec![0];
        let result = normalize_sort_sample_indices(Some(&indices), 5);
        assert!(result.is_err());
    }

    #[test]
    fn normalize_index_exceeds_sample_count_returns_error() {
        let indices = vec![6];
        let result = normalize_sort_sample_indices(Some(&indices), 5);
        assert!(result.is_err());
    }

    // --- label_from_path ---

    #[test]
    fn label_from_path_use_stem_false_returns_filename_with_extension() {
        let path = Path::new("/some/dir/sample.bw");
        let label = label_from_path(path, false);
        assert_eq!(label, "sample.bw");
    }

    #[test]
    fn label_from_path_use_stem_true_returns_filename_without_extension() {
        let path = Path::new("/some/dir/sample.bw");
        let label = label_from_path(path, true);
        assert_eq!(label, "sample");
    }

    // --- infer_region_format ---

    #[test]
    fn infer_format_bed_extension() {
        assert_eq!(infer_region_format(Path::new("regions.bed")), RegionFormat::Bed);
    }

    #[test]
    fn infer_format_gtf_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gtf")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_gtf_gz_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gtf.gz")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_gff_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gff")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_gff_gz_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gff.gz")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_gff3_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gff3")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_gff3_gz_extension() {
        assert_eq!(infer_region_format(Path::new("genes.gff3.gz")), RegionFormat::Gtf);
    }

    #[test]
    fn infer_format_txt_defaults_to_bed() {
        assert_eq!(infer_region_format(Path::new("regions.txt")), RegionFormat::Bed);
    }

    #[test]
    fn infer_format_case_insensitive_bed() {
        assert_eq!(infer_region_format(Path::new("regions.BED")), RegionFormat::Bed);
    }

    // --- next_unique_label ---

    #[test]
    fn next_unique_label_empty_raw_uses_default() {
        let mut seen = HashSet::new();
        let label = next_unique_label("", "genes", &mut seen);
        assert_eq!(label, "genes");
    }

    #[test]
    fn next_unique_label_nonempty_raw_uses_raw() {
        let mut seen = HashSet::new();
        let label = next_unique_label("promoters", "genes", &mut seen);
        assert_eq!(label, "promoters");
    }

    #[test]
    fn next_unique_label_duplicate_appends_suffix() {
        let mut seen = HashSet::new();
        let first = next_unique_label("genes", "default", &mut seen);
        let second = next_unique_label("genes", "default", &mut seen);
        let third = next_unique_label("genes", "default", &mut seen);
        assert_eq!(first, "genes");
        assert_eq!(second, "genes_1");
        assert_eq!(third, "genes_2");
    }

    #[test]
    fn next_unique_label_whitespace_only_raw_uses_default() {
        let mut seen = HashSet::new();
        let label = next_unique_label("   ", "genes", &mut seen);
        assert_eq!(label, "genes");
    }
}
