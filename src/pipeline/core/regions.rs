use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::config::{GeneralOptions, GtfOptions};
use crate::io::{BedReadError, BedRecord, BigWigFile, Group, GroupedBedReader, load_gtf_records};

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
        let default_label = if use_default_genes_label {
            "genes".to_string()
        } else {
            bed_file_label(path)
        };
        match infer_region_format(path) {
            RegionFormat::Bed => {
                let file_groups = parse_grouped_bed(path, default_label.clone())
                    .map_err(anyhow::Error::new)
                    .with_context(|| {
                        format!("Failed to parse regions file '{}'", path.display())
                    })?;
                for mut group in file_groups {
                    group.label = next_unique_label(&group.label, &default_label, &mut seen_labels);
                    groups.push(group);
                }
            }
            RegionFormat::Gtf => {
                let file_groups = parse_grouped_gtf(path, gtf, default_label.clone())
                    .with_context(|| {
                        format!("Failed to parse regions file '{}'", path.display())
                    })?;
                for mut group in file_groups {
                    group.label = next_unique_label(&group.label, &default_label, &mut seen_labels);
                    groups.push(group);
                }
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

// ── chromosome name normalization ──────────────────────────────────────────

/// Return both the original name and a chr-prefix-toggled variant so that
/// names like `chr1` and `1` can be matched against each other.  The first
/// element is the raw name; the second is the toggled version.
fn normalize_chrom_name(name: &str) -> [String; 2] {
    if name.starts_with("chr") {
        [name.to_string(), name[3..].to_string()]
    } else {
        [name.to_string(), format!("chr{}", name)]
    }
}

// ── blacklist helpers ──────────────────────────────────────────────────────

/// Load a blacklist BED file, flatten all groups, sort by (chrom, start), and
/// merge overlapping/adjacent intervals.
pub(crate) fn load_blacklist(path: &Path) -> Result<Vec<(Arc<str>, u32, u32)>> {
    let reader = GroupedBedReader::open(path, "blacklist".to_string())?;
    let mut intervals: Vec<(Arc<str>, u32, u32)> = Vec::new();
    for group in reader {
        let group = group?;
        for record in group.records {
            intervals.push((record.chrom, record.start, record.end));
        }
    }
    if intervals.is_empty() {
        return Ok(Vec::new());
    }
    intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut merged: Vec<(Arc<str>, u32, u32)> = Vec::new();
    for (chrom, start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            if last.0 == chrom && start <= last.2 {
                last.2 = last.2.max(end);
                continue;
            }
        }
        merged.push((chrom, start, end));
    }
    Ok(merged)
}

/// Load chromosome sizes from bigWig score files. Chromosome names are
/// stored as-is (original form from the bigWig header).
pub(crate) fn load_chrom_sizes(scores: &[PathBuf]) -> Result<HashMap<String, u32>> {
    let mut sizes: HashMap<String, u32> = HashMap::new();
    for path in scores {
        let bw = BigWigFile::open_with_block_cache_capacity(path, 0).map_err(|e| {
            anyhow::anyhow!("Failed to open bigWig file '{}': {}", path.display(), e)
        })?;
        for info in bw.chroms() {
            let [canonical_name, _] = normalize_chrom_name(&info.name);
            sizes.entry(canonical_name).or_insert(info.length);
        }
    }
    Ok(sizes)
}

/// Return blacklist intervals for `chrom`, matching both original and
/// chr-prefix-toggled chromosome names. Uses partition_point for true
/// lower/upper bound search in O(log n + k) time, avoiding the
/// binary_search_by pitfall where the match position is not guaranteed
/// to be the first equal element.
fn blacklist_intervals_for_chrom<'a>(
    blacklist: &'a [(Arc<str>, u32, u32)],
    chrom: &str,
) -> &'a [(Arc<str>, u32, u32)] {
    let [ref name_a, name_b] = normalize_chrom_name(chrom);

    // find the first entry >= the lexicographically smaller variant
    let candidate_a = name_a.as_str();
    let candidate_b = name_b.as_str();
    let probe = candidate_a.min(candidate_b);

    let lo = blacklist.partition_point(|(c, _, _)| c.as_ref() < probe);

    // scan forward to find the range matching either variant
    let mut hi = lo;
    while hi < blacklist.len()
        && (blacklist[hi].0.as_ref() == name_a || blacklist[hi].0.as_ref() == name_b)
    {
        hi += 1;
    }
    &blacklist[lo..hi]
}

/// Subtract sorted, non-overlapping blacklist intervals from a genomic span.
/// Returns the resulting allowed intervals.
///
/// Edge cases (matching Python `blSubtract()`):
/// - No overlap: blacklist intervals before or after the span leave it intact.
/// - Partial overlap: start/end portions are trimmed.
/// - Enclosed: blacklist interval inside the span splits it in two.
/// - Adjacent: back-to-back blacklist intervals are handled correctly.
/// - Empty blacklist: the full span is returned.
pub fn subtract_blacklist(span: (u32, u32), blacklist: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut result = Vec::new();
    let mut cursor = span.0;
    for &(bl_start, bl_end) in blacklist {
        if bl_end <= cursor {
            continue;
        }
        if bl_start > cursor {
            result.push((cursor, bl_start.min(span.1)));
        }
        cursor = bl_end;
        if cursor >= span.1 {
            break;
        }
    }
    if cursor < span.1 {
        result.push((cursor, span.1));
    }
    result
}

/// Determine whether a BED record should be dispatched given the blacklist.
///
/// Keep deepTools compatibility: Python subtracts blacklist intervals from
/// mapReduce genome chunks before region dispatch
/// (deeptools/mapReduce.py:87-104,239-263), then computes signal without a
/// blacklist mask (deeptools/heatmapper.py:531-538). This is not the cleanest
/// design, but output parity depends on it.
pub(crate) fn record_passes_blacklist(
    record: &BedRecord,
    blacklist: &[(Arc<str>, u32, u32)],
    chrom_sizes: &HashMap<String, u32>,
) -> bool {
    let [canonical_name, _] = normalize_chrom_name(&record.chrom);
    let chrom_size = match chrom_sizes.get(&canonical_name) {
        Some(s) => *s,
        // Chromosome not in score files: allow dispatch (it will produce NaN
        // rows downstream, matching Python behavior).
        None => return true,
    };

    let bl_slice = blacklist_intervals_for_chrom(blacklist, &canonical_name);
    if bl_slice.is_empty() {
        return true;
    }

    let bl_intervals: Vec<(u32, u32)> = bl_slice.iter().map(|(_, s, e)| (*s, *e)).collect();
    let allowed = subtract_blacklist((0, chrom_size), &bl_intervals);

    // A record is dispatched if its interval overlaps any allowed interval.
    allowed
        .iter()
        .any(|&(a_start, a_end)| record.start < a_end && a_start < record.end)
}

fn parse_grouped_bed(path: &Path, default_label: String) -> Result<Vec<Group>, BedReadError> {
    let reader = GroupedBedReader::open(path, default_label)?;
    reader.collect()
}

fn parse_grouped_gtf(
    path: &Path,
    options: &GtfOptions,
    default_label: String,
) -> Result<Vec<Group>> {
    let records = load_gtf_records(path, options)?;
    if records.is_empty() {
        bail!("no data records found in GTF file '{}'", path.display());
    }
    // Raw label (undeduplicated); load_groups() will deduplicate.
    Ok(vec![Group {
        label: default_label,
        records,
    }])
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
        assert_eq!(
            infer_region_format(Path::new("regions.bed")),
            RegionFormat::Bed
        );
    }

    #[test]
    fn infer_format_gtf_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gtf")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_gtf_gz_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gtf.gz")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_gff_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gff")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_gff_gz_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gff.gz")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_gff3_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gff3")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_gff3_gz_extension() {
        assert_eq!(
            infer_region_format(Path::new("genes.gff3.gz")),
            RegionFormat::Gtf
        );
    }

    #[test]
    fn infer_format_txt_defaults_to_bed() {
        assert_eq!(
            infer_region_format(Path::new("regions.txt")),
            RegionFormat::Bed
        );
    }

    #[test]
    fn infer_format_case_insensitive_bed() {
        assert_eq!(
            infer_region_format(Path::new("regions.BED")),
            RegionFormat::Bed
        );
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

    // ── normalize_chrom_name ───────────────────────────────────────────────

    #[test]
    fn normalize_chrom_with_prefix() {
        let [a, b] = normalize_chrom_name("chr1");
        assert_eq!(a, "chr1");
        assert_eq!(b, "1");
    }

    #[test]
    fn normalize_chrom_without_prefix() {
        let [a, b] = normalize_chrom_name("1");
        assert_eq!(a, "1");
        assert_eq!(b, "chr1");
    }

    #[test]
    fn normalize_chrom_empty() {
        let [a, b] = normalize_chrom_name("");
        assert_eq!(a, "");
        assert_eq!(b, "chr");
    }

    // ── subtract_blacklist ─────────────────────────────────────────────────

    #[test]
    fn subtract_blacklist_no_overlap_before() {
        // blacklist entirely before span
        let result = subtract_blacklist((20, 30), &[(0, 10)]);
        assert_eq!(result, vec![(20, 30)]);
    }

    #[test]
    fn subtract_blacklist_no_overlap_after() {
        // blacklist entirely after span
        let result = subtract_blacklist((0, 10), &[(20, 30)]);
        assert_eq!(result, vec![(0, 10)]);
    }

    #[test]
    fn subtract_blacklist_partial_overlap_left() {
        // blacklist overlaps start of span
        let result = subtract_blacklist((10, 30), &[(5, 15)]);
        assert_eq!(result, vec![(15, 30)]);
    }

    #[test]
    fn subtract_blacklist_partial_overlap_right() {
        // blacklist overlaps end of span
        let result = subtract_blacklist((10, 30), &[(25, 35)]);
        assert_eq!(result, vec![(10, 25)]);
    }

    #[test]
    fn subtract_blacklist_enclosed() {
        // blacklist fully inside span → splits span in two
        let result = subtract_blacklist((10, 30), &[(15, 20)]);
        assert_eq!(result, vec![(10, 15), (20, 30)]);
    }

    #[test]
    fn subtract_blacklist_adjacent() {
        let result = subtract_blacklist((0, 30), &[(0, 10), (10, 20)]);
        assert_eq!(result, vec![(20, 30)]);
    }

    #[test]
    fn subtract_blacklist_empty() {
        let result = subtract_blacklist((0, 10), &[]);
        assert_eq!(result, vec![(0, 10)]);
    }

    #[test]
    fn subtract_blacklist_span_fully_covered() {
        // blacklist covers entire span
        let result = subtract_blacklist((10, 20), &[(0, 30)]);
        assert!(result.is_empty());
    }

    #[test]
    fn subtract_blacklist_multiple_intervals() {
        let result = subtract_blacklist((0, 100), &[(10, 20), (40, 60), (90, 95)]);
        assert_eq!(result, vec![(0, 10), (20, 40), (60, 90), (95, 100)]);
    }

    // ── blacklist_intervals_for_chrom ────────────────────────────────────────

    fn make_blacklist(entries: Vec<(&str, u32, u32)>) -> Vec<(Arc<str>, u32, u32)> {
        entries
            .into_iter()
            .map(|(c, s, e)| (Arc::from(c), s, e))
            .collect()
    }

    #[test]
    fn blacklist_intervals_single_interval_per_chrom() {
        let bl = make_blacklist(vec![("chr1", 10, 20), ("chr2", 30, 40)]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, 10);
        assert_eq!(result[0].2, 20);
    }

    #[test]
    fn blacklist_intervals_multiple_intervals_same_chrom() {
        // Regression: binary_search_by does not guarantee the first match.
        // With 3 intervals on chr1, the search must return all 3, not a
        // suffix starting from a later interval.
        let bl = make_blacklist(vec![
            ("chr1", 10, 20),
            ("chr1", 40, 60),
            ("chr1", 90, 95),
            ("chr2", 30, 40),
        ]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result.len(), 3);
        let starts: Vec<u32> = result.iter().map(|(_, s, _)| *s).collect();
        assert_eq!(starts, vec![10, 40, 90]);
    }

    #[test]
    fn blacklist_intervals_chrom_not_in_blacklist_returns_empty() {
        let bl = make_blacklist(vec![("chr1", 10, 20)]);
        let result = blacklist_intervals_for_chrom(&bl, "chrX");
        assert!(result.is_empty());
    }

    #[test]
    fn blacklist_intervals_chr_prefix_toggle_match() {
        let bl = make_blacklist(vec![("1", 10, 20)]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, 10);
    }
}
