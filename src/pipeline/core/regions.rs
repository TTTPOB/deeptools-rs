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

    // Python's heatmapper.py always calls smartLabels() for sample labels
    // regardless of the --smartLabels CLI flag, so we always strip extensions.
    Ok(paths
        .iter()
        .map(|path| label_from_path(path, true))
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

/// Normalize a chromosome name to a canonical form for HashMap indexing.
/// Strips `chr` prefix if present, then maps `M` → `MT` so both `chrM` and `MT` converge.
fn normalize_chrom(name: &str) -> String {
    let stripped = if let Some(rest) = name.strip_prefix("chr") {
        rest
    } else {
        name
    };
    if stripped == "M" {
        "MT".to_string()
    } else {
        stripped.to_string()
    }
}

// ── blacklist helpers ──────────────────────────────────────────────────────

/// Load a blacklist BED file, flatten all groups, and return a HashMap
/// keyed by normalized chromosome name with sorted, merged intervals.
pub(crate) fn load_blacklist(path: &Path) -> Result<HashMap<String, Vec<(u32, u32)>>> {
    let reader = GroupedBedReader::open(path, "blacklist".to_string())?;
    let mut map: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for group in reader {
        let group = group?;
        for record in group.records {
            let key = normalize_chrom(&record.chrom);
            map.entry(key).or_default().push((record.start, record.end));
        }
    }
    // Sort and merge per chromosome
    for intervals in map.values_mut() {
        intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged: Vec<(u32, u32)> = Vec::new();
        for &(start, end) in intervals.iter() {
            if let Some(last) = merged.last_mut() {
                if start <= last.1 {
                    last.1 = last.1.max(end);
                    continue;
                }
            }
            merged.push((start, end));
        }
        *intervals = merged;
    }
    Ok(map)
}

/// Load chromosome sizes from bigWig score files. Chromosome names are
/// normalized via `normalize_chrom()` so lookups are alias-agnostic.
pub(crate) fn load_chrom_sizes(scores: &[PathBuf]) -> Result<HashMap<String, u32>> {
    let mut sizes: HashMap<String, u32> = HashMap::new();
    for path in scores {
        let bw = BigWigFile::open_with_block_cache_capacity(path, 0).map_err(|e| {
            anyhow::anyhow!("Failed to open bigWig file '{}': {}", path.display(), e)
        })?;
        for info in bw.chroms() {
            let key = normalize_chrom(&info.name);
            sizes.entry(key).or_insert(info.length);
        }
    }
    Ok(sizes)
}

/// Return blacklist intervals for `chrom` via normalized HashMap lookup.
/// The returned slice is already sorted and merged.
fn blacklist_intervals_for_chrom<'a>(
    blacklist: &'a HashMap<String, Vec<(u32, u32)>>,
    chrom: &str,
) -> &'a [(u32, u32)] {
    let key = normalize_chrom(chrom);
    blacklist.get(&key).map(|v| v.as_slice()).unwrap_or(&[])
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

/// Precompute allowed (non-blacklisted) intervals per chromosome.
/// This avoids recomputing `subtract_blacklist` for every record.
pub(crate) fn precompute_allowed_intervals(
    blacklist: &HashMap<String, Vec<(u32, u32)>>,
    chrom_sizes: &HashMap<String, u32>,
) -> HashMap<String, Vec<(u32, u32)>> {
    let mut allowed_map = HashMap::new();
    for (chrom, size) in chrom_sizes {
        let bl_intervals = blacklist_intervals_for_chrom(blacklist, chrom);
        if bl_intervals.is_empty() {
            continue;
        }
        let allowed = subtract_blacklist((0, *size), bl_intervals);
        if !allowed.is_empty() {
            allowed_map.insert(chrom.clone(), allowed);
        }
    }
    allowed_map
}

/// Determine whether a BED record should be dispatched given precomputed
/// allowed intervals and the set of chromosomes that have blacklist entries.
///
/// Keep deepTools compatibility: Python subtracts blacklist intervals from
/// mapReduce genome chunks before region dispatch
/// (deeptools/mapReduce.py:87-104,239-263), then computes signal without a
/// blacklist mask (deeptools/heatmapper.py:531-538). This is not the cleanest
/// design, but output parity depends on it.
pub(crate) fn record_passes_blacklist(
    record: &BedRecord,
    blacklist: &HashMap<String, Vec<(u32, u32)>>,
    allowed_intervals: &HashMap<String, Vec<(u32, u32)>>,
    chrom_sizes: &HashMap<String, u32>,
) -> bool {
    let key = normalize_chrom(&record.chrom);
    if !chrom_sizes.contains_key(&key) {
        return true;
    }

    let has_blacklist = blacklist.contains_key(&key);
    if !has_blacklist {
        return true;
    }

    let allowed = match allowed_intervals.get(&key) {
        Some(intervals) => intervals.as_slice(),
        // Chromosome fully covered by blacklist
        None => return false,
    };

    // Python uses findOverlaps(..., trimOverlap=True) which drops any
    // region whose start falls before the allowed chunk start. Match this
    // by requiring record.start to be contained within an allowed interval.
    let i = allowed.partition_point(|&(_, a_end)| a_end <= record.start);
    i < allowed.len() && allowed[i].0 <= record.start
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

    // ── normalize_chrom ────────────────────────────────────────────────────

    #[test]
    fn normalize_chrom_strips_chr_prefix() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("chrX"), "X");
    }

    #[test]
    fn normalize_chrom_no_prefix_unchanged() {
        assert_eq!(normalize_chrom("1"), "1");
        assert_eq!(normalize_chrom("X"), "X");
    }

    #[test]
    fn normalize_chrom_chrm_and_mt_converge() {
        assert_eq!(normalize_chrom("chrM"), "MT");
        assert_eq!(normalize_chrom("MT"), "MT");
    }

    #[test]
    fn normalize_chrom_empty_string() {
        assert_eq!(normalize_chrom(""), "");
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

    fn make_blacklist_map(entries: Vec<(&str, u32, u32)>) -> HashMap<String, Vec<(u32, u32)>> {
        let mut map: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        for (c, s, e) in entries {
            let key = normalize_chrom(c);
            map.entry(key).or_default().push((s, e));
        }
        for intervals in map.values_mut() {
            intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            // merge overlapping/adjacent
            let mut merged: Vec<(u32, u32)> = Vec::new();
            for &(start, end) in intervals.iter() {
                if let Some(last) = merged.last_mut() {
                    if start <= last.1 {
                        last.1 = last.1.max(end);
                        continue;
                    }
                }
                merged.push((start, end));
            }
            *intervals = merged;
        }
        map
    }

    #[test]
    fn blacklist_intervals_single_interval_per_chrom() {
        let bl = make_blacklist_map(vec![("chr1", 10, 20), ("chr2", 30, 40)]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result, &[(10, 20)]);
    }

    #[test]
    fn blacklist_intervals_multiple_intervals_same_chrom() {
        let bl = make_blacklist_map(vec![
            ("chr1", 10, 20),
            ("chr1", 40, 60),
            ("chr1", 90, 95),
            ("chr2", 30, 40),
        ]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result, &[(10, 20), (40, 60), (90, 95)]);
    }

    #[test]
    fn blacklist_intervals_chrom_not_in_blacklist_returns_empty() {
        let bl = make_blacklist_map(vec![("chr1", 10, 20)]);
        let result = blacklist_intervals_for_chrom(&bl, "chrX");
        assert!(result.is_empty());
    }

    #[test]
    fn blacklist_intervals_chr_prefix_toggle_match() {
        let bl = make_blacklist_map(vec![("1", 10, 20)]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        assert_eq!(result, &[(10, 20)]);
    }

    #[test]
    fn blacklist_mixed_chr_alias_merges_correctly() {
        // Regression: entries with mixed `1` and `chr1` must merge into one bucket
        let bl = make_blacklist_map(vec![("1", 100, 200), ("chr1", 150, 250), ("1", 300, 400)]);
        let result = blacklist_intervals_for_chrom(&bl, "chr1");
        // (100,200) and (150,250) overlap → merged to (100,250); (300,400) separate
        assert_eq!(result, &[(100, 250), (300, 400)]);
    }

    #[test]
    fn blacklist_chrm_mt_alias_merges() {
        let bl = make_blacklist_map(vec![("chrM", 100, 200), ("MT", 150, 250)]);
        let result = blacklist_intervals_for_chrom(&bl, "chrM");
        assert_eq!(result, &[(100, 250)]);
        let result2 = blacklist_intervals_for_chrom(&bl, "MT");
        assert_eq!(result2, &[(100, 250)]);
    }

    // ── record_passes_blacklist ─────────────────────────────────────────────

    fn make_record(chrom: &str, start: u32, end: u32) -> BedRecord {
        use crate::io::readers::bed::Strand;
        BedRecord {
            chrom: Arc::from(chrom),
            start,
            end,
            name: None,
            score: None,
            score_raw: None,
            strand: Strand::Unstranded,
            strand_raw: None,
            extra_fields: Vec::new(),
        }
    }

    fn make_chrom_sizes(entries: Vec<(&str, u32)>) -> HashMap<String, u32> {
        entries
            .into_iter()
            .map(|(c, s)| (normalize_chrom(c), s))
            .collect()
    }

    #[test]
    fn passes_blacklist_start_in_allowed_interval() {
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 100, 150);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn fails_blacklist_start_inside_blacklisted_region() {
        // Regression: ch1 115-150 + blacklist ch1 110-130.
        // record.start=115 falls in [110,130) which is blacklisted.
        // Python drops this via trimOverlap=True; Rust must match.
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 115, 150);
        assert!(!record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn fails_blacklist_start_at_blacklist_start_boundary() {
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 110, 150);
        assert!(!record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn passes_blacklist_start_at_blacklist_end_boundary() {
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 130, 150);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn passes_blacklist_no_blacklist_for_chrom() {
        let bl = make_blacklist_map(vec![("ch2", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400), ("ch2", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 115, 150);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn passes_blacklist_chrom_not_in_scores() {
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("chrUn", 0, 100);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn fails_blacklist_fully_covered_chrom() {
        let bl = make_blacklist_map(vec![("ch1", 0, 400)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 100, 200);
        assert!(!record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn passes_blacklist_start_just_before_blacklist() {
        let bl = make_blacklist_map(vec![("ch1", 110, 130)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 109, 150);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    #[test]
    fn passes_blacklist_start_in_gap_between_two_blacklist_intervals() {
        let bl = make_blacklist_map(vec![("ch1", 10, 20), ("ch1", 40, 60)]);
        let cs = make_chrom_sizes(vec![("ch1", 400)]);
        let ai = precompute_allowed_intervals(&bl, &cs);
        let record = make_record("ch1", 30, 50);
        assert!(record_passes_blacklist(&record, &bl, &ai, &cs));
    }

    // ── load_groups cross-file label deduplication ─────────────────────────

    #[test]
    fn load_groups_deduplicates_same_label_across_files() {
        use crate::config::GtfOptions;
        use std::io::Write;
        use tempfile::NamedTempFile;

        // The `# label` line in BED format *terminates* (names) the preceding
        // group of records. Records are written first, then the label line.
        let mut file1 = NamedTempFile::with_suffix(".bed").unwrap();
        writeln!(file1, "chr1\t100\t200").unwrap();
        writeln!(file1, "chr1\t300\t400").unwrap();
        writeln!(file1, "# promoters").unwrap();
        file1.flush().unwrap();

        let mut file2 = NamedTempFile::with_suffix(".bed").unwrap();
        writeln!(file2, "chr2\t500\t600").unwrap();
        writeln!(file2, "# promoters").unwrap();
        file2.flush().unwrap();

        let paths = vec![file1.path().to_path_buf(), file2.path().to_path_buf()];
        let groups = load_groups(&paths, &GtfOptions::default()).unwrap();

        assert_eq!(groups.len(), 2, "expected two groups");
        assert_eq!(
            groups[0].label, "promoters",
            "first group should keep original label"
        );
        assert_eq!(
            groups[1].label, "promoters_1",
            "second group should get _1 suffix"
        );
        assert_eq!(
            groups[0].records.len(),
            2,
            "first group should have 2 records"
        );
        assert_eq!(
            groups[1].records.len(),
            1,
            "second group should have 1 record"
        );
        assert_eq!(groups[0].records[0].chrom.as_ref(), "chr1");
        assert_eq!(groups[0].records[0].start, 100);
        assert_eq!(groups[0].records[0].end, 200);
        assert_eq!(groups[1].records[0].chrom.as_ref(), "chr2");
        assert_eq!(groups[1].records[0].start, 500);
        assert_eq!(groups[1].records[0].end, 600);
    }
}
