use std::path::PathBuf;
use std::process::Command;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn compute_matrix_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_compute_matrix_rs"))
}

fn compare_matrix_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_compare_matrix"))
}

fn data_root() -> PathBuf {
    project_root().join("deeptools/deeptools/test/test_heatmapper")
}

fn test_data_root() -> PathBuf {
    project_root().join("deeptools/deeptools/test/test_data")
}

fn blacklist_data_root() -> PathBuf {
    project_root().join("tests/data")
}

fn corner_case_root() -> PathBuf {
    project_root().join("tests/data/corner_cases")
}

fn run_compute_and_compare_corner(reference_mat: &str, args: &[&str], tolerance: f64) {
    let reference_path = corner_case_root().join(reference_mat);
    run_compute_and_compare_at(reference_path, args, tolerance);
}

/// Run compute_matrix_rs with given args, write output to a temp file,
/// then compare the output against `reference_mat` using compare_matrix diff.
///
/// `reference_mat` is resolved relative to `data_root()`.
fn run_compute_and_compare(reference_mat: &str, args: &[&str], tolerance: f64) {
    let reference_path = data_root().join(reference_mat);
    run_compute_and_compare_at(reference_path, args, tolerance);
}

/// Like `run_compute_and_compare`, but `reference_mat` is resolved relative to `blacklist_data_root()`.
fn run_compute_and_compare_blacklist(reference_mat: &str, args: &[&str], tolerance: f64) {
    let reference_path = blacklist_data_root().join(reference_mat);
    run_compute_and_compare_at(reference_path, args, tolerance);
}

fn run_compute_and_compare_at(reference_path: PathBuf, args: &[&str], tolerance: f64) {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let output_path = tmp.path().to_path_buf();

    // Run compute_matrix_rs
    let mut cmd = Command::new(compute_matrix_bin());
    cmd.args(args);
    cmd.arg("-o").arg(&output_path);
    let status = cmd.status().expect("failed to run compute_matrix_rs");
    assert!(status.success(), "compute_matrix_rs failed with {status}");

    // Run compare_matrix diff
    // Ignore "proc number" (thread count varies between machines) and "scale"
    // (Python emits int for the default value vs our float — the difference is
    // harmless because no downstream tool reads the scale field back).
    let status = Command::new(compare_matrix_bin())
        .arg("diff")
        .arg(&output_path)
        .arg(&reference_path)
        .arg("--tolerance")
        .arg(format!("{tolerance}"))
        .arg("--ignore")
        .arg("proc number")
        .arg("--ignore")
        .arg("scale")
        .status()
        .expect("failed to run compare_matrix");

    assert!(
        status.success(),
        "Matrix mismatch for {}! compare_matrix diff exited with {status}",
        reference_path.display()
    );
}

#[test]
fn reference_point_basic() {
    let dr = data_root();
    run_compute_and_compare(
        "master.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn reference_point_center() {
    let dr = data_root();
    run_compute_and_compare(
        "master_center.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--referencePoint",
            "center",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn reference_point_tes() {
    let dr = data_root();
    run_compute_and_compare(
        "master_TES.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--referencePoint",
            "TES",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn reference_point_missing_data_as_zero() {
    let dr = data_root();
    run_compute_and_compare(
        "master_nan_to_zero.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--missingDataAsZero",
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_basic() {
    let dr = data_root();
    run_compute_and_compare(
        "master_scale_reg.mat",
        &[
            "scale-regions",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "-m",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn multiple_bed() {
    let dr = data_root();
    run_compute_and_compare(
        "master_multibed.mat",
        &[
            "reference-point",
            "-R",
            dr.join("group1.bed").to_str().unwrap(),
            dr.join("group2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn region_extend_beyond_chr() {
    let dr = data_root();
    run_compute_and_compare(
        "master_extend_beyond_chr_size.mat",
        &[
            "reference-point",
            "-R",
            dr.join("group1.bed").to_str().unwrap(),
            dr.join("group2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "500",
            "--bs",
            "1",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_unscaled() {
    let dr = data_root();
    run_compute_and_compare(
        "master_unscaled.mat",
        &[
            "scale-regions",
            "-R",
            dr.join("unscaled.bed").to_str().unwrap(),
            "-S",
            dr.join("unscaled.bigWig").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "100",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn gtf_input() {
    let tdr = test_data_root();
    run_compute_and_compare(
        "master_gtf.mat",
        &[
            "scale-regions",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "20",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
        ],
        5e-6,
    );
}

// ── blacklist parity tests ──────────────────────────────────────────────────

#[test]
fn reference_point_blacklist() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_blacklist.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn reference_point_blacklist_missing_data_as_zero() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_blacklist_nan_to_zero.mat",
        &[
            "reference-point",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--missingDataAsZero",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_blacklist() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_scale_reg_blacklist.mat",
        &[
            "scale-regions",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "-m",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_blacklist_missing_data_as_zero() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_scale_reg_blacklist_nan_to_zero.mat",
        &[
            "scale-regions",
            "-R",
            dr.join("test2.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "-m",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--missingDataAsZero",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

// ── blacklist + empty group parity (non-keep sortRegions) ──────────────────
// When blacklist empties a group, Python drops it from the header for
// no/ascend/descend modes. These tests verify the group remap produces
// matching headers and values.

#[test]
fn blacklist_empty_group_sort_no() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_blacklist_empty_group_no.mat",
        &[
            "reference-point",
            "-R",
            bdr.join("test_blacklist_empty_group.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--sortRegions",
            "no",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn blacklist_empty_group_sort_descend() {
    let dr = data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_blacklist_empty_group_descend.mat",
        &[
            "reference-point",
            "-R",
            bdr.join("test_blacklist_empty_group.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "1",
            "-p",
            "1",
            "--sortRegions",
            "descend",
            "--blackListFileName",
            bdr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

// ── metagene + intron blacklist parity ──────────────────────────────────────
// Spec requires: "scale-regions, BED12/GTF metagene: Blacklist falls in an
// intron of a multi-exon region." The blacklist at 3R:200-300 falls inside
// introns of both transcripts (exon gaps 50-400 and 150-500), so neither
// metagene record is filtered.

#[test]
fn metagene_blacklist_intron() {
    let tdr = test_data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_metagene_blacklist_intron.mat",
        &[
            "scale-regions",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "20",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
            "--blackListFileName",
            bdr.join("test_blacklist_intron.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn metagene_blacklist_intron_missing_data_as_zero() {
    let tdr = test_data_root();
    let bdr = blacklist_data_root();
    run_compute_and_compare_blacklist(
        "master_metagene_blacklist_intron_nan_to_zero.mat",
        &[
            "scale-regions",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "20",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
            "--missingDataAsZero",
            "--blackListFileName",
            bdr.join("test_blacklist_intron.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn metagene() {
    let tdr = test_data_root();
    run_compute_and_compare(
        "master_metagene.mat",
        &[
            "scale-regions",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "20",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
        ],
        5e-6,
    );
}

#[test]
fn metagene_missing_data_as_zero() {
    let tdr = test_data_root();
    run_compute_and_compare_blacklist(
        "master_metagene_nan_to_zero.mat",
        &[
            "scale-regions",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "-a",
            "300",
            "-b",
            "500",
            "--unscaled5prime",
            "20",
            "--unscaled3prime",
            "50",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
            "--missingDataAsZero",
        ],
        5e-6,
    );
}

// ── metagene reference-point parity tests ─────────────────────────────────

#[test]
fn corner_case_metagene_reference_point() {
    let tdr = test_data_root();
    // With included_intervals now filtering out intronic signal, Rust matches
    // Python's exon-fragment coverage fetching to within floating-point tolerance.
    run_compute_and_compare_corner(
        "master_metagene_refpoint.mat",
        &[
            "reference-point",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "--referencePoint",
            "TSS",
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
        ],
        5e-6,
    );
}

#[test]
fn corner_case_metagene_reference_point_center() {
    let tdr = test_data_root();
    // With included_intervals now filtering out intronic signal, Rust matches
    // Python's exon-fragment coverage fetching to within floating-point tolerance.
    run_compute_and_compare_corner(
        "master_metagene_center.mat",
        &[
            "reference-point",
            "-R",
            tdr.join("test.gtf").to_str().unwrap(),
            "-S",
            tdr.join("test1.bw.bw").to_str().unwrap(),
            "--referencePoint",
            "center",
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--metagene",
        ],
        5e-6,
    );
}

// ── Corner case integration tests ─────────────────────────────────────────

#[test]
fn corner_case_short_body_scale_regions() {
    let dr = data_root();
    run_compute_and_compare_corner(
        "master_short_body.mat",
        &[
            "scale-regions",
            "-R",
            corner_case_root().join("short_body.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-m",
            "100",
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
        ],
        5e-6,
    );
}

#[test]
fn corner_case_short_body_missing_data_as_zero() {
    let dr = data_root();
    run_compute_and_compare_corner(
        "master_short_body_nan_to_zero.mat",
        &[
            "scale-regions",
            "-R",
            corner_case_root().join("short_body.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-m",
            "100",
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--missingDataAsZero",
        ],
        5e-6,
    );
}

#[test]
fn corner_case_scale_with_max_threshold() {
    let dr = data_root();
    // Python checks thresholds pre-scale: raw max 3.0 < 4 → keep all rows.
    // Old Rust checked post-scale: scaled max 6.0 >= 4 → would drop rows.
    run_compute_and_compare_corner(
        "master_scale_threshold.mat",
        &[
            "reference-point",
            "-R",
            dr.join("group1.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--scale",
            "2",
            "--maxThreshold",
            "4",
        ],
        5e-6,
    );
}

#[test]
fn corner_case_scale_zero_with_max_threshold() {
    let dr = data_root();
    // --scale 0 --maxThreshold 3: ch3 raw max 3.0 >= 3 → filtered out,
    // ch1/ch2 kept with all values scaled to 0.  Old Rust scaled first
    // (all values became 0), then threshold check saw 0 < 3 → kept all rows.
    run_compute_and_compare_corner(
        "master_scale_zero_threshold.mat",
        &[
            "reference-point",
            "-R",
            dr.join("group1.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--scale",
            "0",
            "--maxThreshold",
            "3",
        ],
        5e-6,
    );
}

#[test]
fn corner_case_skipzeros_removes_zero_region() {
    let dr = data_root();
    run_compute_and_compare_corner(
        "master_skipzeros.mat",
        &[
            "reference-point",
            "-R",
            corner_case_root().join("skipzeros.bed").to_str().unwrap(),
            "-S",
            dr.join("test.bw").to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--bs",
            "10",
            "-p",
            "1",
            "--missingDataAsZero",
            "--skipZeros",
        ],
        5e-6,
    );
}

#[test]
fn outfilename_matrix_header_consistent_after_skipzeros() {
    let dr = data_root();
    let mat_tmp = tempfile::NamedTempFile::new().unwrap();
    let tab_tmp = tempfile::NamedTempFile::new().unwrap();

    let mut cmd = Command::new(compute_matrix_bin());
    cmd.args([
        "reference-point",
        "-R",
        corner_case_root().join("skipzeros.bed").to_str().unwrap(),
        "-S",
        dr.join("test.bw").to_str().unwrap(),
        "-b",
        "100",
        "-a",
        "100",
        "--bs",
        "10",
        "-p",
        "1",
        "--missingDataAsZero",
        "--skipZeros",
        "-o",
        mat_tmp.path().to_str().unwrap(),
        "--outFileNameMatrix",
        tab_tmp.path().to_str().unwrap(),
    ]);
    let status = cmd.status().unwrap();
    assert!(status.success());

    let tab = std::fs::read_to_string(tab_tmp.path()).unwrap();
    let lines: Vec<&str> = tab.lines().collect();
    assert!(
        lines.len() >= 4,
        "expected >=3 header lines + 1 data row, got {}",
        lines.len()
    );

    // Parse group counts from line 1: "#Group1:N\tGroup2:M"
    let line1_stripped = lines[0].strip_prefix('#').unwrap_or_else(|| {
        panic!("line 1 should start with '#': {}", lines[0]);
    });
    let line1_counts: Vec<usize> = line1_stripped
        .split('\t')
        .map(|p| p.split(':').nth(1).unwrap().parse::<usize>().unwrap())
        .collect();
    let line1_total: usize = line1_counts.iter().sum();

    // Parse group counts from line 3 first N entries
    let line3_parts: Vec<&str> = lines[2].split('\t').collect();
    let num_groups = line1_counts.len();
    let line3_counts: Vec<usize> = line3_parts[..num_groups]
        .iter()
        .map(|p| p.split(':').nth(1).unwrap().parse::<usize>().unwrap())
        .collect();

    let data_rows = lines.len() - 3; // 3 header lines

    assert_eq!(
        line1_counts, line3_counts,
        "line 1 and line 3 group counts diverge in:\n{tab}"
    );
    assert_eq!(
        line1_total, data_rows,
        "header total {line1_total} != data rows {data_rows} in:\n{tab}"
    );
}

/// Runtime empty group: two groups, one entirely filtered by --skipZeros.
/// Verifies that the zero-count group appears explicitly in both the gzip
/// header (`group_boundaries`) and the `--outFileNameMatrix` tab header.
#[test]
fn runtime_empty_group_preserved_with_zero_count() {
    use std::io::Write;

    let dr = data_root();

    // Group "genes": ch1:100-150 — has signal in test.bw, survives skipZeros.
    let mut genes_bed = tempfile::NamedTempFile::new().unwrap();
    writeln!(genes_bed, "ch1\t100\t150\thas_signal\t0\t+").unwrap();

    // Group "peaks": ch1:300-350 — no signal in test.bw, filtered by skipZeros + missingDataAsZero.
    let mut peaks_bed = tempfile::NamedTempFile::new().unwrap();
    writeln!(peaks_bed, "ch1\t300\t350\tno_signal\t0\t+").unwrap();

    let mat_tmp = tempfile::NamedTempFile::new().unwrap();
    let tab_tmp = tempfile::NamedTempFile::new().unwrap();

    let mut cmd = Command::new(compute_matrix_bin());
    cmd.args([
        "reference-point",
        "-R",
        genes_bed.path().to_str().unwrap(),
        peaks_bed.path().to_str().unwrap(),
        "--smartLabels",
        "-S",
        dr.join("test.bw").to_str().unwrap(),
        "-b",
        "100",
        "-a",
        "100",
        "--bs",
        "10",
        "-p",
        "1",
        "--missingDataAsZero",
        "--skipZeros",
        "-o",
        mat_tmp.path().to_str().unwrap(),
        "--outFileNameMatrix",
        tab_tmp.path().to_str().unwrap(),
    ]);
    let status = cmd.status().unwrap();
    assert!(status.success());

    // ── Verify tab file ──────────────────────────────────────────────
    let tab = std::fs::read_to_string(tab_tmp.path()).unwrap();
    let lines: Vec<&str> = tab.lines().collect();
    assert!(
        lines.len() >= 4,
        "expected >=3 header lines + 1 data row, got {}",
        lines.len()
    );

    // Line 1: #genes:1\tpeaks:0
    let line1_stripped = lines[0]
        .strip_prefix('#')
        .unwrap_or_else(|| panic!("line 1 should start with '#': {}", lines[0]));
    let line1_parts: Vec<&str> = line1_stripped.split('\t').collect();
    assert_eq!(line1_parts.len(), 2);
    let genes_count: usize = line1_parts[0].split(':').nth(1).unwrap().parse().unwrap();
    let peaks_count: usize = line1_parts[1].split(':').nth(1).unwrap().parse().unwrap();
    assert_eq!(genes_count, 1, "genes group should have 1 surviving row");
    assert_eq!(
        peaks_count, 0,
        "peaks group should be empty after skipZeros"
    );

    // Line 3: genes:1\tpeaks:0\ttest\ttest\t...
    let line3_parts: Vec<&str> = lines[2].split('\t').collect();
    let line3_genes: usize = line3_parts[0].split(':').nth(1).unwrap().parse().unwrap();
    let line3_peaks: usize = line3_parts[1].split(':').nth(1).unwrap().parse().unwrap();
    assert_eq!(line3_genes, 1, "line 3 genes count should match line 1");
    assert_eq!(line3_peaks, 0, "line 3 peaks count should match line 1");

    // Data rows: only the genes row
    let data_rows = lines.len() - 3;
    assert_eq!(
        data_rows, 1,
        "should have exactly 1 data row, got {data_rows}"
    );

    // ── Verify gzip header ───────────────────────────────────────────
    use std::io::BufRead;
    let mat_bytes = std::fs::read(mat_tmp.path()).unwrap();
    let decoder = flate2::read::GzDecoder::new(&mat_bytes[..]);
    let mut reader = std::io::BufReader::new(decoder);
    let mut header_line = String::new();
    reader.read_line(&mut header_line).unwrap();
    // Strip the leading '@' and trailing newline (plus padding spaces).
    let json = header_line
        .strip_prefix('@')
        .unwrap()
        .trim_end_matches(|c: char| c == '\n' || c == ' ');
    let header: serde_json::Value = serde_json::from_str(json).unwrap();
    let boundaries: Vec<usize> = header["group_boundaries"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    // genes:1, peaks:0 → boundaries should be [0, 1, 1]
    assert_eq!(boundaries, vec![0, 1, 1]);

    let labels = header["group_labels"].as_array().unwrap();
    assert_eq!(labels.len(), 2, "should have 2 group labels");
    // With --smartLabels, labels are derived from BED file names.  Just
    // verify they're non-empty — exact names depend on temp file paths.
    assert!(!labels[0].as_str().unwrap().is_empty());
    assert!(!labels[1].as_str().unwrap().is_empty());
}
