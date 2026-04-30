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
    // Ignore "proc number" (thread count varies between machines).
    let status = Command::new(compare_matrix_bin())
        .arg("diff")
        .arg(&output_path)
        .arg(&reference_path)
        .arg("--tolerance")
        .arg(format!("{tolerance}"))
        .arg("--ignore")
        .arg("proc number")
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
