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

/// Run compute_matrix_rs with given args, write output to a temp file,
/// then compare the output against `reference_mat` using compare_matrix diff.
///
/// `reference_mat` is resolved relative to `data_root()`.
fn run_compute_and_compare(reference_mat: &str, args: &[&str], tolerance: f64) {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let output_path = tmp.path().to_path_buf();

    // Run compute_matrix_rs
    let mut cmd = Command::new(compute_matrix_bin());
    cmd.args(args);
    cmd.arg("-o").arg(&output_path);
    let status = cmd.status().expect("failed to run compute_matrix_rs");
    assert!(status.success(), "compute_matrix_rs failed with {status}");

    // Run compare_matrix diff
    let reference_path = data_root().join(reference_mat);
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
        "Matrix mismatch for {reference_mat}! compare_matrix diff exited with {status}"
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
    run_compute_and_compare(
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
            dr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn reference_point_blacklist_missing_data_as_zero() {
    let dr = data_root();
    run_compute_and_compare(
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
            dr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_blacklist() {
    let dr = data_root();
    run_compute_and_compare(
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
            dr.join("test_blacklist.bed").to_str().unwrap(),
        ],
        5e-6,
    );
}

#[test]
fn scale_regions_blacklist_missing_data_as_zero() {
    let dr = data_root();
    run_compute_and_compare(
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
            dr.join("test_blacklist.bed").to_str().unwrap(),
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
