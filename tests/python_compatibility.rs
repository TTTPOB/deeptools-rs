use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use compute_matrix_rs::matrix_compare::diff::full_diff;
use compute_matrix_rs::matrix_compare::parse::load_matrix;
use serde::Deserialize;
use tempfile::NamedTempFile;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn compute_matrix_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_compute_matrix_rs"))
}

#[derive(Debug, Deserialize)]
struct CaseManifest {
    paths: HashMap<String, String>,
    comparison: ComparisonConfig,
    cases: Vec<CaseConfig>,
}

#[derive(Debug, Deserialize)]
struct ComparisonConfig {
    default_tolerance: f64,
    ignore_header_keys: Vec<String>,
    max_diffs: usize,
}

#[derive(Debug, Deserialize)]
struct CaseConfig {
    id: String,
    tags: Vec<String>,
    mode: String,
    region_files: Vec<String>,
    score_files: Vec<String>,
    options: Vec<String>,
    reference: Option<ReferenceConfig>,
}

#[derive(Debug, Deserialize)]
struct ReferenceConfig {
    matrix: String,
}

fn load_manifest() -> CaseManifest {
    let config_dir = project_root().join("scripts/configs");
    let common_path = config_dir.join("common.json");
    let compat_dir = config_dir.join("compat");

    let common_raw = std::fs::read_to_string(&common_path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", common_path.display()));
    #[derive(Debug, Deserialize)]
    struct CommonConfig {
        paths: HashMap<String, String>,
        comparison: ComparisonConfig,
    }

    #[derive(Debug, Deserialize)]
    struct CasesConfig {
        cases: Vec<CaseConfig>,
    }

    let common: CommonConfig = serde_json::from_str(&common_raw)
        .unwrap_or_else(|err| panic!("failed to parse {}: {err}", common_path.display()));
    let mut compat_files: Vec<PathBuf> = std::fs::read_dir(&compat_dir)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", compat_dir.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().and_then(|ext| ext.to_str()) == Some("json"))
        .collect();
    compat_files.sort();

    let mut all_cases = Vec::new();
    for path in compat_files {
        let raw = std::fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
        let cases: CasesConfig = serde_json::from_str(&raw)
            .unwrap_or_else(|err| panic!("failed to parse {}: {err}", path.display()));
        all_cases.extend(cases.cases);
    }

    CaseManifest {
        paths: common.paths,
        comparison: common.comparison,
        cases: all_cases,
    }
}

fn resolve_path(manifest: &CaseManifest, raw: &str) -> PathBuf {
    let mut resolved = raw.to_owned();
    for (key, value) in &manifest.paths {
        resolved = resolved.replace(&format!("{{{key}}}"), value);
    }
    project_root().join(resolved)
}

fn resolve_arg(manifest: &CaseManifest, raw: &str) -> String {
    let mut resolved = raw.to_owned();
    for (key, value) in &manifest.paths {
        let replacement = project_root().join(value).display().to_string();
        resolved = resolved.replace(&format!("{{{key}}}"), &replacement);
    }
    resolved
}

fn build_compute_args(manifest: &CaseManifest, case: &CaseConfig, output: &Path) -> Vec<String> {
    let mut args = Vec::new();
    args.push(case.mode.clone());

    args.push("-R".to_owned());
    args.extend(
        case.region_files
            .iter()
            .map(|path| resolve_path(manifest, path).display().to_string()),
    );

    args.push("-S".to_owned());
    args.extend(
        case.score_files
            .iter()
            .map(|path| resolve_path(manifest, path).display().to_string()),
    );

    args.extend(case.options.iter().map(|arg| resolve_arg(manifest, arg)));
    args.push("-o".to_owned());
    args.push(output.display().to_string());
    args
}

fn run_case(manifest: &CaseManifest, case: &CaseConfig) {
    let reference = case
        .reference
        .as_ref()
        .unwrap_or_else(|| panic!("case {} has no reference matrix", case.id));
    let reference_path = resolve_path(manifest, &reference.matrix);
    let output = NamedTempFile::new().unwrap();
    let output_path = output.path().to_path_buf();

    let args = build_compute_args(manifest, case, &output_path);
    let status = Command::new(compute_matrix_bin())
        .args(&args)
        .status()
        .unwrap_or_else(|err| panic!("failed to run compute_matrix_rs for {}: {err}", case.id));
    assert!(status.success(), "case {} failed with {status}", case.id);

    let candidate = load_matrix(&output_path)
        .unwrap_or_else(|err| panic!("failed to read candidate for {}: {err:#}", case.id));
    let expected = load_matrix(&reference_path).unwrap_or_else(|err| {
        panic!(
            "failed to read reference {} for {}: {err:#}",
            reference_path.display(),
            case.id
        )
    });

    let diff = full_diff(
        &candidate.header_json,
        &expected.header_json,
        &candidate.rows,
        &expected.rows,
        manifest.comparison.default_tolerance,
        &manifest.comparison.ignore_header_keys,
        manifest.comparison.max_diffs,
    );

    assert!(
        diff.matches,
        "case {} mismatched against {}\nheader_diffs={:?}\nbed_diff={:?}\nvalue_result={:?}",
        case.id,
        reference_path.display(),
        diff.header_diffs,
        diff.bed_diff,
        diff.value_result
    );
}

#[test]
fn manifest_compat_cases_match_references() {
    let manifest = load_manifest();
    let mut count = 0usize;

    for case in &manifest.cases {
        if case.tags.iter().any(|tag| tag == "compat") {
            count += 1;
            run_case(&manifest, case);
        }
    }

    assert!(count > 0, "manifest contains no compat cases");
}

#[test]
fn outfilename_matrix_header_consistent_after_skipzeros() {
    let manifest = load_manifest();
    let case = manifest
        .cases
        .iter()
        .find(|case| case.id == "corner_case_skipzeros_removes_zero_region")
        .expect("missing skipZeros case");

    let mat_tmp = tempfile::NamedTempFile::new().unwrap();
    let tab_tmp = tempfile::NamedTempFile::new().unwrap();

    let mut args = build_compute_args(&manifest, case, mat_tmp.path());
    args.push("--outFileNameMatrix".to_owned());
    args.push(tab_tmp.path().display().to_string());

    let status = Command::new(compute_matrix_bin())
        .args(&args)
        .status()
        .unwrap();
    assert!(status.success());

    let tab = std::fs::read_to_string(tab_tmp.path()).unwrap();
    let lines: Vec<&str> = tab.lines().collect();
    assert!(
        lines.len() >= 4,
        "expected >=3 header lines + 1 data row, got {}",
        lines.len()
    );

    let line1_stripped = lines[0].strip_prefix('#').unwrap_or_else(|| {
        panic!("line 1 should start with '#': {}", lines[0]);
    });
    let line1_counts: Vec<usize> = line1_stripped
        .split('\t')
        .map(|p| p.split(':').nth(1).unwrap().parse::<usize>().unwrap())
        .collect();
    let line1_total: usize = line1_counts.iter().sum();

    let line3_parts: Vec<&str> = lines[2].split('\t').collect();
    let num_groups = line1_counts.len();
    let line3_counts: Vec<usize> = line3_parts[..num_groups]
        .iter()
        .map(|p| p.split(':').nth(1).unwrap().parse::<usize>().unwrap())
        .collect();

    let data_rows = lines.len() - 3;

    assert_eq!(
        line1_counts, line3_counts,
        "line 1 and line 3 group counts diverge in:\n{tab}"
    );
    assert_eq!(
        line1_total, data_rows,
        "header total {line1_total} != data rows {data_rows} in:\n{tab}"
    );
}

#[test]
fn scale_zero_skipzeros_fails_when_all_rows_filter_out() {
    let manifest = load_manifest();
    let output = NamedTempFile::new().unwrap();
    let signal = resolve_path(&manifest, "{heatmapper}/test.bw");
    let regions = resolve_path(&manifest, "{heatmapper}/group1.bed");

    let result = Command::new(compute_matrix_bin())
        .args([
            "reference-point",
            "-R",
            regions.to_str().unwrap(),
            "-S",
            signal.to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--binSize",
            "10",
            "-p",
            "1",
            "--scale",
            "0",
            "--skipZeros",
            "-o",
            output.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        !result.status.success(),
        "--scale 0 --skipZeros should fail after filtering every row"
    );

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("No regions remain after runtime row filtering"),
        "unexpected stderr:\n{stderr}"
    );
}

#[test]
fn runtime_empty_group_preserved_with_zero_count() {
    use std::io::Write;

    let manifest = load_manifest();
    let mut genes_bed = tempfile::NamedTempFile::new().unwrap();
    writeln!(genes_bed, "ch1\t100\t150\thas_signal\t0\t+").unwrap();

    let mut peaks_bed = tempfile::NamedTempFile::new().unwrap();
    writeln!(peaks_bed, "ch1\t300\t350\tno_signal\t0\t+").unwrap();

    let mat_tmp = tempfile::NamedTempFile::new().unwrap();
    let tab_tmp = tempfile::NamedTempFile::new().unwrap();
    let signal = resolve_path(&manifest, "{heatmapper}/test.bw");

    let status = Command::new(compute_matrix_bin())
        .args([
            "reference-point",
            "-R",
            genes_bed.path().to_str().unwrap(),
            peaks_bed.path().to_str().unwrap(),
            "--smartLabels",
            "-S",
            signal.to_str().unwrap(),
            "-b",
            "100",
            "-a",
            "100",
            "--binSize",
            "10",
            "-p",
            "1",
            "--missingDataAsZero",
            "--skipZeros",
            "-o",
            mat_tmp.path().to_str().unwrap(),
            "--outFileNameMatrix",
            tab_tmp.path().to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let tab = std::fs::read_to_string(tab_tmp.path()).unwrap();
    let lines: Vec<&str> = tab.lines().collect();
    assert!(
        lines.len() >= 4,
        "expected >=3 header lines + 1 data row, got {}",
        lines.len()
    );

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

    let line3_parts: Vec<&str> = lines[2].split('\t').collect();
    let line3_genes: usize = line3_parts[0].split(':').nth(1).unwrap().parse().unwrap();
    let line3_peaks: usize = line3_parts[1].split(':').nth(1).unwrap().parse().unwrap();
    assert_eq!(line3_genes, 1, "line 3 genes count should match line 1");
    assert_eq!(line3_peaks, 0, "line 3 peaks count should match line 1");

    let data_rows = lines.len() - 3;
    assert_eq!(
        data_rows, 1,
        "should have exactly 1 data row, got {data_rows}"
    );

    let mat = load_matrix(mat_tmp.path()).unwrap();
    let boundaries: Vec<usize> = mat.header_json["group_boundaries"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    assert_eq!(boundaries, vec![0, 1, 1]);

    let labels = mat.header_json["group_labels"].as_array().unwrap();
    assert_eq!(labels.len(), 2, "should have 2 group labels");
    assert!(!labels[0].as_str().unwrap().is_empty());
    assert!(!labels[1].as_str().unwrap().is_empty());
}
