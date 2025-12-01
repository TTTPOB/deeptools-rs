# Python Compatibility Verification System for computeMatrix

## Overview

This document outlines a systematic approach to ensure the Rust reimplementation of `computeMatrix` produces **identical results** to the Python DeepTools implementation. The system extracts and automates all 13 Python test cases from `deeptools/test/test_heatmapper.py`, compares against reference matrix files, and integrates with the existing regression harness.

**Goal**: Achieve 100% compatibility with Python DeepTools computeMatrix across all documented test scenarios, with numerical matrix accuracy and downstream tool interoperability.

**Note**: This plan aligns with the project's validation strategy (overall_status.md) to maintain a unified regression harness that drives both `reference-point` and `scale-regions` modes while sharing datasets under `target/compute-matrix-datasets/`.

---

## DeepTools Python Test Coverage Analysis

Based on analysis of `deeptools/deeptools/test/test_heatmapper.py`, the Python test suite includes **13 functional tests** covering computeMatrix functionality:

### Core Functionality Tests (8 tests)

| Test Name | Description | Reference Matrix |
|-----------|-------------|------------------|
| `test_computeMatrix_reference_point` | Basic reference-point mode (TSS default) | `master.mat` |
| `test_computeMatrix_reference_point_center` | Reference-point with `--referencePoint center` | `master_center.mat` |
| `test_computeMatrix_reference_point_tes` | Reference-point with `--referencePoint TES` | `master_TES.mat` |
| `test_computeMatrix_reference_point_missing_data_as_zero` | Reference-point with `--missingDataAsZero` | `master_nan_to_zero.mat` |
| `test_computeMatrix_scale_regions` | Basic scale-regions mode | `master_scale_reg.mat` |
| `test_computeMatrix_multiple_bed` | Multiple BED files as input (group handling) | `master_multibed.mat` |
| `test_computeMatrix_region_extend_over_chr_end` | Regions extending beyond chromosome boundaries | `master_extend_beyond_chr_size.mat` |
| `test_computeMatrix_unscaled` | Scale-regions with `--unscaled5prime` and `--unscaled3prime` | `master_unscaled.mat` |

### Advanced Functionality Tests (2 tests)

| Test Name | Description | Reference Matrix |
|-----------|-------------|------------------|
| `test_computeMatrix_gtf` | GTF file input with scale-regions | `master_gtf.mat` |
| `test_computeMatrix_metagene` | Metagene mode with `--metagene` flag | `master_metagene.mat` |

### Low-level Unit Tests (4 tests - Zone Chopping)

| Test Name | Description |
|-----------|-------------|
| `test_chopRegions_body` | Region chopping for body mode (scale-regions) |
| `test_chopRegions_TSS` | Region chopping for TSS reference point |
| `test_chopRegions_TES` | Region chopping for TES reference point |
| `test_chopRegionsFromMiddle` | Region chopping for center reference point |

### Test Data Assets

Located in `deeptools/deeptools/test/test_heatmapper/`:
- **Signal files**: `test.bw`, `unscaled.bigWig`
- **Region files**: `test2.bed`, `group1.bed`, `group2.bed`, `unscaled.bed`
- **Reference matrices**: 10 `.mat` files serving as ground truth
- **GTF data**: `../test_data/test.gtf`, `../test_data/test1.bw.bw`

---

## Existing Infrastructure Analysis

### Current Regression Script (`scripts/compute_matrix_regression.py`)

The existing 970-line script provides robust infrastructure that should be **preserved and extended**:

#### Key Components to Retain

1. **`CommandTiming` dataclass** (lines 155-175): Tracks wall/user/system time
2. **`CachedResult` dataclass** (lines 200-235): Stores command hash, output path, timing, timestamp
3. **`run_command_cached()`** (lines 320-380): Hash-based execution caching with JSON persistence
4. **`Matrix` / `MatrixRow` classes** (lines 35-55): Matrix parsing and representation
5. **`compare_matrices()`** (lines 580-610): Tolerance-based value comparison
6. **`compare_headers()`** (lines 500-520): JSON header field comparison
7. **`ensure_downloaded()`** / `prepare_dataset()`**: Dataset download and caching

#### Cache Directory Structure (Existing)
```
target/<mode>-regression/
├── .cache/
│   ├── reference/           # Python command timing cache
│   │   └── <hash>.json
│   └── rust/                # Rust command timing cache
│       └── <hash>.json
├── <mode>_reference_<hash>.mat.gz
└── <mode>_rust_<hash>.mat.gz
```

---

## Proposed Three-Tier Validation System

### Tier 1: Direct Reference Matrix Comparison

Compare Rust output directly against the 10 pre-computed Python reference matrices stored in `test_heatmapper/`.

```python
def validate_against_reference(test_name: str, rust_output: Path, expected_matrix: Path) -> ValidationResult:
    """
    Compare Rust-generated matrix against Python reference matrix.
    This is the primary validation tier - if this passes, the implementation is correct.
    """
    rust_matrix = load_matrix(rust_output)
    ref_matrix = load_matrix(expected_matrix)
    
    return ValidationResult(
        test_name=test_name,
        header_match=compare_headers(ref_matrix.header, rust_matrix.header),
        value_match=compare_matrix_values(ref_matrix, rust_matrix, tolerance=0),
        row_count_match=len(ref_matrix.rows) == len(rust_matrix.rows)
    )
```

### Tier 2: Live Python-Rust Cross-Validation

Run both implementations with identical arguments and compare outputs.

```python
def cross_validate_implementations(test_scenario: TestScenario) -> CrossValidationResult:
    """
    Execute both Python and Rust implementations and compare outputs.
    Useful for scenarios not covered by reference matrices or when testing new features.
    """
    python_output = run_python_command(test_scenario.python_command)
    rust_output = run_rust_command(test_scenario.rust_command)
    
    return CrossValidationResult(
        test_name=test_scenario.name,
        python_output=python_output,
        rust_output=rust_output,
        matrices_identical=compare_matrices(python_output, rust_output, tolerance=1e-10),
        python_timing=get_timing(python_output),
        rust_timing=get_timing(rust_output)
    )
```

### Tier 3: Downstream Tool Compatibility

Verify that Rust-generated matrices work correctly with downstream DeepTools utilities.

```python
def verify_downstream_compatibility(rust_matrix: Path) -> DownstreamResult:
    """
    Ensure generated matrices are compatible with plotHeatmap, plotProfile, etc.
    """
    results = {}
    
    # Test plotHeatmap
    results['plotHeatmap'] = run_command([
        'pixi', 'run', 'plotHeatmap', '-m', str(rust_matrix),
        '-out', '/tmp/test_heatmap.png'
    ])
    
    # Test plotProfile  
    results['plotProfile'] = run_command([
        'pixi', 'run', 'plotProfile', '-m', str(rust_matrix),
        '-out', '/tmp/test_profile.png'
    ])
    
    return DownstreamResult(compatibility=all(r.success for r in results.values()))
```

---

## Implementation Architecture

### Module Structure

```
scripts/
├── compute_matrix_regression.py      # Main entry point (enhanced)
├── regression/
│   ├── __init__.py
│   ├── core/
│   │   ├── __init__.py
│   │   ├── timing.py                 # CommandTiming, CachedResult (extracted)
│   │   ├── cache.py                  # run_command_cached, load_cache, save_cache
│   │   └── matrix.py                 # Matrix, MatrixRow, load_matrix
│   ├── comparison/
│   │   ├── __init__.py
│   │   ├── header_compare.py         # compare_headers()
│   │   ├── value_compare.py          # compare_matrix_values(), almost_equal()
│   │   └── reporter.py               # Comparison result formatting
│   ├── test_extraction/
│   │   ├── __init__.py
│   │   ├── python_test_parser.py     # Parse test_heatmapper.py
│   │   └── scenario_generator.py     # Generate test scenarios from parsed tests
│   ├── datasets/
│   │   ├── __init__.py
│   │   ├── downloader.py             # ensure_downloaded(), prepare_dataset()
│   │   └── test_data_manager.py      # Manage test_heatmapper/ assets
│   └── runners/
│       ├── __init__.py
│       ├── python_runner.py          # Execute pixi-based Python commands
│       └── rust_runner.py            # Execute cargo-based Rust commands
└── config/
    └── python_compatibility.yaml     # Test scenario definitions
```

### Test Scenario Configuration

```yaml
# config/python_compatibility.yaml
test_suites:
  python_compatibility:
    description: "All 10 computeMatrix tests from test_heatmapper.py"
    data_root: "deeptools/deeptools/test/test_heatmapper"
    test_data_root: "deeptools/deeptools/test/test_data"
    
    scenarios:
      # Core Functionality Tests
      - name: reference_point_basic
        python_test: test_computeMatrix_reference_point
        reference_matrix: master.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/test2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: reference_point_center
        python_test: test_computeMatrix_reference_point_center
        reference_matrix: master_center.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/test2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 --referencePoint center -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: reference_point_tes
        python_test: test_computeMatrix_reference_point_tes
        reference_matrix: master_TES.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/test2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 --referencePoint TES -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: reference_point_missing_data_as_zero
        python_test: test_computeMatrix_reference_point_missing_data_as_zero
        reference_matrix: master_nan_to_zero.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/test2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 -bs 1 -p 1 --missingDataAsZero
          --outFileName {output}
        tolerance: 0
        
      - name: scale_regions_basic
        python_test: test_computeMatrix_scale_regions
        reference_matrix: master_scale_reg.mat
        rust_command: >
          cargo run --release --quiet -- scale-regions
          -R {data_root}/test2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 -m 100 -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: multiple_bed
        python_test: test_computeMatrix_multiple_bed
        reference_matrix: master_multibed.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/group1.bed {data_root}/group2.bed
          -S {data_root}/test.bw
          -b 100 -a 100 -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: region_extend_beyond_chr
        python_test: test_computeMatrix_region_extend_over_chr_end
        reference_matrix: master_extend_beyond_chr_size.mat
        rust_command: >
          cargo run --release --quiet -- reference-point
          -R {data_root}/group1.bed {data_root}/group2.bed
          -S {data_root}/test.bw
          -b 100 -a 500 -bs 1 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: scale_regions_unscaled
        python_test: test_computeMatrix_unscaled
        reference_matrix: master_unscaled.mat
        rust_command: >
          cargo run --release --quiet -- scale-regions
          -R {data_root}/unscaled.bed
          -S {data_root}/unscaled.bigWig
          -a 300 -b 500 --unscaled5prime 100 --unscaled3prime 50 -bs 10 -p 1
          --outFileName {output}
        tolerance: 0

      # Advanced Functionality Tests
      - name: gtf_input
        python_test: test_computeMatrix_gtf
        reference_matrix: master_gtf.mat
        rust_command: >
          cargo run --release --quiet -- scale-regions
          -R {test_data_root}/test.gtf
          -S {test_data_root}/test1.bw.bw
          -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 -bs 10 -p 1
          --outFileName {output}
        tolerance: 0
        
      - name: metagene
        python_test: test_computeMatrix_metagene
        reference_matrix: master_metagene.mat
        rust_command: >
          cargo run --release --quiet -- scale-regions
          -R {test_data_root}/test.gtf
          -S {test_data_root}/test1.bw.bw
          -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 -bs 10 -p 1 --metagene
          --outFileName {output}
        tolerance: 0

  # Extended testing with ENCODE data (existing infrastructure)
  encode_k562_atac:
    description: "Large-scale performance testing with ENCODE K562 ATAC-seq data"
    data_root: "target/compute-matrix-datasets/encode_k562_atac"
    scenarios:
      - name: reference_point_center_encode
        cross_validate: true  # Run both Python and Rust, compare outputs
        tolerance: 1e-5
        # ... (existing ENCODE-based tests)
```

---

## Caching Strategy

### Preserved Cache Mechanisms

The existing caching system will be **fully retained** with these enhancements:

```python
# Extended cache structure
target/python-compatibility/
├── .cache/
│   ├── rust/                          # Rust command timing cache
│   │   └── <scenario_hash>.json       # {command_hash, output_path, timing, timestamp}
│   └── reference/                     # Optional: regenerated Python reference cache
│       └── <scenario_hash>.json
├── outputs/
│   └── <scenario_name>_rust.mat.gz    # Rust-generated matrices
└── reports/
    └── compatibility_<date>.json      # Test results summary
```

### Cache Key Generation

```python
def compute_scenario_hash(scenario: TestScenario) -> str:
    """
    Generate a stable hash for caching based on:
    - Rust command (excluding output filename)
    - Reference matrix path
    - Test data file hashes (optional, for integrity)
    """
    hasher = hashlib.sha256()
    
    # Hash command parameters (excluding output)
    filtered_cmd = filter_output_arg(scenario.rust_command)
    hasher.update(filtered_cmd.encode())
    
    # Hash reference matrix identifier
    hasher.update(scenario.reference_matrix.encode())
    
    return hasher.hexdigest()[:16]
```

### Cache Invalidation Rules

1. **Reference Matrix Change**: Invalidate if reference `.mat` file is modified
2. **Rust Binary Change**: Invalidate on `cargo build` if binary timestamp changes
3. **Command Parameter Change**: Automatic via hash-based key generation
4. **Manual Override**: `--no-cache` flag to force re-execution

---

## CLI Interface

### New Command-Line Options

```bash
# Run all Python compatibility tests against reference matrices
python scripts/compute_matrix_regression.py --mode python-compatibility

# Run specific test suite
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --suite python_compatibility

# Run single test scenario
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --test reference_point_center

# Generate detailed JSON report
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --output reports/compatibility_$(date +%Y%m%d).json

# Cross-validate with live Python execution (Tier 2)
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --cross-validate

# Verify downstream tool compatibility (Tier 3)
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --verify-downstream

# Force re-execution (bypass cache)
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --no-cache

# Keep existing Rust outputs (reuse cached)
python scripts/compute_matrix_regression.py --mode python-compatibility \
    --keep-rust
```

### Backward Compatibility

Existing CLI interface remains fully functional:

```bash
# Original reference-point regression (unchanged)
python scripts/compute_matrix_regression.py --mode reference-point

# Original scale-regions regression (unchanged)
python scripts/compute_matrix_regression.py --mode scale-regions
```

---

## Reporting System

### JSON Report Structure

```json
{
  "compatibility_report": {
    "timestamp": "2025-12-01T10:30:00Z",
    "total_tests": 10,
    "passed": 10,
    "failed": 0,
    "skipped": 0,
    "byte_for_byte_matches": 10,
    
    "test_results": {
      "reference_point_basic": {
        "status": "PASS",
        "reference_matrix": "master.mat",
        "rust_output": "reference_point_basic_rust.mat.gz",
        "header_match": true,
        "value_match": true,
        "row_count": 6,
        "max_delta": 0.0,
        "rust_timing": {
          "wall_seconds": 0.15,
          "user_seconds": 0.12,
          "system_seconds": 0.02
        },
        "from_cache": false
      },
      "reference_point_center": {
        "status": "PASS",
        "reference_matrix": "master_center.mat",
        "header_match": true,
        "value_match": true,
        "row_count": 6,
        "max_delta": 0.0,
        "rust_timing": {
          "wall_seconds": 0.14,
          "user_seconds": 0.11,
          "system_seconds": 0.02
        },
        "from_cache": true
      }
      // ... remaining tests
    },
    
    "performance_summary": {
      "average_rust_time": 0.145,
      "total_rust_time": 1.45,
      "cached_tests": 3,
      "fresh_tests": 7
    }
  }
}
```

### Console Output Format

```
================================================================================
📊 PYTHON COMPATIBILITY TEST RESULTS
================================================================================

Test Suite: python_compatibility (10 tests)

  ✅ reference_point_basic ............ PASS (0.15s)
  ✅ reference_point_center ........... PASS (cached)
  ✅ reference_point_tes .............. PASS (0.14s)
  ✅ reference_point_missing_data ..... PASS (0.16s)
  ✅ scale_regions_basic .............. PASS (0.18s)
  ✅ multiple_bed ..................... PASS (0.13s)
  ✅ region_extend_beyond_chr ......... PASS (0.15s)
  ✅ scale_regions_unscaled ........... PASS (0.17s)
  ✅ gtf_input ........................ PASS (0.22s)
  ✅ metagene ......................... PASS (0.21s)

================================================================================
SUMMARY: 10/10 tests passed (100% compatibility)
         10/10 byte-for-byte matches
         Total time: 1.51s (3 cached, 7 fresh)
================================================================================
```

---

## Implementation Plan

### Phase 1: Test Extraction & Configuration (Week 1)

**Tasks:**
1. Extract module structure from existing `compute_matrix_regression.py`
   - Move `CommandTiming`, `CachedResult` to `regression/core/timing.py`
   - Move `Matrix`, `MatrixRow`, `load_matrix` to `regression/core/matrix.py`
   - Move cache functions to `regression/core/cache.py`
   - Move comparison functions to `regression/comparison/`

2. Create test scenario configuration
   - Write `config/python_compatibility.yaml` with all 10 test scenarios
   - Implement YAML loader in `regression/test_extraction/scenario_generator.py`

3. Implement reference matrix comparison
   - Add Tier 1 validation in `regression/comparison/value_compare.py`
   - Handle uncompressed `.mat` files from `test_heatmapper/`

**Success Criteria:**
- All 10 test scenarios defined in YAML
- Reference matrix loading works for all `.mat` files
- Cache mechanism preserved and functional

### Phase 2: Core Validation Implementation (Week 1-2)

**Tasks:**
1. Implement test runner
   - Add `--mode python-compatibility` CLI option
   - Create scenario execution loop with caching
   - Integrate with existing `run_command_cached()`

2. Implement Tier 1 validation
   - Compare Rust output against reference matrices
   - Report header mismatches, value deltas, row count differences

3. Add reporting system
   - JSON report generation
   - Console output formatting
   - Detailed error messages for failures

**Success Criteria:**
- All 10 tests can be executed via CLI
- Caching works correctly for repeated runs
- JSON and console reports generated

### Phase 3: Cross-Validation & Downstream Testing (Week 2)

**Tasks:**
1. Implement Tier 2 cross-validation (optional)
   - Add `--cross-validate` flag
   - Execute Python via pixi, compare with Rust output
   - Track performance differences

2. Implement Tier 3 downstream validation (optional)
   - Add `--verify-downstream` flag
   - Test `plotHeatmap`, `plotProfile` with Rust matrices

3. Integration with existing ENCODE tests
   - Merge python-compatibility mode with existing regression tests
   - Unified reporting for all test types

**Success Criteria:**
- Cross-validation functional for all scenarios
- Downstream tools work with Rust-generated matrices
- Single unified test harness for all validation types

### Phase 4: CI/CD Integration (Week 3)

**Tasks:**
1. Create GitHub Actions workflow
   - Run python-compatibility tests on PR/push
   - Cache test data and reference matrices
   - Upload test reports as artifacts

2. Add pre-commit hook (optional)
   - Run subset of fast tests before commit
   - Block commits with failing compatibility tests

3. Documentation
   - Update `AGENTS.md` with new testing workflow
   - Document all CLI options and configuration format

**Success Criteria:**
- Automated testing on every PR
- Clear pass/fail indication in CI
- Comprehensive documentation

---

## Success Metrics

### Primary Compatibility Criteria

| Metric | Target | Measurement |
|--------|--------|-------------|
| Test Pass Rate | 100% | All 10 Python test scenarios pass |
| Numerical Accuracy | 100% | Zero tolerance deviations (tolerance=1e-5) |
| Header Compatibility | 100% | All JSON header fields match |
| Row Count Accuracy | 100% | Exact row counts in all scenarios |

### Secondary Performance Criteria

| Metric | Target | Measurement |
|--------|--------|-------------|
| Speedup vs Python | ≥2x | Average wall time ratio |
| CI Runtime | <5 min | Total GitHub Actions runtime |

### Downstream Compatibility

| Tool | Target | Validation |
|------|--------|------------|
| plotHeatmap | Works | Generates valid PNG output |
| plotProfile | Works | Generates valid PNG output |
| computeMatrixOperations | Works | Can read and manipulate matrices |

---

## Conclusion

This Python Compatibility Verification System provides a **systematic, automated approach** to ensure the Rust reimplementation produces identical results to Python DeepTools. By:

1. **Extracting all 13 Python test cases** into a structured configuration
2. **Leveraging the 10 reference matrices** as ground truth for numerical validation
3. **Preserving the existing cache mechanism** for efficient repeated testing
4. **Implementing a three-tier validation system** (reference, cross-validation, downstream)
5. **Integrating with CI/CD** for automated verification on every code change

The system guarantees compatibility while maintaining the performance advantages of the Rust implementation. The modular architecture allows for easy extension as new test scenarios are identified or as the Rust implementation adds new features.