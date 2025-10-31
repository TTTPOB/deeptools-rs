# Comprehensive Regression Test System for computeMatrix

## Overview

This document outlines a plan to rework the current monolithic regression script (`scripts/compute_matrix_regression.py`) into a modular, configurable test system that can handle many different computeMatrix scenarios. The new system will be split into reusable components and support configurable command lists for comprehensive testing.

**Note**: This plan aligns with the project's validation strategy (line 112-114 in overall_status.md) to maintain a unified regression harness that drives both `reference-point` and `scale-regions` modes while sharing datasets under `target/compute-matrix-datasets/`.
## DeepTools Python Test Coverage Analysis

Based on analysis of the DeepTools Python implementation (`deeptools/deeptools/test/test_heatmapper.py`), the current test suite includes **13 functional tests** covering computeMatrix functionality:

### Core Functionality Tests (8 tests)
1. **`test_computeMatrix_reference_point()`** - Basic reference-point mode with default parameters
2. **`test_computeMatrix_reference_point_center()`** - Reference-point with center reference point
3. **`test_computeMatrix_reference_point_tes()`** - Reference-point with TES reference point
4. **`test_computeMatrix_reference_point_missing_data_as_zero()`** - Reference-point with `--missingDataAsZero` flag
5. **`test_computeMatrix_scale_regions()`** - Basic scale-regions mode
6. **`test_computeMatrix_multiple_bed()`** - Multiple BED files as input
7. **`test_computeMatrix_region_extend_over_chr_end()`** - Regions extending beyond chromosome end
8. **`test_computeMatrix_unscaled()`** - Scale-regions with unscaled 5' and 3' regions

### Advanced Functionality Tests (2 tests)
9. **`test_computeMatrix_gtf()`** - GTF file input with scale-regions
10. **`test_computeMatrix_metagene()`** - Metagene mode with GTF input

### Low-level Unit Tests (3 tests)
11. **`test_chopRegions_body()`** - Tests region chopping for body mode
12. **`test_chopRegions_TSS()`** - Tests region chopping for TSS reference point
13. **`test_chopRegions_TES()`** - Tests region chopping for TES reference point
14. **`test_chopRegionsFromMiddle()`** - Tests region chopping from middle (center reference point)

### Test Data Coverage
The test suite uses **10 different reference matrix files**:
- `master.mat` - Basic reference-point
- `master_center.mat` - Reference-point with center
- `master_TES.mat` - Reference-point with TES
- `master_nan_to_zero.mat` - Missing data as zero
- `master_scale_reg.mat` - Scale-regions
- `master_multibed.mat` - Multiple BED files
- `master_extend_beyond_chr_size.mat` - Regions beyond chromosome end
- `master_unscaled.mat` - Scale-regions with unscaled regions
- `master_gtf.mat` - GTF input
- `master_metagene.mat` - Metagene mode

### Test Parameters Covered
- **Modes**: `reference-point`, `scale-regions`
- **Reference Points**: `center`, `TSS`, `TES`
- **Parameters**: upstream/downstream distances, bin size, region body length, unscaled regions
- **Data Types**: BED files, bigWig files, GTF files
- **Special Cases**: Missing data handling, chromosome boundary conditions, multiple input files
- **Strand Handling**: Both positive and negative strands in region chopping tests

### Gaps in Current Test Coverage
Based on the analysis, the following areas have limited or no test coverage:
1. **Parameter Matrix Testing**: No systematic testing of parameter combinations
2. **Performance Testing**: No large-scale performance benchmarks
3. **Error Handling**: Limited testing of error conditions and edge cases
4. **Output Formats**: Limited testing of different output formats (matrix, tabulated, BED)
5. **Sorting and Filtering**: No tests for `--sortRegions`, `--skipZeros`, thresholding
6. **Multiple Signal Files**: Limited testing with multiple bigWig inputs
7. **Different Average Types**: No testing of mean, median, min, max, sum, std
8. **Chromosome Name Variations**: No testing of chromosome name normalization


## Current Limitations

The existing regression script has several limitations:

1. **Monolithic Design**: All functionality is in a single 969-line file
2. **Limited Test Scenarios**: Only supports basic `reference-point` and `scale-regions` modes
3. **Hardcoded Configuration**: Test parameters are embedded in the code
4. **Fixed Dataset**: Only uses one specific dataset (ENCODE K562 ATAC-seq)
5. **Limited Extensibility**: Adding new test scenarios requires code modifications

## Proposed Modular Architecture

The new system will be split into four main modules:

### 1. data_downloader Module

**Purpose**: Handle downloading, caching, and preparation of test datasets

**Responsibilities**:
- Download test data from various sources (ENCODE, custom URLs, local files)
- Manage dataset caching and versioning under `target/compute-matrix-datasets/`
- Support multiple dataset types (bigWig, BED, GTF) with focus on current ENCODE K562 ATAC-seq data
- Prepare region files for testing (no file splitting)
- Validate dataset integrity
- Integrate with existing pixi environment for DeepTools reference data

**Key Classes**:
- `DatasetManager`: Orchestrates dataset operations
- `DataDownloader`: Handles URL downloads with retry logic
- `DatasetCache`: Manages local caching of datasets

### 2. task_runner Module

**Purpose**: Execute computeMatrix commands and manage results

**Responsibilities**:
- Execute provided command strings directly (no complex CLI argument mapping)
- Execute both reference (Python/pixi) and candidate (Rust/cargo) implementations
- Track execution timing and resource usage (extending existing `CommandTiming` class)
- Manage result caching (enhancing existing caching mechanism)
- Support both reference and candidate command strings
- Handle parallel execution of multiple test scenarios
- Support the existing `--keep-ref` and `--keep-rust` flags

**Key Classes**:
- `TaskRunner`: Main execution orchestrator
- `CommandBuilder`: Constructs command lines from parameters
- `ExecutionCache`: Manages result caching
- `PerformanceTracker`: Tracks timing and resource usage

### 3. matrix_comparer Module

**Purpose**: Compare matrix outputs and report differences

**Responsibilities**:
- Load and parse matrix files (extending existing `Matrix`, `MatrixRow` classes)
- Perform detailed comparisons (headers, values, metadata) with byte-for-byte accuracy
- Generate comprehensive difference reports
- Support various comparison tolerances and modes (maintaining existing tolerance logic)
- Export comparison results in multiple formats
- Ensure compatibility with downstream tools (`plotHeatmap`, `plotProfile`, `computeMatrixOperations`)
- Compare outputs from any command string pairs

**Key Classes**:
- `MatrixComparer`: Main comparison engine
- `MatrixLoader`: Loads and parses matrix files
- `DifferenceReporter`: Generates comparison reports
- `ComparisonResult`: Stores and formats comparison results

### 4. test_config Module

**Purpose**: Define and manage test scenarios and configurations

**Responsibilities**:
- Define test scenarios in configuration files
- Support parameter variations and combinations
- Manage test suites and test groups
- Provide configuration validation
- Support test scenario inheritance
- Support arbitrary command strings for both reference and candidate implementations
- Validate basic command structure

**Key Classes**:
- `TestConfig`: Loads and validates configuration files
- `TestScenario`: Represents a single test scenario
- `TestSuite`: Manages collections of test scenarios
- `ParameterGenerator`: Generates parameter combinations

## Configuration System

The new system will use YAML configuration files to define test scenarios:

```yaml
# test_scenarios.yaml
test_suites:
  basic_functionality:
    description: "Basic functionality tests"
    datasets:
      - encode_k562_atac
    scenarios:
      - name: "reference_point_center"
        reference_command: "pixi run computeMatrix reference-point -R target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed -S target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig --beforeRegionStartLength 100 --afterRegionStartLength 100 --referencePoint center --binSize 10 --numberOfProcessors 4 --outFileName {output}"
        candidate_command: "cargo run --release --quiet -- reference-point -R target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed -S target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig --beforeRegionStartLength 100 --afterRegionStartLength 100 --referencePoint center --binSize 10 --numberOfProcessors 4 --outFileName {output}"
        tolerance: 1e-5
      
      - name: "scale_regions_default"
        reference_command: "pixi run computeMatrix scale-regions -R target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed -S target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig --regionBodyLength 200 --beforeRegionStartLength 100 --afterRegionStartLength 100 --unscaled5prime 50 --unscaled3prime 50 --binSize 10 --numberOfProcessors 4 --outFileName {output}"
        candidate_command: "cargo run --release --quiet -- scale-regions -R target/compute-matrix-datasets/encode_k562_atac/ENCFF333TAT.bed -S target/compute-matrix-datasets/encode_k562_atac/ENCFF093IIW.bigWig target/compute-matrix-datasets/encode_k562_atac/ENCFF019IPA.bigWig --regionBodyLength 200 --beforeRegionStartLength 100 --afterRegionStartLength 100 --unscaled5prime 50 --unscaled3prime 50 --binSize 10 --numberOfProcessors 4 --outFileName {output}"
        tolerance: 1e-5

  custom_commands:
    description: "Custom command tests"
    scenarios:
      - name: "custom_test_1"
        reference_command: "pixi run computeMatrix reference-point -R /path/to/custom.bed -S /path/to/custom.bw --beforeRegionStartLength 500 --afterRegionStartLength 500 --referencePoint TSS --binSize 25 --outFileName {output}"
        candidate_command: "cargo run --release --quiet -- reference-point -R /path/to/custom.bed -S /path/to/custom.bw --beforeRegionStartLength 500 --afterRegionStartLength 500 --referencePoint TSS --binSize 25 --outFileName {output}"
        tolerance: 1e-6

datasets:
  encode_k562_atac:
    description: "ENCODE K562 ATAC-seq data"
    signals:
      - url: "https://www.encodeproject.org/files/ENCFF093IIW/@@download/ENCFF093IIW.bigWig"
        name: "k562_1"
      - url: "https://www.encodeproject.org/files/ENCFF019IPA/@@download/ENCFF019IPA.bigWig"
        name: "k562_2"
      - url: "https://www.encodeproject.org/files/ENCFF656ZKM/@@download/ENCFF656ZKM.bigWig"
        name: "k562_3"
    regions:
      - url: "https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz"
        name: "k562_peaks"
```

## Comprehensive Test Scenarios

The new system will support testing many different scenarios:

### 1. Basic Functionality Tests
- Both `reference-point` and `scale-regions` modes
- Different reference points (TSS, TES, center)
- Various upstream/downstream lengths
- Different bin sizes

### 2. Edge Case Tests
- Very small regions (smaller than bin size)
- Very large regions (megabase-scale)
- Regions at chromosome boundaries
- Empty region files
- Single-region files

### 3. Parameter Combination Tests
- All valid parameter combinations
- Boundary value testing
- Invalid parameter rejection
- Parameter inheritance and overrides

### 4. Performance Tests
- Large datasets (many regions, many signals)
- Memory usage profiling
- Execution time benchmarks
- Scalability testing

### 5. Data Format Tests
- Different bigWig formats
- Various BED formats (BED3, BED6, BED12)
- GTF format support
- Compressed input files

### 6. Error Handling Tests
- Missing files
- Corrupted data files
- Network failures
- Insufficient permissions

## Implementation Plan

### Phase 1: Core Module Development
1. **data_downloader Module**
   - Implement `DatasetManager` class
   - Add support for multiple data sources
   - Implement caching and versioning
   - Add data validation

2. **task_runner Module**
   - Implement `TaskRunner` class
   - Add command building logic
   - Implement execution caching
   - Add performance tracking

3. **matrix_comparer Module**
   - Implement `MatrixComparer` class
   - Add detailed comparison logic
   - Implement various report formats
   - Add tolerance handling

### Phase 2: Configuration System
1. **test_config Module**
   - Implement YAML configuration parser
   - Add parameter matrix expansion
   - Implement test scenario inheritance
   - Add configuration validation

2. **Configuration Files**
   - Create basic test scenario configurations
   - Add dataset definitions
   - Implement command string templates
   - Add output path substitution

### Phase 3: Integration and Testing
1. **Main Script Integration**
   - Rewrite main script to use new modules
   - Add command-line interface for configuration
   - Implement parallel test execution
   - Add progress reporting
   - Support command string execution

2. **Test Suite Expansion**
   - Create comprehensive test scenarios
   - Add edge case tests
   - Implement performance tests
   - Add error handling tests
   - Support custom command strings

### Phase 4: Advanced Features
1. **Reporting and Visualization**
   - Add HTML report generation
   - Implement trend analysis
   - Add performance graphs
   - Implement regression detection

2. **CI/CD Integration**
   - Add GitHub Actions integration
   - Implement automated test execution
   - Add result archiving
   - Implement notification system

## Directory Structure

```
scripts/
├── compute_matrix_regression.py          # Main entry point
├── regression/
│   ├── __init__.py
│   ├── data_downloader/
│   │   ├── __init__.py
│   │   ├── dataset_manager.py
│   │   ├── data_downloader.py
│   │   ├── dataset_cache.py
│   │   └── region_splitter.py
│   ├── task_runner/
│   │   ├── __init__.py
│   │   ├── task_runner.py
│   │   ├── command_builder.py
│   │   ├── execution_cache.py
│   │   └── performance_tracker.py
│   ├── matrix_comparer/
│   │   ├── __init__.py
│   │   ├── matrix_comparer.py
│   │   ├── matrix_loader.py
│   │   ├── difference_reporter.py
│   │   └── comparison_result.py
│   ├── test_config/
│   │   ├── __init__.py
│   │   ├── test_config.py
│   │   ├── test_scenario.py
│   │   ├── test_suite.py
│   │   └── parameter_generator.py
│   └── utils/
│       ├── __init__.py
│       ├── command_timing.py
│       ├── file_utils.py
│       └── validation.py
├── config/
│   ├── test_scenarios.yaml              # Main test configuration
│   ├── datasets.yaml                    # Dataset definitions
│   └── test_environments.yaml           # Environment-specific settings
└── reports/                             # Generated test reports
    ├── templates/
    │   ├── html_report.html
    │   └── json_schema.json
    └── output/
```

## Usage Examples

### Basic Usage
```bash
# Run all test scenarios
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml

# Run specific test suite
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --suite basic_functionality

# Run with custom tolerance
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --tolerance 1e-6

# Generate detailed report
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --report-format html --output reports/latest.html
```

### Advanced Usage
```bash
# Run with parallel execution
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --parallel 4

# Run performance benchmarks
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --suite performance --benchmark

# Run with custom dataset
python scripts/compute_matrix_regression.py --config config/test_scenarios.yaml --dataset custom_dataset.yaml

# Generate parameter matrix
python scripts/compute_matrix_regression.py --generate-param-matrix --output config/generated_scenarios.yaml
```

## Benefits of the New System

1. **Modularity**: Each component has a single responsibility and can be tested independently
2. **Extensibility**: New test scenarios can be added without code changes
3. **Reusability**: Components can be reused for other testing purposes
4. **Maintainability**: Smaller, focused modules are easier to maintain
5. **Flexibility**: Configuration-driven approach supports many testing scenarios
6. **Scalability**: Parallel execution and efficient caching support large test suites
7. **Reporting**: Comprehensive reports help identify and diagnose issues
8. **CI/CD Integration**: Designed for automated testing in continuous integration

## Migration Strategy

1. **Backward Compatibility**: Maintain compatibility with existing command-line interface
2. **Incremental Migration**: Implement modules incrementally while keeping the original script functional
3. **Parallel Development**: New modules can be developed in parallel with existing functionality
4. **Testing**: Each module will have comprehensive unit tests
5. **Documentation**: Detailed documentation for each module and the overall system

## Success Metrics

1. **Code Quality**: Reduced complexity and improved maintainability
2. **Test Coverage**: Increased test scenario coverage
3. **Execution Time**: Improved parallel execution and caching
4. **Extensibility**: Easier to add new test scenarios
5. **Reporting**: Better visibility into test results and trends
6. **CI/CD Integration**: Automated testing in continuous integration

## Conclusion

This modular regression testing system will significantly improve the testing capabilities for the computeMatrix reimplementation. By separating concerns into focused modules and using a configuration-driven approach, we can create a comprehensive, extensible testing framework that ensures the reliability and correctness of the Rust implementation.