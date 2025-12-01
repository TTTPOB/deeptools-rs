#!/usr/bin/env python3
"""Generate regression data for computeMatrix.

This harness uses the deepTools reference implementation (via pixi) to
produce gzipped matrix outputs for a fixed BED/bigWig corpus. It then runs the
Rust reimplementation with the same arguments and compares the resulting
matrices field-by-field, reporting the maximum absolute deviation observed.

Both `reference-point` and `scale-regions` subcommands are supported.

Python Compatibility Mode:
  The `--mode python-compatibility` option runs all 10 test cases from
  test_heatmapper.py against pre-computed reference matrices, validating
  that the Rust implementation produces identical output.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import List

from regression.comparison import (
    CompatibilityReport,
    TestResult,
    compare_matrices,
)
from regression.comparison.reporter import (
    create_report,
    format_summary,
    format_test_result,
    generate_json_report,
)

# Import from refactored modules
from regression.core import (
    CommandError,
    CommandTiming,
    compute_params_hash,
    load_cache,
    load_matrix,
    run_command_cached,
)
from regression.core.cache import load_timing_if_kept
from regression.core.timing import format_timing
from regression.datasets import (
    prepare_dataset,
)
from regression.datasets.test_data_manager import get_test_data_paths
from regression.test_extraction import (
    TestScenario,
    TestSuite,
    load_reference_matrix,
    load_test_config,
)
from regression.test_extraction.scenario_generator import get_default_config_path


DEFAULT_TOLERANCE = 5e-6


def _require_multiple(value: int, divisor: int, flag: str) -> None:
    if value % divisor != 0:
        raise CommandError(
            f"{flag} ({value}) must be a multiple of the bin size ({divisor})"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pixi",
        default="pixi",
        help="Path to the pixi executable (default: %(default)s)",
    )
    parser.add_argument(
        "--cargo",
        default="cargo",
        help="Path to the cargo executable (default: %(default)s)",
    )
    parser.add_argument(
        "--mode",
        choices=["reference-point", "scale-regions", "python-compatibility"],
        default="reference-point",
        help="computeMatrix subcommand or test mode (default: %(default)s)",
    )
    parser.add_argument(
        "--suite",
        default="python_compatibility",
        help="Test suite to run in python-compatibility mode (default: %(default)s)",
    )
    parser.add_argument(
        "--test",
        default=None,
        help="Run only a specific test scenario by name",
    )
    parser.add_argument(
        "--reference-point",
        choices=["TSS", "TES", "center"],
        default="center",
        help="Reference point passed to computeMatrix (default: %(default)s)",
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=100,
        help="Bases upstream of the reference point (default: %(default)s)",
    )
    parser.add_argument(
        "--downstream",
        type=int,
        default=100,
        help="Bases downstream of the reference point (default: %(default)s)",
    )
    parser.add_argument(
        "--region-body-length",
        type=int,
        default=200,
        help="Body length used for scale-regions mode (default: %(default)s)",
    )
    parser.add_argument(
        "--unscaled-5-prime",
        type=int,
        default=50,
        help="Unscaled bases at the 5' end for scale-regions mode (default: %(default)s)",
    )
    parser.add_argument(
        "--unscaled-3-prime",
        type=int,
        default=50,
        help="Unscaled bases at the 3' end for scale-regions mode (default: %(default)s)",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=10,
        help="Bin size passed to computeMatrix (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write artefacts (default: target/<mode>-regression)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path for JSON report output (python-compatibility mode)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=(
            "Absolute tolerance for matrix value comparison (default/max: %(default)s; "
            "smaller values are allowed, larger values are clamped)"
        ),
    )
    parser.add_argument(
        "--keep-ref",
        dest="keep_ref",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Reuse existing reference outputs and cached timing (use --no-keep-ref to regenerate)",
    )
    parser.add_argument(
        "--keep-rust",
        dest="keep_rust",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Reuse existing Rust outputs and cached timing (use --no-keep-rust to regenerate)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print executed commands and comparison details",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Skip cache and always run commands fresh",
    )
    parser.add_argument(
        "--cross-validate",
        action="store_true",
        help="Also run Python implementation and compare (python-compatibility mode)",
    )
    parser.add_argument(
        "--verify-downstream",
        action="store_true",
        help="Verify compatibility with downstream tools like plotHeatmap (python-compatibility mode)",
    )
    return parser.parse_args()


def ensure_commands_available(args: argparse.Namespace) -> None:
    for exe in (args.pixi, args.cargo):
        if shutil.which(exe) is None:
            raise CommandError(f"Required executable '{exe}' not found on PATH")


def resolve_tolerance(*candidates: float) -> float:
    """Clamp comparison tolerance to the strictest positive value (<= DEFAULT_TOLERANCE).

    Accepts heterogeneous inputs (float, int, str) and ignores invalid/zero/None.
    """

    def _as_float(value):
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    positives = []
    for candidate in candidates:
        numeric = _as_float(candidate)
        if numeric is not None and numeric > 0:
            positives.append(numeric)

    if not positives:
        return DEFAULT_TOLERANCE

    return min(min(positives), DEFAULT_TOLERANCE)


def run_python_compatibility_mode(args: argparse.Namespace) -> int:
    """Run tests against Python DeepTools reference matrices."""
    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or (repo_root / "target" / "python-compatibility")
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = output_dir / ".cache" / "rust"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Load test configuration
    config_path = get_default_config_path()
    if not config_path.exists():
        print(f"❌ Configuration file not found: {config_path}", file=sys.stderr)
        return 2

    suites = load_test_config(config_path)
    if args.suite not in suites:
        print(f"❌ Test suite '{args.suite}' not found", file=sys.stderr)
        print(f"   Available suites: {list(suites.keys())}", file=sys.stderr)
        return 2

    suite = suites[args.suite]
    test_paths = get_test_data_paths(repo_root)

    # Filter scenarios if a specific test is requested
    scenarios = suite.scenarios
    if args.test:
        scenarios = [s for s in scenarios if s.name == args.test]
        if not scenarios:
            print(f"❌ Test '{args.test}' not found in suite '{args.suite}'", file=sys.stderr)
            available = [s.name for s in suite.scenarios]
            print(f"   Available tests: {available}", file=sys.stderr)
            return 2

    print("=" * 70)
    print("📊 PYTHON COMPATIBILITY TEST SUITE")
    print("=" * 70)
    print(f"Suite: {args.suite}")
    print(f"Tests: {len(scenarios)}")
    print(f"Data root: {test_paths.data_root}")
    print("=" * 70)
    print()

    results: List[TestResult] = []

    for scenario in scenarios:
        if not scenario.enabled:
            result = TestResult(
                name=scenario.name,
                status="SKIP",
                reference_matrix=scenario.reference_matrix,
            )
            results.append(result)
            print(format_test_result(result))
            continue

        try:
            result = run_single_scenario(
                scenario=scenario,
                repo_root=repo_root,
                data_root=test_paths.data_root,
                test_data_root=test_paths.test_data_root,
                output_dir=output_dir,
                cache_dir=cache_dir,
                args=args,
            )
        except Exception as exc:
            result = TestResult(
                name=scenario.name,
                status="ERROR",
                reference_matrix=scenario.reference_matrix,
                error_message=str(exc),
            )
            if args.verbose:
                import traceback
                traceback.print_exc()

        results.append(result)
        print(format_test_result(result))

    # Generate report
    report = create_report(results)
    print(format_summary(report))

    # Write JSON report if requested
    if args.output:
        generate_json_report(report, args.output)
        print(f"\n📄 Report written to: {args.output}")

    # Return appropriate exit code
    if report.failed > 0 or report.errors > 0:
        return 1
    return 0


def run_single_scenario(
    scenario: TestScenario,
    repo_root: Path,
    data_root: Path,
    test_data_root: Path,
    output_dir: Path,
    cache_dir: Path,
    args: argparse.Namespace,
) -> TestResult:
    """Run a single test scenario and return the result."""

    # Locate reference matrix
    reference_path = load_reference_matrix(scenario.reference_matrix, data_root)

    # Build Rust command
    rust_output = output_dir / "outputs" / f"{scenario.name}_rust.mat.gz"
    rust_output.parent.mkdir(parents=True, exist_ok=True)

    rust_cmd = scenario.get_rust_command(
        cargo=args.cargo,
        data_root=data_root,
        test_data_root=test_data_root,
        output_path=rust_output,
    )

    # Check if we should reuse existing output
    cache_key = compute_params_hash(rust_cmd)

    if args.keep_rust and rust_output.exists():
        timing, from_cache = load_timing_if_kept(
            rust_output, cache_key, cache_dir, args.keep_rust, args.no_cache, args.verbose
        )
        if timing is None:
            # Output exists but no cached timing - load from cache or run
            timing, from_cache = run_command_cached(
                rust_cmd,
                rust_output,
                cache_dir,
                cwd=repo_root,
                verbose=args.verbose,
                skip_cache=args.no_cache,
                cache_key=cache_key,
                quiet=True,
            )
    else:
        # Remove existing output if not keeping
        if rust_output.exists():
            rust_output.unlink()

        timing, from_cache = run_command_cached(
            rust_cmd,
            rust_output,
            cache_dir,
            cwd=repo_root,
            verbose=args.verbose,
            skip_cache=args.no_cache,
            cache_key=cache_key,
            quiet=True,
        )

    # Load and compare matrices
    ref_matrix = load_matrix(reference_path)
    rust_matrix = load_matrix(rust_output)

    tolerance = resolve_tolerance(scenario.tolerance, args.tolerance)
    ok, max_delta, issues = compare_matrices(ref_matrix, rust_matrix, tolerance, args.verbose)

    # Determine header/value match status
    header_issues = [i for i in issues if "header" in i.lower() or "key" in i.lower()]
    value_issues = [i for i in issues if "value" in i.lower() or "delta" in i.lower() or "row" in i.lower()]

    return TestResult(
        name=scenario.name,
        status="PASS" if ok else "FAIL",
        reference_matrix=scenario.reference_matrix,
        rust_output=str(rust_output),
        header_match=len(header_issues) == 0,
        value_match=len(value_issues) == 0,
        row_count=len(rust_matrix.rows),
        max_delta=max_delta,
        rust_timing=timing,
        from_cache=from_cache,
        issues=issues if not ok else [],
    )


def run_encode_regression_mode(args: argparse.Namespace) -> int:
    """Run the original ENCODE-based regression test."""
    repo_root = Path(__file__).resolve().parents[1]
    mode = args.mode
    output_dir = args.output_dir or (repo_root / "target" / f"{mode}-regression")
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_root = output_dir / ".cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    reference_cache_dir = cache_root / "reference"
    reference_cache_dir.mkdir(parents=True, exist_ok=True)
    rust_cache_dir = cache_root / "rust"
    rust_cache_dir.mkdir(parents=True, exist_ok=True)

    dataset_root = repo_root / "target" / "compute-matrix-datasets" / "encode_k562_atac"
    regions, signals = prepare_dataset(dataset_root, verbose=args.verbose)

    if not regions:
        raise FileNotFoundError("No region files available for regression harness")
    if not signals:
        raise FileNotFoundError(
            "No bigWig signal files available for regression harness"
        )

    _require_multiple(args.upstream, args.bin_size, "--beforeRegionStartLength")
    _require_multiple(args.downstream, args.bin_size, "--afterRegionStartLength")
    if mode == "scale-regions":
        _require_multiple(args.region_body_length, args.bin_size, "--regionBodyLength")
        _require_multiple(args.unscaled_5_prime, args.bin_size, "--unscaled5prime")
        _require_multiple(args.unscaled_3_prime, args.bin_size, "--unscaled3prime")

    if mode == "reference-point":
        command_common = [
            "--beforeRegionStartLength",
            str(args.upstream),
            "--afterRegionStartLength",
            str(args.downstream),
            "--referencePoint",
            args.reference_point,
            "--binSize",
            str(args.bin_size),
            "--numberOfProcessors",
            "4",
        ]
    else:
        command_common = [
            "--regionBodyLength",
            str(args.region_body_length),
            "--beforeRegionStartLength",
            str(args.upstream),
            "--afterRegionStartLength",
            str(args.downstream),
            "--unscaled5prime",
            str(args.unscaled_5_prime),
            "--unscaled3prime",
            str(args.unscaled_3_prime),
            "--binSize",
            str(args.bin_size),
            "--numberOfProcessors",
            "4",
        ]

    region_args: List[str] = ["-R"]
    region_args.extend(str(region) for region in regions)

    signal_args: List[str] = ["-S"]
    signal_args.extend(str(signal) for signal in signals)

    # Build command without output filename to compute hash
    base_reference_cmd = (
        [
            args.pixi,
            "run",
            "computeMatrix",
            mode,
        ]
        + region_args
        + signal_args
        + command_common
    )

    base_candidate_cmd = (
        [
            args.cargo,
            "run",
            "--release",
            "--quiet",
            "--",
            mode,
        ]
        + region_args
        + signal_args
        + command_common
    )

    # Compute hash from parameters (excluding output filename)
    reference_hash = compute_params_hash(base_reference_cmd)

    # Use reference hash in output filenames (shared across implementations)
    reference_output = output_dir / f"{mode}_reference_{reference_hash}.mat.gz"
    candidate_output = output_dir / f"{mode}_rust_{reference_hash}.mat.gz"

    if not args.keep_ref and reference_output.exists():
        reference_output.unlink()
    if not args.keep_rust and candidate_output.exists():
        candidate_output.unlink()

    # Build complete commands with output filenames
    reference_cmd = base_reference_cmd + ["--outFileName", str(reference_output)]
    candidate_cmd = base_candidate_cmd + ["--outFileName", str(candidate_output)]

    if args.verbose:
        print("Using regions:")
        for region in regions:
            print(f"  - {region}")
        print("Using signals:")
        for signal in signals:
            print(f"  - {signal}")
        print(f"Writing outputs to: {output_dir}")

    reference_timing = None
    reference_from_cache = False
    if args.keep_ref and reference_output.exists():
        reference_timing, reference_from_cache = load_timing_if_kept(
            reference_output,
            reference_hash,
            reference_cache_dir,
            args.keep_ref,
            args.no_cache,
            args.verbose,
        )

    if reference_timing is None:
        reference_timing, reference_from_cache = run_command_cached(
            reference_cmd,
            reference_output,
            reference_cache_dir,
            cwd=repo_root,
            verbose=args.verbose,
            skip_cache=args.no_cache,
            enable_cache=True,
            cache_key=reference_hash,
        )

    candidate_timing = None
    candidate_from_cache = False
    if args.keep_rust and candidate_output.exists():
        candidate_timing, candidate_from_cache = load_timing_if_kept(
            candidate_output,
            reference_hash,
            rust_cache_dir,
            args.keep_rust,
            args.no_cache,
            args.verbose,
        )

    if candidate_timing is None:
        candidate_timing, candidate_from_cache = run_command_cached(
            candidate_cmd,
            candidate_output,
            rust_cache_dir,
            cwd=repo_root,
            verbose=args.verbose,
            skip_cache=args.no_cache,
            enable_cache=True,
            cache_key=reference_hash,
        )

    reference_matrix = load_matrix(reference_output)
    candidate_matrix = load_matrix(candidate_output)

    tolerance = resolve_tolerance(args.tolerance)
    ok, max_delta, issues = compare_matrices(
        reference_matrix, candidate_matrix, tolerance, args.verbose
    )

    print("\n" + "=" * 60)
    print("📊 PERFORMANCE SUMMARY")
    print("=" * 60)
    if reference_timing:
        print(format_timing("Python (pixi)", reference_timing, reference_from_cache))
    if candidate_timing:
        print(format_timing("Rust (cargo)", candidate_timing, candidate_from_cache))

    if (
        reference_timing
        and candidate_timing
        and not reference_from_cache
        and not candidate_from_cache
    ):
        speedup = reference_timing.wall_seconds / candidate_timing.wall_seconds
        print(f"Speedup: {speedup:.2f}x (wall time)")
    print("=" * 60)

    if ok:
        print(f"✅ Matrices match within tolerance (≤ {tolerance:.1e})")
        print(
            f"   Samples: {reference_matrix.sample_count}, bins/sample: {reference_matrix.bin_count // max(reference_matrix.sample_count, 1)}"
        )
        print(f"   Rows compared: {len(reference_matrix.rows)}")
        print(f"   Max abs delta: {max_delta:.3e}")
        return 0

    print(f"❌ Matrices differ (tolerance {tolerance:.1e})")
    for issue in issues:
        print(f" - {issue}")
    print(f"Max abs delta observed: {max_delta:.3e}")
    return 1


def main() -> int:
    args = parse_args()
    ensure_commands_available(args)

    if args.mode == "python-compatibility":
        return run_python_compatibility_mode(args)
    else:
        return run_encode_regression_mode(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CommandError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
