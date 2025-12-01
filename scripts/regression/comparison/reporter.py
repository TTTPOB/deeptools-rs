"""
Reporting utilities for test results.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from ..core.timing import CommandTiming


@dataclass
class TestResult:
    """Result of a single test scenario."""

    name: str
    status: str  # "PASS", "FAIL", "SKIP", "ERROR"
    reference_matrix: str
    rust_output: Optional[str] = None
    header_match: bool = True
    value_match: bool = True
    row_count: int = 0
    max_delta: float = 0.0
    rust_timing: Optional[CommandTiming] = None
    from_cache: bool = False
    error_message: Optional[str] = None
    issues: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        result = {
            "status": self.status,
            "reference_matrix": self.reference_matrix,
            "header_match": self.header_match,
            "value_match": self.value_match,
            "row_count": self.row_count,
            "max_delta": self.max_delta,
            "from_cache": self.from_cache,
        }
        if self.rust_output:
            result["rust_output"] = self.rust_output
        if self.rust_timing:
            result["rust_timing"] = self.rust_timing.to_dict()
        if self.error_message:
            result["error_message"] = self.error_message
        if self.issues:
            result["issues"] = self.issues
        return result


@dataclass
class CompatibilityReport:
    """Complete test report for a test suite."""

    timestamp: str
    total_tests: int = 0
    passed: int = 0
    failed: int = 0
    skipped: int = 0
    errors: int = 0
    byte_for_byte_matches: int = 0
    test_results: dict = field(default_factory=dict)
    performance_summary: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "compatibility_report": {
                "timestamp": self.timestamp,
                "total_tests": self.total_tests,
                "passed": self.passed,
                "failed": self.failed,
                "skipped": self.skipped,
                "errors": self.errors,
                "byte_for_byte_matches": self.byte_for_byte_matches,
                "test_results": self.test_results,
                "performance_summary": self.performance_summary,
            }
        }


def format_test_result(result: TestResult) -> str:
    """Format a single test result for console output."""
    status_icons = {
        "PASS": "✅",
        "FAIL": "❌",
        "SKIP": "⏭️",
        "ERROR": "💥",
    }
    icon = status_icons.get(result.status, "❓")

    # Determine timing info
    timing_info = ""
    if result.from_cache:
        timing_info = "(cached)"
    elif result.rust_timing:
        timing_info = f"({result.rust_timing.wall_seconds:.2f}s)"

    # Pad test name for alignment
    name_padded = result.name.ljust(35, ".")

    return f"  {icon} {name_padded} {result.status} {timing_info}"


def format_summary(report: CompatibilityReport) -> str:
    """Format the complete test report summary for console output."""
    lines = [
        "",
        "=" * 70,
        "📊 PYTHON COMPATIBILITY TEST RESULTS",
        "=" * 70,
        "",
    ]

    # Test results
    for name, result_dict in report.test_results.items():
        result = TestResult(name=name, **result_dict) if isinstance(result_dict, dict) else result_dict
        if isinstance(result_dict, dict):
            # Reconstruct TestResult from dict
            result = TestResult(
                name=name,
                status=result_dict.get("status", "UNKNOWN"),
                reference_matrix=result_dict.get("reference_matrix", ""),
                header_match=result_dict.get("header_match", True),
                value_match=result_dict.get("value_match", True),
                row_count=result_dict.get("row_count", 0),
                max_delta=result_dict.get("max_delta", 0.0),
                from_cache=result_dict.get("from_cache", False),
            )
        lines.append(format_test_result(result))

    # Summary
    lines.extend([
        "",
        "=" * 70,
        f"SUMMARY: {report.passed}/{report.total_tests} tests passed "
        f"({100 * report.passed / max(report.total_tests, 1):.0f}% compatibility)",
        f"         {report.byte_for_byte_matches}/{report.total_tests} byte-for-byte matches",
    ])

    # Performance summary
    if report.performance_summary:
        cached = report.performance_summary.get("cached_tests", 0)
        fresh = report.performance_summary.get("fresh_tests", 0)
        total_time = report.performance_summary.get("total_rust_time", 0)
        lines.append(f"         Total time: {total_time:.2f}s ({cached} cached, {fresh} fresh)")

    lines.append("=" * 70)

    return "\n".join(lines)


def generate_json_report(report: CompatibilityReport, output_path: Path) -> None:
    """Write the report to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report.to_dict(), f, indent=2)


def create_report(results: List[TestResult]) -> CompatibilityReport:
    """Create a compatibility report from test results."""
    report = CompatibilityReport(
        timestamp=datetime.now().isoformat(),
        total_tests=len(results),
    )

    total_time = 0.0
    cached_count = 0
    fresh_count = 0

    for result in results:
        report.test_results[result.name] = result.to_dict()

        if result.status == "PASS":
            report.passed += 1
            if result.max_delta == 0.0 and result.header_match:
                report.byte_for_byte_matches += 1
        elif result.status == "FAIL":
            report.failed += 1
        elif result.status == "SKIP":
            report.skipped += 1
        else:
            report.errors += 1

        if result.rust_timing:
            if result.from_cache:
                cached_count += 1
            else:
                fresh_count += 1
                total_time += result.rust_timing.wall_seconds

    report.performance_summary = {
        "total_rust_time": total_time,
        "cached_tests": cached_count,
        "fresh_tests": fresh_count,
    }
    if fresh_count > 0:
        report.performance_summary["average_rust_time"] = total_time / fresh_count

    return report
