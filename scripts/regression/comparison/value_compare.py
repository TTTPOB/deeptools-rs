"""
Matrix value comparison utilities.
"""

from __future__ import annotations

import math
from typing import List, Tuple

from ..core.matrix import Matrix, MatrixRow


def almost_equal(a: float, b: float, tol: float) -> bool:
    """
    Check if two floats are approximately equal within tolerance.

    Handles NaN values: two NaN values are considered equal.
    """
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isnan(a) or math.isnan(b):
        return False
    return abs(a - b) <= tol


def compare_rows(
    reference: MatrixRow, candidate: MatrixRow, tol: float
) -> Tuple[float, List[str]]:
    """
    Compare two matrix rows and return differences.

    Args:
        reference: The expected row
        candidate: The actual row to compare
        tol: Tolerance for floating-point comparison

    Returns:
        Tuple of (max_delta, list of difference messages)
    """
    messages = []
    max_delta = 0.0

    if (
        reference.chrom,
        reference.start,
        reference.end,
        reference.name,
        reference.strand,
        reference.extra_fields,
    ) != (
        candidate.chrom,
        candidate.start,
        candidate.end,
        candidate.name,
        candidate.strand,
        candidate.extra_fields,
    ):
        messages.append(
            "Row metadata differs: reference="
            f"{reference.chrom}:{reference.start}-{reference.end} {reference.name} {reference.strand}"
            f", candidate={candidate.chrom}:{candidate.start}-{candidate.end} {candidate.name} {candidate.strand}"
        )

    if reference.score != candidate.score:
        if reference.score is None and candidate.score is None:
            pass
        else:
            messages.append(
                f"Score mismatch: reference={reference.score}, candidate={candidate.score}"
            )

    if len(reference.values) != len(candidate.values):
        messages.append(
            f"Value count mismatch: reference={len(reference.values)}, candidate={len(candidate.values)}"
        )
        return max_delta, messages

    for idx, (ref_value, cand_value) in enumerate(
        zip(reference.values, candidate.values)
    ):
        if not almost_equal(ref_value, cand_value, tol):
            messages.append(
                f"Value delta at bin {idx}: reference={ref_value}, candidate={cand_value}, "
                f"abs diff={abs(ref_value - cand_value):.3e}"
            )
        if not math.isnan(ref_value) and not math.isnan(cand_value):
            max_delta = max(max_delta, abs(ref_value - cand_value))

    return max_delta, messages


def compare_matrices(
    reference: Matrix, candidate: Matrix, tol: float, verbose: bool = False
) -> Tuple[bool, float, List[str]]:
    """
    Compare two matrices and return comparison results.

    Args:
        reference: The expected matrix
        candidate: The actual matrix to compare
        tol: Tolerance for floating-point comparison
        verbose: Include detailed row-level differences

    Returns:
        Tuple of (success, max_abs_delta, list of issues)
    """
    from .header_compare import compare_headers

    issues: List[str] = []

    header_issues = compare_headers(reference.header, candidate.header)
    issues.extend(header_issues)

    if len(reference.rows) != len(candidate.rows):
        issues.append(
            f"Row count mismatch: reference has {len(reference.rows)}, candidate has {len(candidate.rows)}"
        )
        return False, 0.0, issues

    max_abs_delta = 0.0

    for idx, (ref_row, cand_row) in enumerate(zip(reference.rows, candidate.rows)):
        row_delta, row_issues = compare_rows(ref_row, cand_row, tol)
        max_abs_delta = max(max_abs_delta, row_delta)
        if row_issues and verbose:
            issues.append(f"Row {idx} differences:\n  " + "\n  ".join(row_issues))
        elif row_issues:
            issues.append(f"Row {idx} differs (rerun with --verbose for details)")

    success = not issues
    return success, max_abs_delta, issues
