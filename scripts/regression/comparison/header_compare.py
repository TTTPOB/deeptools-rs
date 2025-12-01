"""
Header comparison utilities for matrix validation.
"""

from __future__ import annotations

from typing import List


def compare_headers(reference: dict, candidate: dict) -> List[str]:
    """
    Compare two matrix headers and return a list of differences.

    Args:
        reference: The expected header dictionary
        candidate: The actual header dictionary to compare

    Returns:
        List of difference messages (empty if headers match)
    """
    messages = []
    reference_keys = set(reference.keys())
    candidate_keys = set(candidate.keys())

    missing = reference_keys - candidate_keys
    extra = candidate_keys - reference_keys
    if missing:
        messages.append(f"Missing header keys in candidate: {sorted(missing)}")
    if extra:
        messages.append(f"Unexpected header keys in candidate: {sorted(extra)}")

    shared = reference_keys & candidate_keys
    for key in sorted(shared):
        ref_value = reference[key]
        cand_value = candidate[key]
        if ref_value != cand_value:
            messages.append(
                f"Header mismatch for '{key}': reference={ref_value!r}, candidate={cand_value!r}"
            )

    return messages
