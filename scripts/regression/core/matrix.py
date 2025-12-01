"""
Matrix representation and loading utilities.

Extracted from compute_matrix_regression.py for modularity.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class MatrixRow:
    """Represents a single row in the computeMatrix output."""

    chrom: str
    start: int
    end: int
    name: str
    score: float | None
    strand: str
    extra_fields: List[str]
    values: List[float]


@dataclass
class Matrix:
    """Represents a complete computeMatrix output with header and rows."""

    header: dict
    rows: List[MatrixRow]

    @property
    def bin_count(self) -> int:
        """Total number of bins (columns) in the matrix."""
        boundaries = self.header.get("sample_boundaries")
        if not boundaries:
            return 0
        return int(boundaries[-1])

    @property
    def sample_count(self) -> int:
        """Number of samples in the matrix."""
        boundaries = self.header.get("sample_boundaries")
        if not boundaries:
            return 0
        return max(len(boundaries) - 1, 0)


def parse_value(field: str) -> float:
    """Parse a matrix value field, handling NaN."""
    lowered = field.lower()
    if lowered == "nan":
        return float("nan")
    return float(field)


def load_matrix(path: Path) -> Matrix:
    """
    Load a computeMatrix output file.

    Supports both gzipped (.gz) and uncompressed (.mat) files.

    Args:
        path: Path to the matrix file

    Returns:
        Matrix object with header and rows
    """
    # Determine if file is gzipped based on extension or content
    def open_gzipped(p: Path):
        return gzip.open(p, "rt", encoding="utf-8")

    def open_plain(p: Path):
        return open(p, "r", encoding="utf-8")

    if path.suffix == ".gz":
        opener = open_gzipped
    else:
        opener = open_plain

    with opener(path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"Matrix file '{path}' is empty")
        header_line = header_line.strip()
        if header_line.startswith("@"):
            header_line = header_line[1:]
        header = json.loads(header_line)

        rows: List[MatrixRow] = []
        total_values = (
            int(header["sample_boundaries"][-1])
            if header.get("sample_boundaries")
            else 0
        )

        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Malformed row in '{path}': {line!r}")

            chrom, start, end, name, raw_score, strand = fields[:6]
            start_i = int(start)
            end_i = int(end)
            name = name or "."
            strand = strand or "."

            score: float | None
            if raw_score == ".":
                score = None
            else:
                score = float(raw_score)

            remainder = fields[6:]
            extra_count = max(len(remainder) - total_values, 0)
            extra_fields = remainder[:extra_count]
            value_fields = remainder[extra_count:]

            values = [parse_value(field) for field in value_fields]
            rows.append(
                MatrixRow(
                    chrom=chrom,
                    start=start_i,
                    end=end_i,
                    name=name,
                    score=score,
                    strand=strand,
                    extra_fields=extra_fields,
                    values=values,
                )
            )

    return Matrix(header=header, rows=rows)
