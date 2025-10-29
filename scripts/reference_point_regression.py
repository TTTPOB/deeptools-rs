#!/usr/bin/env python3
"""Generate regression data for computeMatrix reference-point mode.

This script uses the deepTools reference implementation (via pixi) to
produce a gzipped matrix output for a fixed BED/bigWig pair. It then runs the
Rust reimplementation with the same arguments and compares the resulting
matrices field-by-field, reporting the maximum absolute deviation observed.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import subprocess
import sys
import time
import urllib.request
import resource
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


@dataclass
class MatrixRow:
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
    header: dict
    rows: List[MatrixRow]

    @property
    def bin_count(self) -> int:
        boundaries = self.header.get("sample_boundaries")
        if not boundaries:
            return 0
        return int(boundaries[-1])

    @property
    def sample_count(self) -> int:
        boundaries = self.header.get("sample_boundaries")
        if not boundaries:
            return 0
        return max(len(boundaries) - 1, 0)


class CommandError(RuntimeError):
    pass


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
        "--bin-size",
        type=int,
        default=10,
        help="Bin size passed to computeMatrix (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write artefacts (default: target/reference-point-regression)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-5,
        help="Acceptable absolute tolerance when comparing matrix values (default: %(default)s)",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Keep existing outputs instead of regenerating them",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print executed commands and comparison details",
    )
    return parser.parse_args()


def ensure_commands_available(args: argparse.Namespace) -> None:
    for exe in (args.pixi, args.cargo):
        if shutil.which(exe) is None:
            raise CommandError(f"Required executable '{exe}' not found on PATH")


@dataclass
class CommandTiming:
    wall_seconds: float
    user_seconds: float
    system_seconds: float


def run_command(
    command: Sequence[str], *, cwd: Path | None = None, verbose: bool = False
) -> CommandTiming:
    if verbose:
        print("+", " ".join(str(part) for part in command))
    start = time.perf_counter()
    usage_before_self = resource.getrusage(resource.RUSAGE_SELF)
    usage_before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        subprocess.run(command, cwd=cwd, check=True)
    except subprocess.CalledProcessError as exc:
        raise CommandError(
            f"Command failed with exit code {exc.returncode}: {' '.join(command)}"
        ) from exc
    end = time.perf_counter()
    usage_after_self = resource.getrusage(resource.RUSAGE_SELF)
    usage_after_children = resource.getrusage(resource.RUSAGE_CHILDREN)

    user_seconds = (usage_after_children.ru_utime - usage_before_children.ru_utime) + (
        usage_after_self.ru_utime - usage_before_self.ru_utime
    )
    system_seconds = (
        usage_after_children.ru_stime - usage_before_children.ru_stime
    ) + (usage_after_self.ru_stime - usage_before_self.ru_stime)

    return CommandTiming(
        wall_seconds=end - start,
        user_seconds=user_seconds,
        system_seconds=system_seconds,
    )


def format_timing(label: str, timing: CommandTiming) -> str:
    return (
        f"{label}: "
        f"wall={timing.wall_seconds:.2f}s, "
        f"user={timing.user_seconds:.2f}s, "
        f"sys={timing.system_seconds:.2f}s"
    )


def load_matrix(path: Path) -> Matrix:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
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


def parse_value(field: str) -> float:
    lowered = field.lower()
    if lowered == "nan":
        return float("nan")
    return float(field)


def compare_headers(reference: dict, candidate: dict) -> List[str]:
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


def almost_equal(a: float, b: float, tol: float) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    return abs(a - b) <= tol


def compare_rows(
    reference: MatrixRow, candidate: MatrixRow, tol: float
) -> Tuple[float, List[str]]:
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

# update: use larger dataset and more complex conditions
# bw files:
# encode:
# k562_1: https://www.encodeproject.org/files/ENCFF093IIW/@@download/ENCFF093IIW.bigWig
# k562_2: https://www.encodeproject.org/files/ENCFF019IPA/@@download/ENCFF019IPA.bigWig
# k562_3: https://www.encodeproject.org/files/ENCFF656ZKM/@@download/ENCFF656ZKM.bigWig
# bed file:
# k562: https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz
# use peak center for reference point, use whole peak with +- 2 reagion for scale-region
# split the peak file into 2 files in half, to simulate two regions.

REFERENCE_POINT_SIGNALS = [
    (
        "ENCFF093IIW.bigWig",
        "https://www.encodeproject.org/files/ENCFF093IIW/@@download/ENCFF093IIW.bigWig",
    ),
    (
        "ENCFF019IPA.bigWig",
        "https://www.encodeproject.org/files/ENCFF019IPA/@@download/ENCFF019IPA.bigWig",
    ),
    (
        "ENCFF656ZKM.bigWig",
        "https://www.encodeproject.org/files/ENCFF656ZKM/@@download/ENCFF656ZKM.bigWig",
    ),
]

REFERENCE_POINT_BED = (
    "ENCFF333TAT.bed",
    "https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz",
)


def prepare_reference_point_dataset(
    dataset_root: Path, *, verbose: bool = False
) -> Tuple[List[Path], List[Path]]:
    dataset_root.mkdir(parents=True, exist_ok=True)

    signals: List[Path] = []
    for filename, url in REFERENCE_POINT_SIGNALS:
        destination = dataset_root / filename
        ensure_downloaded(url, destination, verbose=verbose)
        signals.append(destination)

    bed_name, bed_url = REFERENCE_POINT_BED
    bed_gz = dataset_root / f"{bed_name}.gz"
    bed_plain = dataset_root / bed_name

    ensure_downloaded(bed_url, bed_gz, verbose=verbose)
    if not bed_plain.exists() or bed_plain.stat().st_mtime < bed_gz.stat().st_mtime:
        if verbose:
            print(f"Decompressing {bed_gz.name} -> {bed_plain.name}")
        with gzip.open(bed_gz, "rb") as src, open(bed_plain, "wb") as dst:
            shutil.copyfileobj(src, dst)

    region_dir = dataset_root / "regions"
    region_dir.mkdir(exist_ok=True)
    region_paths = split_bed_file(bed_plain, region_dir, verbose=verbose)

    return region_paths, signals


def ensure_downloaded(url: str, destination: Path, *, verbose: bool = False) -> None:
    if destination.exists():
        return
    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    if verbose:
        print(f"Downloading {url} -> {destination}")
    try:
        with urllib.request.urlopen(url) as response, open(tmp_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
    except Exception as exc:  # pragma: no cover - network failure path
        if tmp_path.exists():
            tmp_path.unlink()
        raise CommandError(f"Failed to download {url}: {exc}") from exc
    tmp_path.rename(destination)


def split_bed_file(
    source: Path, target_dir: Path, *, verbose: bool = False
) -> List[Path]:
    with open(source, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    comments = [line for line in lines if line.startswith("#")]
    records = [line for line in lines if line.strip() and not line.startswith("#")]

    if not records:
        raise CommandError(f"BED file '{source}' does not contain any records")

    midpoint = (len(records) + 1) // 2
    parts = [
        ("part1", records[:midpoint]),
        ("part2", records[midpoint:]),
    ]

    region_paths: List[Path] = []
    for suffix, subset in parts:
        if not subset:
            continue
        region_path = target_dir / f"{source.stem}_{suffix}.bed"
        if verbose:
            print(f"Writing {region_path} with {len(subset)} records")
        with open(region_path, "w", encoding="utf-8") as handle:
            handle.writelines(comments)
            handle.writelines(subset)
        region_paths.append(region_path)

    return region_paths


def main() -> int:
    args = parse_args()
    ensure_commands_available(args)

    repo_root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir or (
        repo_root / "target" / "reference-point-regression"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_root = output_dir / "datasets" / "encode_k562_atac"
    regions, signals = prepare_reference_point_dataset(
        dataset_root, verbose=args.verbose
    )

    if not regions:
        raise FileNotFoundError("No region files available for regression harness")
    if not signals:
        raise FileNotFoundError(
            "No bigWig signal files available for regression harness"
        )

    suffix = args.reference_point.lower()
    reference_output = output_dir / f"reference_{suffix}.mat.gz"
    candidate_output = output_dir / f"rust_{suffix}.mat.gz"

    if not args.keep:
        if reference_output.exists():
            reference_output.unlink()
        if candidate_output.exists():
            candidate_output.unlink()

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

    region_args: List[str] = ["-R"]
    region_args.extend(str(region) for region in regions)

    signal_args: List[str] = ["-S"]
    signal_args.extend(str(signal) for signal in signals)

    reference_cmd = (
        [
            args.pixi,
            "run",
            "computeMatrix",
            "reference-point",
        ]
        + region_args
        + signal_args
        + command_common
        + ["--outFileName", str(reference_output)]
    )

    candidate_cmd = (
        [
            args.cargo,
            "run",
            "--release",
            "--quiet",
            "--",
            "reference-point",
        ]
        + region_args
        + signal_args
        + command_common
        + ["--outFileName", str(candidate_output)]
    )

    if args.verbose:
        print("Using regions:")
        for region in regions:
            print(f"  - {region}")
        print("Using signals:")
        for signal in signals:
            print(f"  - {signal}")
        print(f"Writing outputs to: {output_dir}")

    if not args.keep or not reference_output.exists():
        reference_timing = run_command(reference_cmd, cwd=repo_root, verbose=args.verbose)
    elif args.verbose:
        print(f"Skipping reference generation, file already exists: {reference_output}")
        reference_timing = None

    if not args.keep or not candidate_output.exists():
        candidate_timing = run_command(candidate_cmd, cwd=repo_root, verbose=args.verbose)
    elif args.verbose:
        print(f"Skipping Rust run, file already exists: {candidate_output}")
        candidate_timing = None

    reference_matrix = load_matrix(reference_output)
    candidate_matrix = load_matrix(candidate_output)

    ok, max_delta, issues = compare_matrices(
        reference_matrix, candidate_matrix, args.tolerance, args.verbose
    )

    if ok:
        print("✅ Matrices match within tolerance")
        print(
            f"   Samples: {reference_matrix.sample_count}, bins/sample: {reference_matrix.bin_count // max(reference_matrix.sample_count, 1)}"
        )
        print(f"   Rows compared: {len(reference_matrix.rows)}")
        print(f"   Max abs delta: {max_delta:.3e}")
        if reference_timing:
            print("   " + format_timing("Python (pixi)", reference_timing))
        if candidate_timing:
            print("   " + format_timing("Rust (cargo)", candidate_timing))
        return 0

    print("❌ Matrices differ")
    for issue in issues:
        print(f" - {issue}")
    print(f"Max abs delta observed: {max_delta:.3e}")
    return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CommandError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
