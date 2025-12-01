"""
Dataset download and management utilities.

Extracted from compute_matrix_regression.py for modularity.
"""

from __future__ import annotations

import gzip
import shutil
import urllib.request
from pathlib import Path
from typing import List, Tuple

from ..core.command import CommandError

# ENCODE K562 ATAC-seq dataset URLs
HARNESS_SIGNALS = [
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

HARNESS_BED = (
    "ENCFF333TAT.bed",
    "https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz",
)


def ensure_downloaded(url: str, destination: Path, *, verbose: bool = False) -> None:
    """
    Download a file from URL if it doesn't exist.

    Args:
        url: Source URL to download from
        destination: Local path to save the file
        verbose: Print download progress

    Raises:
        CommandError: If download fails
    """
    if destination.exists():
        return

    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    if verbose:
        print(f"Downloading {url} -> {destination}")

    try:
        with urllib.request.urlopen(url) as response, open(tmp_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
    except Exception as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise CommandError(f"Failed to download {url}: {exc}") from exc

    tmp_path.rename(destination)


def split_bed_file(
    source: Path, target_dir: Path, *, verbose: bool = False
) -> List[Path]:
    """
    Split a BED file into two parts for multi-group testing.

    Args:
        source: Source BED file path
        target_dir: Directory to write split files
        verbose: Print progress

    Returns:
        List of paths to the split BED files
    """
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


def prepare_dataset(
    dataset_root: Path, *, verbose: bool = False
) -> Tuple[List[Path], List[Path]]:
    """
    Prepare the ENCODE K562 ATAC-seq dataset for regression testing.

    Downloads files if needed and splits the BED file into two groups.

    Args:
        dataset_root: Directory to store the dataset
        verbose: Print progress

    Returns:
        Tuple of (region_paths, signal_paths)
    """
    dataset_root.mkdir(parents=True, exist_ok=True)

    signals: List[Path] = []
    for filename, url in HARNESS_SIGNALS:
        destination = dataset_root / filename
        ensure_downloaded(url, destination, verbose=verbose)
        signals.append(destination)

    bed_name, bed_url = HARNESS_BED
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
