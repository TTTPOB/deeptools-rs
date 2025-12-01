"""
Python command runner.

Provides utilities for executing the Python DeepTools computeMatrix implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from ..core.cache import run_command_cached
from ..core.timing import CommandTiming


def build_python_command(
    pixi: str,
    mode: str,
    regions: List[Path],
    signals: List[Path],
    output_path: Path,
    extra_args: List[str],
) -> List[str]:
    """
    Build the complete Python computeMatrix command.

    Args:
        pixi: Path to pixi executable
        mode: computeMatrix mode (reference-point or scale-regions)
        regions: List of region file paths
        signals: List of signal file paths
        output_path: Path for the output file
        extra_args: Additional command arguments

    Returns:
        Complete command list
    """
    cmd = [pixi, "run", "computeMatrix", mode]

    cmd.append("-R")
    cmd.extend(str(r) for r in regions)

    cmd.append("-S")
    cmd.extend(str(s) for s in signals)

    cmd.extend(extra_args)
    cmd.extend(["--outFileName", str(output_path)])

    return cmd


def run_python_command(
    command: List[str],
    output_path: Path,
    cache_dir: Path,
    *,
    cwd: Path | None = None,
    verbose: bool = False,
    skip_cache: bool = False,
    cache_key: str | None = None,
    quiet: bool = False,
) -> tuple[CommandTiming, bool]:
    """
    Execute the Python computeMatrix command with caching.

    Args:
        command: Complete command list
        output_path: Path where output will be written
        cache_dir: Directory for cache files
        cwd: Working directory
        verbose: Enable verbose output
        skip_cache: Skip cache lookup
        cache_key: Optional cache key override
        quiet: Suppress progress output

    Returns:
        Tuple of (timing, from_cache)
    """
    return run_command_cached(
        command,
        output_path,
        cache_dir,
        cwd=cwd,
        verbose=verbose,
        skip_cache=skip_cache,
        cache_key=cache_key,
        quiet=quiet,
    )
