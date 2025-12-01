"""
Rust command runner.

Provides utilities for executing the Rust computeMatrix implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

from ..core.cache import run_command_cached
from ..core.timing import CommandTiming


def build_rust_command(
    cargo: str,
    rust_args: List[str],
    output_path: Path,
) -> List[str]:
    """
    Build the complete Rust command.

    Args:
        cargo: Path to cargo executable
        rust_args: Arguments for the Rust implementation
        output_path: Path for the output file

    Returns:
        Complete command list
    """
    cmd = [cargo, "run", "--release", "--quiet", "--"]
    cmd.extend(rust_args)
    cmd.extend(["--outFileName", str(output_path)])
    return cmd


def run_rust_command(
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
    Execute the Rust computeMatrix command with caching.

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
