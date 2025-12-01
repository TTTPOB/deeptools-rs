"""
Caching utilities for command execution.

Extracted from compute_matrix_regression.py for modularity.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Sequence

from .command import run_command
from .timing import CachedResult, CommandTiming, format_timing


def compute_params_hash(command: Sequence[str]) -> str:
    """
    Compute a hash of command parameters excluding the output filename.

    This filters out --outFileName and its argument to create a stable hash
    based only on the computation parameters.
    """
    filtered_cmd = []
    skip_next = False

    for i, part in enumerate(command):
        if skip_next:
            skip_next = False
            continue

        if part == "--outFileName":
            skip_next = True  # Skip the next argument (the filename)
            continue

        filtered_cmd.append(str(part))

    cmd_str = " ".join(filtered_cmd)
    return hashlib.sha256(cmd_str.encode()).hexdigest()[
        :16
    ]  # Use shorter hash for filenames


def load_cache(cache_file: Path) -> CachedResult | None:
    """Load a single cache entry from disk."""
    if not cache_file.exists():
        return None
    try:
        with open(cache_file, "r", encoding="utf-8") as f:
            data = json.load(f)
            return CachedResult.from_dict(data)
    except (json.JSONDecodeError, KeyError, ValueError) as exc:
        print(
            f"Warning: Failed to load cache from {cache_file}: {exc}", file=sys.stderr
        )
        return None


def save_cache(cached_result: CachedResult, cache_file: Path) -> None:
    """Save a single cache entry to disk."""
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cached_result.to_dict(), f, indent=2)


def run_command_cached(
    command: Sequence[str],
    output_path: Path,
    cache_dir: Path,
    *,
    cwd: Path | None = None,
    verbose: bool = False,
    skip_cache: bool = False,
    enable_cache: bool = True,
    cache_key: str | None = None,
    quiet: bool = False,
) -> tuple[CommandTiming, bool]:
    """
    Run a command with caching support.

    Args:
        command: Command and arguments to execute
        output_path: Path where the command writes its output
        cache_dir: Directory to store cache files
        cwd: Working directory for the command
        verbose: Print detailed progress information
        skip_cache: If True, skip cache lookup but still save result
        enable_cache: If False, disable all caching
        cache_key: Optional override for the cache filename and stored command hash
        quiet: If True, suppress progress output

    Returns:
        tuple of (CommandTiming, from_cache: bool)
    """
    cmd_hash = cache_key or compute_params_hash(command)
    cache_file = cache_dir / f"{cmd_hash}.json"
    use_cache = enable_cache and not skip_cache

    if use_cache:
        cached = load_cache(cache_file)

        if cached is not None:
            # Check if output file still exists and is newer than cache entry
            if Path(cached.output_path).exists():
                if not quiet:
                    print("🔄 Using cached result")
                if verbose:
                    print(f"   Cache file: {cache_file}")
                    print(f"   Cached output: {cached.output_path}")
                if not quiet:
                    print("   " + format_timing("Cached timing", cached.timing))
                return cached.timing, True
            else:
                if not quiet:
                    print(f"⚠️  Cache hit but output file missing: {cached.output_path}")

    # Run the command
    if not quiet:
        if enable_cache:
            print("🚀 Running command")
            if verbose:
                print(f"   Cache file: {cache_file}")
        else:
            print("🚀 Running command (cache disabled)")

    timing = run_command(command, cwd=cwd, verbose=verbose)

    if not quiet:
        print("   " + format_timing("Execution timing", timing))

    # Save to cache only if caching is enabled
    if use_cache:
        cached_result = CachedResult(
            command_hash=cmd_hash,
            output_path=str(output_path.resolve()),
            timing=timing,
            timestamp=time.time(),
            command=[str(part) for part in command],
        )
        save_cache(cached_result, cache_file)
        if verbose:
            print(f"💾 Cached result to {cache_file}")

    return timing, False


def load_timing_if_kept(
    output: Path,
    cache_key: str,
    cache_dir: Path,
    keep: bool,
    no_cache: bool,
    verbose: bool,
) -> tuple[CommandTiming | None, bool]:
    """
    Load cached timing if output exists and keeping is enabled.

    Returns:
        tuple of (CommandTiming or None, from_cache: bool)
    """
    if not keep or not output.exists():
        return None, False

    if verbose:
        print(f"Output exists (keeping): {output}")

    if no_cache:
        return None, False

    cache_file = cache_dir / f"{cache_key}.json"
    cached = load_cache(cache_file)
    if cached is not None:
        if verbose:
            print(f"  Loaded cached timing info from {cache_file}")
        return cached.timing, True

    return None, False
