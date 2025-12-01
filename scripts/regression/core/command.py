"""
Command execution utilities.

Extracted from compute_matrix_regression.py for modularity.
"""

from __future__ import annotations

import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence

from .timing import CommandTiming, create_timing


class CommandError(RuntimeError):
    """Raised when a command fails to execute."""

    pass


def run_command(
    command: Sequence[str], *, cwd: Path | None = None, verbose: bool = False
) -> CommandTiming:
    """
    Execute a command and return timing information.

    Args:
        command: Command and arguments to execute
        cwd: Working directory for the command
        verbose: Print the command before execution

    Returns:
        CommandTiming with wall, user, and system time

    Raises:
        CommandError: If the command fails
    """
    if verbose:
        print("+", " ".join(str(part) for part in command))

    start = time.perf_counter()
    usage_before_self = resource.getrusage(resource.RUSAGE_SELF)
    usage_before_children = resource.getrusage(resource.RUSAGE_CHILDREN)

    try:
        subprocess.run(command, cwd=cwd, check=True, capture_output=not verbose)
    except subprocess.CalledProcessError as exc:
        end = time.perf_counter()
        usage_after_self = resource.getrusage(resource.RUSAGE_SELF)
        usage_after_children = resource.getrusage(resource.RUSAGE_CHILDREN)

        timing = create_timing(
            start,
            end,
            usage_before_self,
            usage_before_children,
            usage_after_self,
            usage_after_children,
        )

        # Print performance info even on failure
        from .timing import format_timing

        print(f"⏱️  Failed command timing: {format_timing('', timing)}", file=sys.stderr)

        raise CommandError(
            f"Command failed with exit code {exc.returncode}: {' '.join(str(p) for p in command)}"
        ) from exc

    end = time.perf_counter()
    usage_after_self = resource.getrusage(resource.RUSAGE_SELF)
    usage_after_children = resource.getrusage(resource.RUSAGE_CHILDREN)

    return create_timing(
        start,
        end,
        usage_before_self,
        usage_before_children,
        usage_after_self,
        usage_after_children,
    )
