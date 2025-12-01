#!/usr/bin/env python3
"""Run the computeMatrix Python-compatibility regression suite.

This entry point mirrors deepTools' `test_heatmapper.py` cases and ensures the
Rust implementation matches the Python reference byte-for-byte.
"""
from __future__ import annotations

import sys

from custom_compare import (
    CommandError,
    ensure_commands_available,
    parse_args,
    run_python_compatibility_mode,
)


def main() -> int:
    # Force python-compatibility mode to avoid mixing with custom dataset runs.
    args = parse_args()
    args.mode = "python-compatibility"
    ensure_commands_available(args)
    return run_python_compatibility_mode(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CommandError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
