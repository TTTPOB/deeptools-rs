"""
Test data path management.

Provides utilities for locating test data files within the project.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class TestDataPaths:
    """Container for test data directory paths."""

    repo_root: Path
    data_root: Path
    test_data_root: Path

    @property
    def heatmapper_dir(self) -> Path:
        """Path to the test_heatmapper test data directory."""
        return self.data_root

    @property
    def general_test_data(self) -> Path:
        """Path to the general test_data directory."""
        return self.test_data_root


def get_test_data_paths(repo_root: Path) -> TestDataPaths:
    """
    Get paths to test data directories.

    Args:
        repo_root: Root directory of the repository

    Returns:
        TestDataPaths with all relevant paths
    """
    return TestDataPaths(
        repo_root=repo_root,
        data_root=repo_root / "deeptools" / "deeptools" / "test" / "test_heatmapper",
        test_data_root=repo_root / "deeptools" / "deeptools" / "test" / "test_data",
    )


def validate_test_data_exists(paths: TestDataPaths) -> bool:
    """
    Validate that required test data directories exist.

    Args:
        paths: TestDataPaths to validate

    Returns:
        True if all paths exist, False otherwise
    """
    return paths.data_root.exists() and paths.test_data_root.exists()
