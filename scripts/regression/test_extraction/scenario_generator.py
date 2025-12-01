"""
Test scenario loading and generation.

Loads test configurations from YAML and generates executable test scenarios.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import yaml


@dataclass
class TestScenario:
    """A single test scenario for validation."""

    name: str
    python_test: str
    reference_matrix: str
    rust_args: List[str]
    tolerance: float = 5e-6
    cross_validate: bool = False
    enabled: bool = True

    def get_rust_command(
        self,
        cargo: str,
        data_root: Path,
        test_data_root: Path,
        output_path: Path,
    ) -> List[str]:
        """Build the complete Rust command with substituted paths."""
        cmd = [cargo, "run", "--release", "--quiet", "--"]

        for arg in self.rust_args:
            # Substitute path placeholders
            substituted = arg
            substituted = substituted.replace("{data_root}", str(data_root))
            substituted = substituted.replace("{test_data_root}", str(test_data_root))
            substituted = substituted.replace("{output}", str(output_path))
            cmd.append(substituted)

        # Add output filename
        cmd.extend(["--outFileName", str(output_path)])

        return cmd


@dataclass
class TestSuite:
    """A collection of test scenarios."""

    name: str
    description: str
    data_root: str
    test_data_root: str
    scenarios: List[TestScenario] = field(default_factory=list)

    def get_data_root(self, repo_root: Path) -> Path:
        """Get the absolute data root path."""
        return repo_root / self.data_root

    def get_test_data_root(self, repo_root: Path) -> Path:
        """Get the absolute test data root path."""
        return repo_root / self.test_data_root


def load_test_config(config_path: Path) -> Dict[str, TestSuite]:
    """
    Load test configuration from a YAML file.

    Args:
        config_path: Path to the YAML configuration file

    Returns:
        Dictionary mapping suite names to TestSuite objects
    """
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    suites = {}
    for suite_name, suite_config in config.get("test_suites", {}).items():
        scenarios = []
        for scenario_config in suite_config.get("scenarios", []):
            scenario = TestScenario(
                name=scenario_config.get("name", "unnamed"),
                python_test=scenario_config.get("python_test", ""),
                reference_matrix=scenario_config.get("reference_matrix", ""),
                rust_args=scenario_config.get("rust_args", []),
                tolerance=scenario_config.get("tolerance", 5e-6),
                cross_validate=scenario_config.get("cross_validate", False),
                enabled=scenario_config.get("enabled", True),
            )
            scenarios.append(scenario)

        suite = TestSuite(
            name=suite_name,
            description=suite_config.get("description", ""),
            data_root=suite_config.get("data_root", ""),
            test_data_root=suite_config.get("test_data_root", ""),
            scenarios=scenarios,
        )
        suites[suite_name] = suite

    return suites


def load_reference_matrix(
    reference_name: str,
    data_root: Path,
) -> Path:
    """
    Locate a reference matrix file.

    Args:
        reference_name: Name of the reference matrix file
        data_root: Root directory containing test data

    Returns:
        Path to the reference matrix file

    Raises:
        FileNotFoundError: If the reference matrix is not found
    """
    reference_path = data_root / reference_name

    if not reference_path.exists():
        raise FileNotFoundError(
            f"Reference matrix not found: {reference_path}"
        )

    return reference_path


def get_default_config_path() -> Path:
    """Get the default configuration file path."""
    scripts_dir = Path(__file__).resolve().parents[2]
    return scripts_dir / "config" / "python_compatibility.yaml"
