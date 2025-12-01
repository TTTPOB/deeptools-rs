"""
Timing and cached result dataclasses.

Extracted from compute_matrix_regression.py for modularity.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List


@dataclass
class CommandTiming:
    """Tracks wall, user, and system time for command execution."""

    wall_seconds: float
    user_seconds: float
    system_seconds: float

    def to_dict(self) -> dict:
        return {
            "wall_seconds": self.wall_seconds,
            "user_seconds": self.user_seconds,
            "system_seconds": self.system_seconds,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "CommandTiming":
        return cls(
            wall_seconds=data["wall_seconds"],
            user_seconds=data["user_seconds"],
            system_seconds=data["system_seconds"],
        )


@dataclass
class CachedResult:
    """Stores command hash, output path, timing, and timestamp for caching."""

    command_hash: str
    output_path: str
    timing: CommandTiming
    timestamp: float
    command: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "command_hash": self.command_hash,
            "output_path": self.output_path,
            "timing": self.timing.to_dict(),
            "timestamp": self.timestamp,
            "command": self.command,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "CachedResult":
        raw_cmd = data.get("command", [])
        if isinstance(raw_cmd, str):
            command = [raw_cmd]
        else:
            command = [str(part) for part in raw_cmd]
        return cls(
            command_hash=data["command_hash"],
            output_path=data["output_path"],
            timing=CommandTiming.from_dict(data["timing"]),
            timestamp=data["timestamp"],
            command=command,
        )


def calculate_resource_usage(
    usage_before_self,
    usage_before_children,
    usage_after_self,
    usage_after_children,
) -> tuple[float, float]:
    """Calculate user and system seconds from resource usage snapshots."""
    user_seconds = (usage_after_children.ru_utime - usage_before_children.ru_utime) + (
        usage_after_self.ru_utime - usage_before_self.ru_utime
    )
    system_seconds = (
        usage_after_children.ru_stime - usage_before_children.ru_stime
    ) + (usage_after_self.ru_stime - usage_before_self.ru_stime)
    return user_seconds, system_seconds


def create_timing(
    start: float,
    end: float,
    usage_before_self,
    usage_before_children,
    usage_after_self,
    usage_after_children,
) -> CommandTiming:
    """Create a CommandTiming object from resource usage snapshots."""
    user_seconds, system_seconds = calculate_resource_usage(
        usage_before_self, usage_before_children, usage_after_self, usage_after_children
    )
    return CommandTiming(
        wall_seconds=end - start,
        user_seconds=user_seconds,
        system_seconds=system_seconds,
    )


def format_timing(label: str, timing: CommandTiming, from_cache: bool = False) -> str:
    """Format timing information for display."""
    cache_marker = " [FROM CACHE]" if from_cache else ""
    return (
        f"{label}: "
        f"wall={timing.wall_seconds:.2f}s, "
        f"user={timing.user_seconds:.2f}s, "
        f"sys={timing.system_seconds:.2f}s"
        f"{cache_marker}"
    )
