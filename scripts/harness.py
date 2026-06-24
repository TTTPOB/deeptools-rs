#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = REPO_ROOT / "scripts" / "configs"


class HarnessError(RuntimeError):
    pass


@dataclass(frozen=True)
class Case:
    raw: dict[str, Any]
    manifest: "Manifest"

    @property
    def id(self) -> str:
        return self.raw["id"]

    @property
    def tags(self) -> list[str]:
        return list(self.raw.get("tags", []))

    @property
    def mode(self) -> str:
        return self.raw["mode"]

    @property
    def reference(self) -> dict[str, str]:
        return dict(self.raw.get("reference") or {})

    @property
    def dataset(self) -> str | None:
        return self.raw.get("dataset")

    def has_tag(self, tag: str) -> bool:
        return tag in self.tags

    def command(
        self,
        runner: str,
        output: Path,
        *,
        matrix_values: Path | None = None,
        sorted_regions: Path | None = None,
        release: bool = False,
    ) -> list[str]:
        cmd = self.manifest.runner_prefix(runner, self.mode, release=release)
        cmd.extend(["-R", *self.resolved_paths("region_files")])
        cmd.extend(["-S", *self.resolved_paths("score_files")])
        cmd.extend(self.resolved_options())
        cmd.extend(["--outFileName", str(output)])
        if matrix_values is not None:
            cmd.extend(["--outFileNameMatrix", str(matrix_values)])
        if sorted_regions is not None:
            cmd.extend(["--outFileSortedRegions", str(sorted_regions)])
        return cmd

    def resolved_paths(self, key: str) -> list[str]:
        return [str(self.manifest.resolve_path(value)) for value in self.raw.get(key, [])]

    def resolved_options(self) -> list[str]:
        return [self.manifest.resolve_arg(value) for value in self.raw.get("options", [])]


@dataclass(frozen=True)
class Manifest:
    raw: dict[str, Any]
    repo_root: Path

    @classmethod
    def load(cls, *config_names: str) -> "Manifest":
        with (CONFIG_DIR / "common.json").open("r", encoding="utf-8") as handle:
            raw = json.load(handle)
        raw.setdefault("cases", [])
        raw.setdefault("datasets", {})
        for name in config_names:
            for path in config_paths(name):
                with path.open("r", encoding="utf-8") as handle:
                    data = json.load(handle)
                raw["cases"].extend(data.get("cases", []))
                raw["datasets"].update(data.get("datasets", {}))
        return cls(raw, REPO_ROOT)

    @property
    def cases(self) -> list[Case]:
        return [Case(raw=case, manifest=self) for case in self.raw["cases"]]

    @property
    def default_tolerance(self) -> float:
        return float(self.raw["comparison"]["default_tolerance"])

    @property
    def ignore_header_keys(self) -> list[str]:
        return list(self.raw["comparison"]["ignore_header_keys"])

    @property
    def max_diffs(self) -> int:
        return int(self.raw["comparison"]["max_diffs"])

    def case(self, case_id: str) -> Case:
        for case in self.cases:
            if case.id == case_id:
                return case
        raise HarnessError(f"unknown case: {case_id}")

    def select_cases(
        self,
        *,
        tag: str | None = None,
        case_id: str | None = None,
    ) -> list[Case]:
        cases = self.cases
        if tag:
            cases = [case for case in cases if case.has_tag(tag)]
        if case_id:
            cases = [case for case in cases if case.id == case_id]
        return cases

    def resolve_path(self, value: str) -> Path:
        return self.repo_root / self.resolve_arg(value)

    def resolve_arg(self, value: str) -> str:
        resolved = value
        for key, replacement in self.raw["paths"].items():
            resolved = resolved.replace(f"{{{key}}}", replacement)
        return str(self.repo_root / resolved) if resolved.startswith((".", "/")) else resolved

    def runner_prefix(self, runner: str, mode: str, *, release: bool = False) -> list[str]:
        if runner == "rust":
            return [str(binary_path("compute_matrix_rs", release=release)), mode]
        if runner == "python":
            return ["pixi", "run", "computeMatrix", mode]
        raise HarnessError(f"unknown runner: {runner}")

    def dataset(self, name: str) -> dict[str, Any]:
        try:
            return dict(self.raw["datasets"][name])
        except KeyError as exc:
            raise HarnessError(f"unknown dataset: {name}") from exc


def binary_path(name: str, *, release: bool = False) -> Path:
    release_path = REPO_ROOT / "target" / "release" / name
    debug = REPO_ROOT / "target" / "debug" / name
    if release:
        return release_path
    if debug.exists():
        return debug
    if release_path.exists():
        return release_path
    return release_path


def config_paths(name: str) -> list[Path]:
    path = CONFIG_DIR / name
    if path.is_dir():
        return sorted(path.glob("*.json"))
    if path.suffix != ".json":
        path = path.with_suffix(".json")
    return [path]


def run(cmd: list[str], *, cwd: Path = REPO_ROOT, quiet: bool = False) -> float:
    if not quiet:
        print("+", " ".join(cmd))
    start = time.perf_counter()
    completed = subprocess.run(cmd, cwd=cwd)
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        raise HarnessError(f"command failed with exit {completed.returncode}: {' '.join(cmd)}")
    return elapsed


def compare(
    manifest: Manifest,
    candidate: Path,
    reference: Path,
    *,
    tolerance: float | None = None,
    quiet: bool = False,
) -> None:
    cmd = [
        str(binary_path("compare_matrix")),
        "diff",
        str(candidate),
        str(reference),
        "--tolerance",
        str(tolerance or manifest.default_tolerance),
        "--max-diffs",
        str(manifest.max_diffs),
    ]
    for key in manifest.ignore_header_keys:
        cmd.extend(["--ignore", key])
    run(cmd, quiet=quiet)


def command_compat(args: argparse.Namespace) -> int:
    manifest = Manifest.load("compat")
    ensure_binaries()
    cases = manifest.select_cases(tag="compat", case_id=args.case)
    if not cases:
        raise HarnessError("no cases selected")

    output_dir = args.output_dir or REPO_ROOT / "target" / "harness" / "compat"
    output_dir.mkdir(parents=True, exist_ok=True)
    failures = 0

    for case in cases:
        output = output_dir / f"{case.id}.mat.gz"
        tab = output_dir / f"{case.id}.tab" if "matrix_values" in case.reference else None
        bed = output_dir / f"{case.id}.bed" if "sorted_regions" in case.reference else None
        print(f"[compat] {case.id}")
        run(case.command("rust", output, matrix_values=tab, sorted_regions=bed), quiet=args.quiet)
        compare(
            manifest,
            output,
            manifest.resolve_path(case.reference["matrix"]),
            quiet=args.quiet,
        )
        if tab is not None:
            compare_text(tab, manifest.resolve_path(case.reference["matrix_values"]))
        if bed is not None:
            compare_text(bed, manifest.resolve_path(case.reference["sorted_regions"]))

    print(f"[compat] {len(cases) - failures}/{len(cases)} passed")
    return 0


def command_regen_refs(args: argparse.Namespace) -> int:
    manifest = Manifest.load("artifacts.json", "compat/sweep.json")
    cases = manifest.select_cases(case_id=args.case)
    cases = [case for case in cases if case.reference]
    if not cases:
        raise HarnessError("no reference-generating cases selected")

    for case in cases:
        matrix = manifest.resolve_path(case.reference["matrix"])
        matrix.parent.mkdir(parents=True, exist_ok=True)
        tab = (
            manifest.resolve_path(case.reference["matrix_values"])
            if "matrix_values" in case.reference
            else None
        )
        bed = (
            manifest.resolve_path(case.reference["sorted_regions"])
            if "sorted_regions" in case.reference
            else None
        )
        if tab is not None:
            tab.parent.mkdir(parents=True, exist_ok=True)
        if bed is not None:
            bed.parent.mkdir(parents=True, exist_ok=True)
        print(f"[regen-ref] {case.id}")
        run(case.command("python", matrix, matrix_values=tab, sorted_regions=bed), quiet=args.quiet)
    return 0


def command_verify_refs(args: argparse.Namespace) -> int:
    manifest = Manifest.load("artifacts.json", "compat/sweep.json")
    ensure_binaries()
    cases = manifest.select_cases(case_id=args.case)
    cases = [case for case in cases if case.reference]
    if not cases:
        raise HarnessError("no reference cases selected")

    with tempfile.TemporaryDirectory(prefix="compute-matrix-verify-") as tmp:
        tmpdir = Path(tmp)
        for case in cases:
            print(f"[verify-ref] {case.id}")
            output = tmpdir / f"{case.id}.mat.gz"
            run(case.command("rust", output), quiet=args.quiet)
            compare(
                manifest,
                output,
                manifest.resolve_path(case.reference["matrix"]),
                quiet=args.quiet,
            )
            if "matrix_values" in case.reference:
                tab = tmpdir / f"{case.id}.tab"
                run(case.command("rust", output, matrix_values=tab), quiet=True)
                compare_text(tab, manifest.resolve_path(case.reference["matrix_values"]))
            if "sorted_regions" in case.reference:
                bed = tmpdir / f"{case.id}.bed"
                run(case.command("rust", output, sorted_regions=bed), quiet=True)
                compare_text(bed, manifest.resolve_path(case.reference["sorted_regions"]))
    return 0


def command_bench_smoke(args: argparse.Namespace) -> int:
    manifest = Manifest.load("benchmarks.json")
    ensure_binaries(release=True)
    cases = manifest.select_cases(tag=args.tag, case_id=args.case)
    if not cases:
        raise HarnessError("no benchmark cases selected")

    output_dir = args.output_dir or REPO_ROOT / "target" / "harness" / "bench-smoke"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary: list[dict[str, Any]] = []

    for case in cases:
        output = output_dir / f"{case.id}.mat.gz"
        print(f"[bench-smoke] {case.id}")
        if args.warmup:
            warmup = output_dir / f"{case.id}.warmup.mat.gz"
            run(case.command("rust", warmup, release=True), quiet=True)
        elapsed = run(case.command("rust", output, release=True), quiet=args.quiet)
        summary.append({"case": case.id, "elapsed_seconds": elapsed})
        print(f"[bench-smoke] {case.id}: {elapsed:.3f}s")

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[bench-smoke] wrote {summary_path}")
    return 0


def command_prepare_data(args: argparse.Namespace) -> int:
    manifest = Manifest.load("datasets.json")
    if args.dataset == "all":
        names = sorted(manifest.raw.get("datasets", {}).keys())
    else:
        names = [args.dataset]
    for name in names:
        prepare_dataset(manifest, name, verbose=not args.quiet)
    return 0


def command_encode(args: argparse.Namespace) -> int:
    manifest = Manifest.load("benchmarks.json", "datasets.json")
    prepare_dataset(manifest, "encode_k562_atac", verbose=not args.quiet)
    ensure_binaries(release=True)
    cases = manifest.select_cases(tag="encode", case_id=args.case)
    output_dir = args.output_dir or REPO_ROOT / "target" / "harness" / "encode"
    output_dir.mkdir(parents=True, exist_ok=True)

    for case in cases:
        output = output_dir / f"{case.id}.mat.gz"
        print(f"[encode] rust {case.id}")
        rust_time = run(case.command("rust", output, release=True), quiet=args.quiet)
        print(f"[encode] rust {case.id}: {rust_time:.3f}s")
        if args.cross_validate:
            ref = output_dir / f"{case.id}.python.mat.gz"
            print(f"[encode] python {case.id}")
            python_time = run(case.command("python", ref), quiet=args.quiet)
            compare(manifest, output, ref, quiet=args.quiet)
            print(f"[encode] speedup {case.id}: {python_time / rust_time:.2f}x")
    return 0


def command_profile(args: argparse.Namespace) -> int:
    manifest = Manifest.load("benchmarks.json")
    ensure_binaries(release=True)
    case = manifest.case(args.case)
    output_dir = args.output_dir or REPO_ROOT / "bench_reports"
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d-%H%M%S")
    report = output_dir / f"{timestamp}-{args.name or case.id}.md"
    matrix_output = REPO_ROOT / "target" / "harness" / "profile" / f"{case.id}.mat.gz"
    matrix_output.parent.mkdir(parents=True, exist_ok=True)
    cmd = case.command("rust", matrix_output, release=True)

    if args.warmup:
        warmup = matrix_output.with_suffix(".warmup.mat.gz")
        run(case.command("rust", warmup, release=True), quiet=True)

    lines = [
        f"# Profile: {args.name or case.id}",
        f"Time: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
        f"Case: {case.id}",
        f"Command: {' '.join(cmd)}",
        "",
    ]

    lines.extend(run_capture_section("time -v", ["/usr/bin/time", "-v", *cmd]))

    if shutil.which("perf"):
        lines.extend(run_capture_section("perf stat", ["perf", "stat", "-d", *cmd]))
        perf_data = Path(tempfile.mktemp(suffix=".perf.data"))
        try:
            subprocess.run(
                [
                    "perf",
                    "record",
                    "-e",
                    "cpu-clock",
                    "-g",
                    "--call-graph",
                    "dwarf",
                    "-o",
                    str(perf_data),
                    *cmd,
                ],
                cwd=REPO_ROOT,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            if perf_data.exists():
                lines.extend(
                    run_capture_section(
                        "CPU Hotspots",
                        [
                            "perf",
                            "report",
                            "-i",
                            str(perf_data),
                            "--stdio",
                            "--no-children",
                            "--percent-limit",
                            "0.5",
                        ],
                        limit_lines=80,
                    )
                )
        finally:
            perf_data.unlink(missing_ok=True)
    else:
        lines.extend(["## perf", "", "`perf` not found.", ""])

    if shutil.which("heaptrack"):
        heap_dir = Path(tempfile.mkdtemp(prefix="compute-matrix-heaptrack-"))
        heap_prefix = heap_dir / "heaptrack"
        try:
            subprocess.run(
                ["heaptrack", "-o", str(heap_prefix), *cmd],
                cwd=REPO_ROOT,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            heap_file = next(heap_dir.glob("heaptrack*.zst"), None) or next(
                heap_dir.glob("heaptrack*.gz"), None
            )
            if heap_file and shutil.which("heaptrack_print"):
                lines.extend(
                    run_capture_section(
                        "heaptrack summary",
                        ["heaptrack_print", "-f", str(heap_file)],
                        limit_lines=80,
                    )
                )
        finally:
            shutil.rmtree(heap_dir, ignore_errors=True)
    else:
        lines.extend(["## heaptrack", "", "`heaptrack` not found.", ""])

    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"[profile] wrote {report}")
    return 0


def run_capture_section(
    title: str,
    cmd: list[str],
    *,
    limit_lines: int | None = None,
) -> list[str]:
    completed = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = completed.stdout.splitlines()
    if limit_lines is not None:
        output = output[:limit_lines]
    return [f"## {title}", "```", *output, "```", ""]


def compare_text(left: Path, right: Path) -> None:
    if left.read_bytes() != right.read_bytes():
        raise HarnessError(f"text output differs: {left} vs {right}")


def ensure_binaries(*, release: bool = False) -> None:
    names = ["compute_matrix_rs", "compare_matrix"]
    missing = [name for name in names if not binary_path(name).exists()]
    if missing or release:
        cmd = ["cargo", "build"]
        if release:
            cmd.append("--release")
        run(cmd)


def prepare_dataset(manifest: Manifest, name: str, *, verbose: bool) -> None:
    dataset = manifest.dataset(name)
    root = manifest.resolve_path(dataset["root"])
    root.mkdir(parents=True, exist_ok=True)

    for signal in dataset.get("signals", []):
        download(manifest.resolve_path(signal["path"]), signal["url"], verbose=verbose)

    bed = dataset.get("bed")
    if bed:
        compressed = manifest.resolve_path(bed["compressed_path"])
        plain = manifest.resolve_path(bed["path"])
        download(compressed, bed["url"], verbose=verbose)
        if not plain.exists() or plain.stat().st_mtime < compressed.stat().st_mtime:
            if verbose:
                print(f"[prepare-data] decompress {compressed} -> {plain}")
            with gzip.open(compressed, "rb") as src, plain.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        split_bed(plain, manifest.resolve_path(bed["split_dir"]), verbose=verbose)


def download(destination: Path, url: str, *, verbose: bool) -> None:
    if destination.exists():
        if verbose:
            print(f"[prepare-data] exists {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    if verbose:
        print(f"[prepare-data] download {url} -> {destination}")
    try:
        with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    tmp.rename(destination)


def split_bed(source: Path, target_dir: Path, *, verbose: bool) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    comments = [line for line in lines if line.startswith("#")]
    records = [line for line in lines if line.strip() and not line.startswith("#")]
    if not records:
        raise HarnessError(f"BED file has no records: {source}")
    midpoint = (len(records) + 1) // 2
    for suffix, subset in [("part1", records[:midpoint]), ("part2", records[midpoint:])]:
        if not subset:
            continue
        output = target_dir / f"{source.stem}_{suffix}.bed"
        if verbose:
            print(f"[prepare-data] write {output} ({len(subset)} records)")
        output.write_text("".join(comments + subset), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Unified computeMatrix test and benchmark harness")
    sub = parser.add_subparsers(dest="command", required=True)

    compat = sub.add_parser("compat", help="run Rust output against committed reference matrices")
    compat.add_argument("--case")
    compat.add_argument("--output-dir", type=Path)
    compat.add_argument("--quiet", action="store_true")
    compat.set_defaults(func=command_compat)

    regen = sub.add_parser("regen-refs", help="regenerate reference matrices with deeptools")
    regen.add_argument("--case")
    regen.add_argument("--quiet", action="store_true")
    regen.set_defaults(func=command_regen_refs)

    verify = sub.add_parser("verify-refs", help="verify Rust output against reference matrices")
    verify.add_argument("--case")
    verify.add_argument("--quiet", action="store_true")
    verify.set_defaults(func=command_verify_refs)

    bench = sub.add_parser("bench-smoke", help="run small release-mode performance smoke cases")
    bench.add_argument("--case")
    bench.add_argument("--tag", default="perf_smoke")
    bench.add_argument("--output-dir", type=Path)
    bench.add_argument("--warmup", action=argparse.BooleanOptionalAction, default=True)
    bench.add_argument("--quiet", action="store_true")
    bench.set_defaults(func=command_bench_smoke)

    prepare = sub.add_parser("prepare-data", help="download and prepare benchmark datasets")
    prepare.add_argument("dataset", nargs="?", default="all")
    prepare.add_argument("--quiet", action="store_true")
    prepare.set_defaults(func=command_prepare_data)

    encode = sub.add_parser("encode", help="run ENCODE performance cases")
    encode.add_argument("--case")
    encode.add_argument("--output-dir", type=Path)
    encode.add_argument("--cross-validate", action="store_true")
    encode.add_argument("--quiet", action="store_true")
    encode.set_defaults(func=command_encode)

    profile = sub.add_parser("profile", help="profile one manifest case with time/perf/heaptrack")
    profile.add_argument("case")
    profile.add_argument("--name")
    profile.add_argument("--output-dir", type=Path)
    profile.add_argument("--warmup", action=argparse.BooleanOptionalAction, default=True)
    profile.set_defaults(func=command_profile)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except HarnessError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
