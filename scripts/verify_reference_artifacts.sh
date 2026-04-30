#!/usr/bin/env bash
# verify_reference_artifacts.sh
#
# Reads scripts/config/reference_artifacts.yaml, regenerates every declared
# master_* artifact, compares it against the committed version, and checks
# for orphan master_* files not listed in the manifest.
#
# Exit code: 0 = all pass, 1 = any failure.

set -euo pipefail

# ── Resolve project root (directory containing pixi.toml) ──────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ ! -f "$PROJECT_ROOT/pixi.toml" ]]; then
  echo "ERROR: pixi.toml not found at $PROJECT_ROOT" >&2
  exit 1
fi

MANIFEST="$PROJECT_ROOT/scripts/config/reference_artifacts.yaml"
ARTIFACTS_DIR="$PROJECT_ROOT/tests/data"
TMPDIR_BASE="$(mktemp -d "${TMPDIR:-/tmp}/verify_artifacts.XXXXXXXXXX")"

cleanup() { rm -rf "$TMPDIR_BASE"; }
trap cleanup EXIT

PASS=0
FAIL=0
FAILURES=()

echo "=== Reference Artifact Verification ==="
echo "Project root : $PROJECT_ROOT"
echo "Manifest     : $MANIFEST"
echo "Artifacts dir: $ARTIFACTS_DIR"
echo ""

# ── Use Python to parse the YAML manifest and drive verification ───────
cd "$PROJECT_ROOT"

python3 - "$MANIFEST" "$ARTIFACTS_DIR" "$TMPDIR_BASE" "$PROJECT_ROOT" <<'PYEOF'
import sys, os, subprocess, tempfile, glob, shlex
import yaml

manifest_path = sys.argv[1]
artifacts_dir = sys.argv[2]
tmpdir_base   = sys.argv[3]
project_root  = sys.argv[4]

with open(manifest_path) as f:
    manifest = yaml.safe_load(f)

path_vars = manifest["path_variables"]
artifacts = manifest["artifacts"]

# Collect all declared output filenames for orphan check
declared_outputs = set()
for art in artifacts:
    declared_outputs.add(art["output"])

pass_count = 0
fail_count = 0
failures = []

def resolve_vars(s, extra_vars=None):
    """Replace {var} placeholders with path_variables values and extra_vars."""
    result = s
    for k, v in path_vars.items():
        result = result.replace("{" + k + "}", v)
    if extra_vars:
        for k, v in extra_vars.items():
            result = result.replace("{" + k + "}", v)
    return result

def compare_files(committed, generated, method):
    """Compare two files using the specified method. Returns (ok, detail)."""
    if not os.path.exists(committed):
        return False, f"committed file missing: {committed}"
    if not os.path.exists(generated):
        return False, f"generated file missing: {generated}"

    if method == "gzip_content":
        # Decompress both, then diff
        try:
            c_data = subprocess.run(
                ["gzip", "-cd", committed],
                capture_output=True, check=True
            ).stdout
            g_data = subprocess.run(
                ["gzip", "-cd", generated],
                capture_output=True, check=True
            ).stdout
        except subprocess.CalledProcessError as e:
            return False, f"gzip decompress failed: {e}"
        if c_data == g_data:
            return True, ""
        else:
            # Show first few lines of diff for debugging
            import difflib
            c_lines = c_data.decode("utf-8", errors="replace").splitlines(keepends=True)
            g_lines = g_data.decode("utf-8", errors="replace").splitlines(keepends=True)
            diff = list(difflib.unified_diff(
                c_lines[:50], g_lines[:50],
                fromfile="committed", tofile="generated", lineterm=""
            ))
            detail = "\n".join(diff[:30])
            return False, f"content differs (gzip_content):\n{detail}"
    elif method == "direct":
        try:
            result = subprocess.run(
                ["diff", "-q", committed, generated],
                capture_output=True
            )
            if result.returncode == 0:
                return True, ""
            else:
                # Show actual diff
                result2 = subprocess.run(
                    ["diff", "-u", committed, generated],
                    capture_output=True
                )
                detail = result2.stdout.decode("utf-8", errors="replace")[:2000]
                return False, f"content differs (direct):\n{detail}"
        except Exception as e:
            return False, f"diff failed: {e}"
    else:
        return False, f"unknown comparison method: {method}"

# ── Process artifacts that have commands (primary generators) ──────────
# Group: primary artifacts generate themselves and optionally a paired .tab
primary_artifacts = [a for a in artifacts if "command" in a]
secondary_artifacts = [a for a in artifacts if "generated_by" in a]

for art in primary_artifacts:
    output_name = art["output"]
    method = art["comparison_method"]
    cmd_template = art["command"]
    paired_tab = art.get("paired_tab")

    # Prepare temp output paths
    tmpdir = tempfile.mkdtemp(dir=tmpdir_base)
    mat_out = os.path.join(tmpdir, output_name)

    extra_vars = {"out.mat": mat_out}
    if paired_tab:
        tab_out = os.path.join(tmpdir, paired_tab)
        extra_vars["out.tab"] = tab_out

    cmd_resolved = resolve_vars(cmd_template, extra_vars)

    # Run generation command from project root
    print(f"  Generating {output_name} ...")
    result = subprocess.run(
        cmd_resolved,
        shell=True,
        cwd=project_root,
        capture_output=True
    )
    if result.returncode != 0:
        stderr_text = result.stderr.decode("utf-8", errors="replace")[:1000]
        print(f"  FAIL {output_name} — generation command failed (exit {result.returncode})")
        if stderr_text.strip():
            print(f"       stderr: {stderr_text}")
        fail_count += 1
        failures.append(output_name)
        continue

    # Compare primary .mat
    committed = os.path.join(artifacts_dir, output_name)
    ok, detail = compare_files(committed, mat_out, method)
    if ok:
        print(f"  PASS {output_name}")
        pass_count += 1
    else:
        print(f"  FAIL {output_name} — {detail}")
        fail_count += 1
        failures.append(output_name)

    # Compare paired .tab if present
    if paired_tab:
        committed_tab = os.path.join(artifacts_dir, paired_tab)
        tab_method = "direct"
        ok, detail = compare_files(committed_tab, tab_out, tab_method)
        if ok:
            print(f"  PASS {paired_tab}")
            pass_count += 1
        else:
            print(f"  FAIL {paired_tab} — {detail}")
            fail_count += 1
            failures.append(paired_tab)

# ── Orphan check ──────────────────────────────────────────────────────
print("")
print("--- Orphan Check ---")
master_files = glob.glob(os.path.join(artifacts_dir, "master_*"))
orphans = []
for f in sorted(master_files):
    basename = os.path.basename(f)
    if basename not in declared_outputs:
        orphans.append(basename)

if orphans:
    print(f"  FAIL — {len(orphans)} orphan(s) found not declared in manifest:")
    for o in orphans:
        print(f"    - {o}")
    fail_count += 1
    failures.append("orphan_check")
else:
    print(f"  PASS — all {len(master_files)} master_* files are declared in manifest")
    pass_count += 1

# ── Summary ───────────────────────────────────────────────────────────
print("")
print("=== Summary ===")
print(f"  PASS: {pass_count}")
print(f"  FAIL: {fail_count}")
if failures:
    print(f"  Failed items: {', '.join(failures)}")
    sys.exit(1)
else:
    print("  All checks passed.")
    sys.exit(0)
PYEOF
