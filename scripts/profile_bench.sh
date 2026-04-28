#!/bin/bash
set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 <name> <target> <hot-path> -- <command...>"
    echo "Example: $0 p1-i64-div \"eliminate i128 div\" \"write_matrix_value\" -- target/release/compute_matrix_rs scale-regions ..."
    exit 1
fi

NAME="$1"
TARGET="$2"
HOT_PATH="$3"
shift 3

if [ "$1" = "--" ]; then
    shift
fi

# Reject cargo run — heaptrack/time will include cargo overhead
CMD_BASE=$(basename "$1" 2>/dev/null || echo "$1")
if [ "$CMD_BASE" = "cargo" ]; then
    echo "ERROR: Running via 'cargo run' includes cargo overhead in profiling data." >&2
    echo "       Build first: cargo build --release" >&2
    echo "       Then use:    target/release/<binary> directly." >&2
    exit 1
fi

TIMESTAMP=$(date +%Y-%m-%d-%H%M%S)
REPORT="bench_reports/${TIMESTAMP}-${NAME}.md"
mkdir -p bench_reports

cat > "$REPORT" << EOF
# Profile: ${NAME}
Time: $(date -Iseconds)
Command: $@
Target: ${TARGET}
Hot path: ${HOT_PATH}

EOF

echo "=== Run 1/4: /usr/bin/time -v ===" >&2
echo "## /usr/bin/time -v" >> "$REPORT"
echo '```' >> "$REPORT"
/usr/bin/time -v "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

echo "=== Run 2/4: perf stat ===" >&2
echo "## perf stat" >> "$REPORT"
echo '```' >> "$REPORT"
perf stat -d "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

echo "=== Run 3/4: perf record (cpu-clock sampling) ===" >&2
PERF_DATA=$(mktemp --suffix=.perf.data)
perf record -e cpu-clock -g --call-graph dwarf -o "$PERF_DATA" "$@" 2>/dev/null

echo "## CPU Hotspots (perf report, top functions by self time)" >> "$REPORT"
echo '```' >> "$REPORT"
perf report -i "$PERF_DATA" --stdio --no-children --percent-limit 0.5 2>/dev/null \
    | grep -E '^\s+[0-9]' \
    | head -20 >> "$REPORT" || echo "(no samples captured)" >> "$REPORT"
echo '```' >> "$REPORT"

echo "" >> "$REPORT"
echo "## CPU Call Graph (perf report, top call chains)" >> "$REPORT"
echo '```' >> "$REPORT"
perf report -i "$PERF_DATA" --stdio --percent-limit 2 2>/dev/null \
    | head -120 >> "$REPORT" || echo "(no samples captured)" >> "$REPORT"
echo '```' >> "$REPORT"

rm -f "$PERF_DATA"

echo "=== Run 4/4: heaptrack ===" >&2
HEAPTRACK_PREFIX=$(mktemp -d)/heaptrack
heaptrack -o "$HEAPTRACK_PREFIX" "$@" 2>/dev/null

echo "" >> "$REPORT"
echo "## heaptrack summary" >> "$REPORT"
echo '```' >> "$REPORT"
# Find the actual heaptrack output file
HEAPTRACK_FILE=$(ls -1t "${HEAPTRACK_PREFIX}"*.zst "${HEAPTRACK_PREFIX}"*.gz 2>/dev/null | head -1)
if [ -n "$HEAPTRACK_FILE" ]; then
    heaptrack_print -f "$HEAPTRACK_FILE" 2>/dev/null \
        | grep -E '^(peak heap|total memory|calls to|temporary)' \
        >> "$REPORT" || true
    echo '```' >> "$REPORT"

    echo "" >> "$REPORT"
    echo "## Allocation Hotspots (heaptrack, top 10 by peak consumption)" >> "$REPORT"
    echo '```' >> "$REPORT"
    heaptrack_print -f "$HEAPTRACK_FILE" -p -n 10 2>/dev/null \
        | head -80 >> "$REPORT" || echo "(heaptrack_print failed)" >> "$REPORT"
    echo '```' >> "$REPORT"

    echo "" >> "$REPORT"
    echo "## Temporary Allocation Hotspots (heaptrack, top 10)" >> "$REPORT"
    echo '```' >> "$REPORT"
    heaptrack_print -f "$HEAPTRACK_FILE" -T -n 10 2>/dev/null \
        | head -80 >> "$REPORT" || echo "(heaptrack_print failed)" >> "$REPORT"
    echo '```' >> "$REPORT"

    rm -f "$HEAPTRACK_FILE"
else
    echo "(heaptrack output file not found)" >> "$REPORT"
    echo '```' >> "$REPORT"
fi

rm -rf "$(dirname "$HEAPTRACK_PREFIX")" 2>/dev/null || true

# Cleanup stray profiler files in cwd
rm -f perf.data heaptrack.*.gz heaptrack.*.zst 2>/dev/null

echo "" >> "$REPORT"
echo "## Agent TODO: Hotspot Analysis" >> "$REPORT"
cat >> "$REPORT" << 'AGENTEOF'
<!--
Agent: read the CPU Hotspots and Allocation Hotspots sections above.
Write a structured analysis below covering:
1. Comparison vs baseline/previous report (wall clock, RSS, task-clock, ctx switches)
2. Top 5 remaining CPU hotspots and whether they are targets for optimization
3. Top 3 remaining allocation hotspots and whether they can be reduced
4. Verdict: PASS (improvement or regression ≤5%) / FAIL (regression >5%)
-->
AGENTEOF

echo "" >&2
echo "======================================" >&2
echo "Report written to $REPORT" >&2
echo "======================================" >&2
echo "" >&2
echo "Next: read the report and fill in the 'Agent TODO: Hotspot Analysis' section." >&2
