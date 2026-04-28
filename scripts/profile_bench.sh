#!/bin/bash
set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 <name> <target> <hot-path> -- <command...>"
    echo "Example: $0 p1-i64-div \"eliminate i128 div in write_matrix_value\" \"write_matrix_value -> rint + __divti3\" -- cargo run --release -- scale-regions ..."
    exit 1
fi

NAME="$1"
TARGET="$2"
HOT_PATH="$3"
shift 3

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

# /usr/bin/time -v
echo "## /usr/bin/time -v" >> "$REPORT"
echo '```' >> "$REPORT"
/usr/bin/time -v "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# perf stat
echo "## perf stat" >> "$REPORT"
echo '```' >> "$REPORT"
perf stat -d "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# heaptrack
echo "## heaptrack" >> "$REPORT"
echo '```' >> "$REPORT"
heaptrack "$@" 2>> "$REPORT"
echo '```' >> "$REPORT"

# Cleanup raw profiler files
rm -f perf.data heaptrack.*.gz 2>/dev/null

echo "Report written to $REPORT"
