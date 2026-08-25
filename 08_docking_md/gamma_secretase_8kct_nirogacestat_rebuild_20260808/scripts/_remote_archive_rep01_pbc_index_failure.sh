set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
FAILED="$ROOT/analysis/trajectories/8kct_nirogacestat_native/rep01"
STAMP=$(date +%Y%m%dT%H%M%S%z)
ARCHIVE="$ROOT/audit/postproduction_failures/${STAMP}_rep01_zero_based_analysis_index"
test -d "$ROOT"
test -d "$FAILED"
test -f "$ROOT/builds/analysis.ndx"
mkdir -p "$ARCHIVE"
cp -a "$ROOT/builds/analysis.ndx" "$ARCHIVE/analysis.zero_based.ndx"
cp -a "$ROOT/audit/postproduction_runtime/rep01_primary_pbc.log" "$ARCHIVE/rep01_primary_pbc.log" 2>/dev/null || true
mv "$FAILED" "$ARCHIVE/derived_trajectory_attempt"
sha256sum "$ARCHIVE/analysis.zero_based.ndx" > "$ARCHIVE/analysis.zero_based.ndx.sha256"
printf '%s\n' "$ARCHIVE"
find "$ARCHIVE" -maxdepth 2 -type f -printf '%p %s bytes\n' | sort
