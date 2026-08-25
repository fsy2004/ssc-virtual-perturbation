set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
ARCHIVE="$ROOT/audit/postproduction_failures/20260822T151624+0800_rep01_zero_based_analysis_index"
echo '--- matching processes ---'
pgrep -af 'sha256sum|prepare_primary_pbc_trajectories.py' || true
echo '--- invalid output evidence ---'
find "$ARCHIVE/derived_trajectory_attempt" -maxdepth 1 -type f -printf '%f %s bytes\n' | sort
echo '--- active output directory ---'
if test -d "$ROOT/analysis/trajectories/8kct_nirogacestat_native/rep01"; then
  find "$ROOT/analysis/trajectories/8kct_nirogacestat_native/rep01" -maxdepth 1 -type f -printf '%f %s bytes\n' | sort
else
  echo MISSING
fi
df -h "$ROOT"
