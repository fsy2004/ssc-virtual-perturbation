set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
find "$ROOT/audit/postproduction_failures" -maxdepth 1 -type d -name '*_rep01_incomplete_fragment_cluster_group' -printf '%T@ %p\n' | sort -nr | head -5
find "$ROOT/analysis/trajectories/8kct_nirogacestat_native" -maxdepth 1 -type d -printf '%p\n' 2>/dev/null || true
