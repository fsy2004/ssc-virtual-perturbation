set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
REC="$ROOT/analysis/trajectories/8kct_nirogacestat_native/rep01/02_cluster_complex_if_required.command.json"
sed -n '1,240p' "$REC"
printf '\n--- index head ---\n'
sed -n '1,24p' "$ROOT/builds/analysis.ndx"
printf '\n--- index zero hits ---\n'
grep -nE '(^|[[:space:]])0([[:space:]]|$)' "$ROOT/builds/analysis.ndx" | head -20 || true
