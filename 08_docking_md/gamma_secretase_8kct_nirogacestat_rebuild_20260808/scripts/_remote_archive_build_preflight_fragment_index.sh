set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
GMX=/root/GROMACS-2025.2/bin/gmx
TPR="$ROOT/rep01/work/production.tpr"
FAILED="$ROOT/analysis/trajectories/8kct_nirogacestat_native/rep01"
STAMP=$(date +%Y%m%dT%H%M%S%z)
ARCHIVE="$ROOT/audit/postproduction_failures/${STAMP}_rep01_incomplete_fragment_cluster_group"

test -d "$ROOT"
test -d "$FAILED"
test -f "$FAILED/01_whole.xtc"
test -f "$ROOT/builds/analysis.ndx"
mkdir -p "$ARCHIVE"
cp -a "$ROOT/builds/analysis.ndx" "$ARCHIVE/analysis.protein_only_cluster.ndx"
mv "$FAILED" "$ARCHIVE/derived_trajectory_attempt"

cd "$ROOT"
PYTHONPATH=scripts "$PY" scripts/build_analysis_ndx.py \
  --trajectory-topology "$TPR" \
  --output "$ROOT/builds/analysis.ndx"

PREFLIGHT="$ARCHIVE/fragment_closed_cluster_preflight"
mkdir -p "$PREFLIGHT"
printf '1\n0\n' | "$GMX" trjconv \
  -s "$TPR" \
  -f "$ARCHIVE/derived_trajectory_attempt/01_whole.xtc" \
  -o "$PREFLIGHT/cluster_0_20ps.xtc" \
  -n "$ROOT/builds/analysis.ndx" \
  -pbc cluster -e 20 \
  > "$PREFLIGHT/trjconv.stdout" 2> "$PREFLIGHT/trjconv.stderr"
"$GMX" check -f "$PREFLIGHT/cluster_0_20ps.xtc" \
  > "$PREFLIGHT/gmx_check.stdout" 2> "$PREFLIGHT/gmx_check.stderr"
sha256sum \
  "$ROOT/builds/analysis.ndx" \
  "$ROOT/builds/analysis.ndx.provenance.json" \
  "$PREFLIGHT/cluster_0_20ps.xtc" \
  > "$PREFLIGHT/SHA256SUMS"

echo "ARCHIVE=$ARCHIVE"
cat "$ROOT/builds/analysis.ndx.provenance.json"
cat "$PREFLIGHT/SHA256SUMS"
tail -n 20 "$PREFLIGHT/gmx_check.stderr"
df -h "$ROOT"
