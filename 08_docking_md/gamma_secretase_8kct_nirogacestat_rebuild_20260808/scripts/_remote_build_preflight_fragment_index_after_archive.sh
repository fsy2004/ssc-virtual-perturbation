set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
GMX=/root/GROMACS-2025.2/bin/gmx
TPR="$ROOT/rep01/work/production.tpr"
ARCHIVE="$ROOT/audit/postproduction_failures/20260822T153520+0800_rep01_incomplete_fragment_cluster_group"
WHOLE="$ARCHIVE/derived_trajectory_attempt/01_whole.xtc"
PREFLIGHT="$ARCHIVE/fragment_closed_cluster_preflight"

test -f "$TPR"
test -f "$WHOLE"
mkdir -p "$PREFLIGHT"
cd "$ROOT"
PYTHONPATH=scripts "$PY" scripts/build_analysis_ndx.py \
  --trajectory-topology "$TPR" \
  --output "$ROOT/builds/analysis.ndx"
printf '1\n0\n' | "$GMX" trjconv \
  -s "$TPR" -f "$WHOLE" -o "$PREFLIGHT/cluster_0_20ps.xtc" \
  -n "$ROOT/builds/analysis.ndx" -pbc cluster -e 20 \
  > "$PREFLIGHT/trjconv.stdout" 2> "$PREFLIGHT/trjconv.stderr"
"$GMX" check -f "$PREFLIGHT/cluster_0_20ps.xtc" \
  > "$PREFLIGHT/gmx_check.stdout" 2> "$PREFLIGHT/gmx_check.stderr"
sha256sum "$ROOT/builds/analysis.ndx" "$ROOT/builds/analysis.ndx.provenance.json" "$PREFLIGHT/cluster_0_20ps.xtc" > "$PREFLIGHT/SHA256SUMS"
cat "$ROOT/builds/analysis.ndx.provenance.json"
cat "$PREFLIGHT/SHA256SUMS"
tail -n 24 "$PREFLIGHT/gmx_check.stderr"
