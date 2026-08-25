set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
ARCHIVE="$ROOT/audit/postproduction_failures/20260822T151624+0800_rep01_zero_based_analysis_index"
TARGET="$ARCHIVE/derived_trajectory_attempt/01_whole.xtc"
RESOLVED=$(realpath "$TARGET")
case "$RESOLVED" in
  "$ARCHIVE"/*) ;;
  *) echo "REFUSE: target escaped failure archive: $RESOLVED" >&2; exit 75 ;;
esac
test -f "$RESOLVED"
sha256sum "$RESOLVED" > "$ARCHIVE/derived_trajectory_attempt/01_whole.xtc.sha256"
stat -c '%n %s bytes' "$RESOLVED" >> "$ARCHIVE/derived_trajectory_attempt/01_whole.xtc.sha256"
rm -- "$RESOLVED"
test ! -e "$RESOLVED"
echo "REMOVED_INVALID_DERIVED_XTC=$RESOLVED"
cat "$ARCHIVE/derived_trajectory_attempt/01_whole.xtc.sha256"
df -h "$ROOT"
