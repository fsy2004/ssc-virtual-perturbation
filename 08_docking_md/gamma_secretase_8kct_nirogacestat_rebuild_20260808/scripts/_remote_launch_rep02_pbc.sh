#!/bin/bash
set -euo pipefail
root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
runtime="$root/audit/postproduction_runtime"
log="$runtime/rep02_primary_pbc.log"
pidfile="$runtime/rep02_primary_pbc.pid"
mkdir -p "$runtime"
active_pattern='^/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python scripts/prepare_primary_pbc_trajectories.py '
if pgrep -f "$active_pattern" >/dev/null; then
  echo 'REFUSE: another primary PBC process is active' >&2
  pgrep -af "$active_pattern" >&2
  exit 75
fi
if [[ -e "$pidfile" ]]; then
  oldpid=$(tr -d '[:space:]' < "$pidfile")
  if [[ "$oldpid" =~ ^[0-9]+$ ]] && kill -0 "$oldpid" 2>/dev/null; then
    echo "REFUSE: recorded rep02 PBC PID is still active: $oldpid" >&2
    exit 75
  fi
fi
cd "$root"
nohup env PYTHONUNBUFFERED=1 /root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python \
  scripts/prepare_primary_pbc_trajectories.py \
  --release-root "$root" \
  --replica rep02 \
  --gmx /root/GROMACS-2025.2/bin/gmx \
  >"$log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "$pidfile"
sleep 1
if ! kill -0 "$pid" 2>/dev/null; then
  echo 'LAUNCH FAILED; log follows' >&2
  tail -n 80 "$log" >&2 || true
  exit 1
fi
echo "LAUNCHED rep02 primary PBC pid=$pid log=$log"
