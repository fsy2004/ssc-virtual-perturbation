#!/bin/bash
set -euo pipefail
root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
runtime="$root/audit/postproduction_runtime"
pidfile="$runtime/rep01_primary_pbc.pid"
log="$runtime/rep01_primary_pbc.log"
if [[ -f "$pidfile" ]]; then
  pid=$(tr -d '[:space:]' < "$pidfile")
  ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,rss=,cmd= -p "$pid" || true
  pgrep -P "$pid" | xargs -r ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,rss=,cmd= -p || true
fi
printf '%s\n' '--- log tail ---'
tail -n 60 "$log" 2>/dev/null || true
printf '%s\n' '--- outputs ---'
find "$root/rep01" -maxdepth 3 -type f \( -name '*.xtc' -o -name '*.json' -o -name '*.xvg' \) -printf '%p %s bytes\n' 2>/dev/null | sort | tail -n 30
find "$root/analysis/trajectories/8kct_nirogacestat_native/rep01" -maxdepth 1 -type f -printf '%p %s bytes\n' 2>/dev/null | sort
printf '%s\n' '--- disk ---'
df -h /root/autodl-tmp
