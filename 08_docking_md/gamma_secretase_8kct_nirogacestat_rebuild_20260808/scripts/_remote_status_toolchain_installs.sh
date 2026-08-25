#!/bin/bash
set -euo pipefail

audit=/root/autodl-tmp/o6u_md_release_3x500ns_v4/audit/toolchain_install/20260822
for name in gmx_mmpbsa gorder fatslim; do
  pidfile="$audit/${name}_install.pid"
  log="$audit/${name}_install.log"
  printf '%s\n' "--- $name ---"
  if [[ -f "$pidfile" ]]; then
    pid=$(tr -d '[:space:]' < "$pidfile")
    ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,rss=,cmd= -p "$pid" || true
    pgrep -P "$pid" | xargs -r ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,rss=,cmd= -p || true
  else
    echo 'not launched'
  fi
  tail -n 20 "$log" 2>/dev/null || true
done
df -h /root/autodl-tmp
