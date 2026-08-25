#!/bin/bash
set -euo pipefail

path=/root/autodl-tmp/tools/gorder_1_5_0_20260822/rustup-init
if [[ -e "$path" ]]; then
  stat --printf='rustup_init_bytes=%s\nmtime=%y\n' "$path"
else
  echo 'rustup_init=MISSING'
fi
pidfile=/root/autodl-tmp/o6u_md_release_3x500ns_v4/audit/toolchain_install/20260822/gorder_install.pid
if [[ -f "$pidfile" ]]; then
  parent=$(tr -d '[:space:]' < "$pidfile")
  ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,wchan=,cmd= -p "$parent" || true
  pgrep -P "$parent" | xargs -r ps -o pid=,ppid=,etime=,stat=,%cpu=,%mem=,wchan=,cmd= -p || true
fi
