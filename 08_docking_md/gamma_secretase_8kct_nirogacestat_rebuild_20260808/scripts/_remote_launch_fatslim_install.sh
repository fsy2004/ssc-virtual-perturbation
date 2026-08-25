#!/bin/bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
audit="$root/audit/toolchain_install/20260822"
base=/root/autodl-tmp/tools/fatslim_0_2_2_20260822
installer="$root/scripts/install_fatslim_0_2_2.sh"
log="$audit/fatslim_install.log"
pidfile="$audit/fatslim_install.pid"

mkdir -p "$audit" /root/autodl-tmp/tools
[[ -f "$installer" ]]
[[ ! -e "$base" ]] || { echo "refusing existing prefix: $base" >&2; exit 3; }
[[ ! -e "$pidfile" ]] || { echo "refusing existing pidfile: $pidfile" >&2; exit 4; }

nohup bash "$installer" "$base" >"$log" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > "$pidfile"
printf 'pid=%s\nbase=%s\nlog=%s\n' "$pid" "$base" "$log"
