#!/bin/bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
audit="$root/audit/toolchain_install/20260822"
prefix=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822
installer="$root/scripts/install_gmx_mmpbsa_1_6_5_cpu.sh"
log="$audit/gmx_mmpbsa_install.log"
pidfile="$audit/gmx_mmpbsa_install.pid"

mkdir -p "$audit" /root/autodl-tmp/envs
[[ -x /root/miniconda3/bin/conda ]]
[[ -f "$installer" ]]
[[ ! -e "$prefix" ]] || { echo "refusing existing prefix: $prefix" >&2; exit 3; }
[[ ! -e "$pidfile" ]] || { echo "refusing existing pidfile: $pidfile" >&2; exit 4; }

nohup env CONDA_EXE=/root/miniconda3/bin/conda \
  bash "$installer" "$prefix" >"$log" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > "$pidfile"
printf 'pid=%s\nprefix=%s\nlog=%s\n' "$pid" "$prefix" "$log"
