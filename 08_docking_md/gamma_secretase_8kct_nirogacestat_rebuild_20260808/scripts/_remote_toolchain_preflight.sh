#!/bin/bash
set -euo pipefail

date -u '+PRECHECK_UTC=%Y-%m-%dT%H:%M:%SZ'
hostname | sed 's/^/HOST=/'
nproc | sed 's/^/CPU=/'
free -h
df -h /root/autodl-tmp
nvidia-smi --query-gpu=name,temperature.gpu,utilization.gpu,memory.used,memory.total,power.draw --format=csv,noheader

printf '%s\n' '--- tools ---'
for tool in conda mamba micromamba git curl wget cargo rustc; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '%s=%s\n' "$tool" "$(command -v "$tool")"
  else
    printf '%s=MISSING\n' "$tool"
  fi
done

printf '%s\n' '--- PBC processes ---'
ps -eo pid,ppid,etime,%cpu,%mem,stat,args \
  | grep -E 'prepare_primary_pbc|gmx trjconv' \
  | grep -v grep || true

printf '%s\n' '--- install processes ---'
ps -eo pid,ppid,etime,%cpu,%mem,stat,args \
  | grep -E 'gmx_mmpbsa|gorder|FATSLiM|conda create' \
  | grep -v grep || true
