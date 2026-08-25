#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
bash "$root/scripts/_remote_retry_gmx_mmpbsa_v5.sh"
bash "$root/scripts/_remote_retry_gorder_v4.sh"

echo '=== PROCESS SNAPSHOT ==='
ps -eo pid,ppid,stat,etime,pcpu,pmem,rss,args | awk \
  'NR==1 || /prepare_primary_pbc_trajectories|gmx trjconv|install_gmx_mmpbsa|install_gorder|conda.*create|cargo.*install/ {print}'

echo '=== REP02 LOG TAIL ==='
tail -n 20 "$root/audit/postproduction_runtime/rep02_primary_pbc.log" || true

echo '=== REP02 FILES ==='
find "$root/analysis/trajectories/8kct_nirogacestat_native/rep02" -maxdepth 1 -type f \
  -printf '%f\t%s\t%TY-%Tm-%TdT%TH:%TM:%TSZ\n' | sort

echo '=== FATSLIM V4 RECEIPTS ==='
fatslim=/root/autodl-tmp/tools/fatslim_0_2_2_20260822_v4
find "$fatslim" -maxdepth 1 -type f -printf '%f\t%s\n' | sort
if [[ -x "$fatslim/env/bin/fatslim" ]]; then
  "$fatslim/env/bin/fatslim" version
fi
for digest in "$fatslim"/*.sha256; do
  [[ -f "$digest" ]] || continue
  (cd "$(dirname "$digest")" && sha256sum -c "$(basename "$digest")")
done

echo '=== RESOURCES ==='
free -h
df -h /root/autodl-tmp
for f in /sys/fs/cgroup/memory.events /sys/fs/cgroup/memory/memory.oom_control; do
  [[ -r "$f" ]] && { echo "$f"; cat "$f"; }
done
