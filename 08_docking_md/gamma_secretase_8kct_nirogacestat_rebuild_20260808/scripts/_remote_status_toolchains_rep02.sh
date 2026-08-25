#!/usr/bin/env bash
set -u

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
gmx_audit="$root/audit/toolchain_install/20260822_gmx_mmpbsa_retry_v5/20260822T112115Z"
gorder_audit="$root/audit/toolchain_install/20260822_gorder_retry_v4/20260822T112115Z"

echo '=== PROCESSES ==='
ps -eo pid,ppid,stat,etime,pcpu,pmem,rss,args | awk \
  'NR==1 || /prepare_primary_pbc_trajectories|gmx trjconv|install_gmx_mmpbsa|install_gorder|conda.*create|pip.*mpi4py|cargo.*install|rustc/ {print}'

echo '=== GMX_MMPBSA V5 ==='
tail -n 40 "$gmx_audit/gmx_mmpbsa_install.log" 2>/dev/null || true

echo '=== GORDER V4 ==='
tail -n 40 "$gorder_audit/gorder_install.log" 2>/dev/null || true

echo '=== REP02 ==='
tail -n 20 "$root/audit/postproduction_runtime/rep02_primary_pbc.log" 2>/dev/null || true
find "$root/analysis/trajectories/8kct_nirogacestat_native/rep02" -maxdepth 1 -type f \
  -printf '%f\t%s\t%TY-%Tm-%TdT%TH:%TM:%TSZ\n' 2>/dev/null | sort

echo '=== RESOURCES ==='
df -h /root/autodl-tmp
cat /sys/fs/cgroup/memory.events 2>/dev/null || true
exit 0
