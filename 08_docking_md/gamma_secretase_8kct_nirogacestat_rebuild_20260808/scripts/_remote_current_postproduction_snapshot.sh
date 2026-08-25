#!/usr/bin/env bash
set -u

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
traj="$root/analysis/trajectories/8kct_nirogacestat_native"
gorder_audit="$root/audit/toolchain_install/20260822_gorder_retry_v5/20260822T113333Z"
gmx_audit="$root/audit/toolchain_install/20260822_gmx_mmpbsa_retry_v5/20260822T112115Z"

echo '=== FROZEN HASHES ==='
sha256sum \
  /root/autodl-tmp/o6u_md_release_3x500ns_v4.tgz \
  "$root/RELEASE_MANIFEST.json" \
  "$root/CANARY_VALIDATION.json" \
  "$root/PRODUCTION_TPR_RELEASE.json" \
  "$root/rep01/work/production.tpr" \
  "$root/rep02/work/production.tpr" \
  "$root/rep03/work/production.tpr"

echo '=== PROCESSES ==='
ps -eo pid,ppid,stat,etime,pcpu,pmem,rss,args | awk \
  'NR==1 || /prepare_primary_pbc_trajectories|gmx trjconv|gmx mindist|install_gorder|cargo.*install|rustc/ {print}'

echo '=== REP01 RECOVERY ==='
sha256sum "$traj/rep01/trajectory_provenance.pre_qc.json"
find "$root/audit/pbc_resume/20260822T093916Z/rep01" -maxdepth 1 -type f -printf '%f\t%s\n' | sort

echo '=== REP02 ==='
tail -n 20 "$root/audit/postproduction_runtime/rep02_primary_pbc.log" 2>/dev/null || true
find "$traj/rep02" -maxdepth 1 -type f -printf '%f\t%s\t%TY-%Tm-%TdT%TH:%TM:%TSZ\n' 2>/dev/null | sort

echo '=== TOOLCHAINS ==='
(cd "$gmx_audit" && sha256sum -c TOOLCHAIN_ACCEPTANCE.tsv.sha256) || true
tail -n 50 "$gorder_audit/gorder_install.log" 2>/dev/null || true

echo '=== RESOURCES ==='
free -h
df -h /root/autodl-tmp
cat /sys/fs/cgroup/memory.events 2>/dev/null || true
exit 0
