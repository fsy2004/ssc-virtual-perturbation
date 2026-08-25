#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
tool_base="$root/audit/toolchain_install/20260822_cpu_clone_resume"
latest=$(find "$tool_base" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort | tail -1)
tool_audit="$tool_base/$latest"
pbc_audit="$root/audit/pbc_resume/$latest"
echo "STAMP=$latest"
cat "$tool_audit/RESUME_LAUNCH.tsv"

echo PROCESSES
tail -n +2 "$tool_audit/RESUME_LAUNCH.tsv" | while IFS=$'\t' read -r work pid path log; do
  if [[ -r "/proc/$pid/stat" ]]; then
    state=$(awk '{print $3}' "/proc/$pid/stat")
    rss_kib=$(awk '/VmRSS:/ {print $2}' "/proc/$pid/status")
    printf '%s\t%s\trunning\tstate=%s\trss_kib=%s\n' "$work" "$pid" "$state" "${rss_kib:-0}"
  else
    printf '%s\t%s\texited\n' "$work" "$pid"
  fi
done

echo MATCHING_CHILDREN
pgrep -af 'resume_primary_pbc_trajectories|gmx[[:space:]]+(mindist|check)|conda.*create|cargo.*install|rustup-init|install_gmx_mmpbsa|install_gorder|install_fatslim' || true

echo TOP_CPU
ps -eo pid=,ppid=,stat=,pcpu=,pmem=,rss=,etime=,args= --sort=-pcpu | head -18

echo LOG_TAILS
for log in \
  "$tool_audit/gmx_mmpbsa_install.log" \
  "$tool_audit/gorder_install.log" \
  "$tool_audit/fatslim_install.log" \
  "$pbc_audit/rep01_runner.log"
do
  echo "FILE=$log"
  tail -12 "$log" 2>/dev/null || true
done

echo PBC_FILES
find "$root/analysis/trajectories/8kct_nirogacestat_native/rep01" -maxdepth 1 -type f \
  \( -name '09_*' -o -name '10_*' -o -name '11_*' -o -name 'trajectory_provenance*' \) \
  -printf '%f\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS%Tz\n' | sort
find "$pbc_audit/rep01" -maxdepth 1 -type f -printf 'ARCHIVE\t%f\t%s\n' 2>/dev/null | sort || true

echo RESOURCES
uptime
free -h
df -h /root/autodl-tmp
cat /sys/fs/cgroup/memory.events 2>/dev/null || true
