#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
run_id="${O6U_PRIMARY_QC_RUN_ID:-20260822T144701Z_primary_qc_v1}"
out="$root/analysis/primary_postprocessing/$run_id"
runtime="$root/audit/postproduction_runtime"

echo "=== PRIMARY_QC_PROCESS ==="
ps -eo pid,ppid,stat,etime,pcpu,pmem,rss,args | awk \
  'NR==1 || /run_20260822T144701Z_primary_qc_v1|recover_20260822T144701Z_primary_qc_v1_after_legend_fix|analyze_primary_structure_mdanalysis|gmx_energy_qc.py|analyze_membrane_qc_mdanalysis|validate_primary_postprocessing.py/ {print}'

echo "=== STATUS_TSV ==="
cat "$runtime/run_${run_id}.status.tsv" 2>/dev/null || true

echo "=== RECOVERY_STATUS_TSV ==="
cat "$runtime/recover_${run_id}_after_legend_fix.status.tsv" 2>/dev/null || true

echo "=== LOG_TAILS ==="
for f in \
  "$out/logs/structural.log" \
  "$out/logs/energy.log" \
  "$out/logs/membrane.log" \
  "$out/logs/validate_primary.log" \
  "$runtime/run_${run_id}.log" \
  "$runtime/recover_${run_id}_after_legend_fix.log"
do
  echo "--- $f"
  tail -n 30 "$f" 2>/dev/null || true
done

echo "=== OUTPUT_TREE_HEAD ==="
find "$out" -maxdepth 3 -type f -printf '%P\t%s\t%TY-%Tm-%TdT%TH:%TM:%TSZ\n' 2>/dev/null | sort | head -n 80

echo "=== RESOURCES ==="
free -h
df -h /root/autodl-tmp
cat /sys/fs/cgroup/memory.events 2>/dev/null || true
