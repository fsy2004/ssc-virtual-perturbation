#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
tool_audit="$root/audit/toolchain_install/20260822_cpu_clone_resume/$stamp"
pbc_audit="$root/audit/pbc_resume/$stamp"
mkdir -p "$tool_audit" "$pbc_audit" /root/autodl-tmp/envs /root/autodl-tmp/tools

{
  pgrep -af 'prepare_primary_pbc_trajectories|resume_primary_pbc_trajectories|gmx[[:space:]]+mindist|install_gmx_mmpbsa|install_gorder|install_fatslim' || true
} | awk -v self="$$" -v parent="$PPID" '
  $1 != self && $1 != parent && $0 !~ /pgrep -af/
' > "$tool_audit/process_guard.txt"
if [[ -s "$tool_audit/process_guard.txt" ]]; then
  echo "refusing resume because matching work is active" >&2
  cat "$tool_audit/process_guard.txt" >&2
  exit 3
fi

old_audit="$root/audit/toolchain_install/20260822"
mkdir -p "$tool_audit/prior_logs"
for old_file in \
  "$old_audit/gmx_mmpbsa_install.log" "$old_audit/gmx_mmpbsa_install.pid" \
  "$old_audit/gorder_install.log" "$old_audit/gorder_install.pid" \
  "$old_audit/fatslim_install.log" "$old_audit/fatslim_install.pid"
do
  if [[ -e "$old_file" ]]; then
    mv "$old_file" "$tool_audit/prior_logs/"
  fi
done

for old_prefix in \
  /root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822 \
  /root/autodl-tmp/tools/gorder_1_5_0_20260822 \
  /root/autodl-tmp/tools/fatslim_0_2_2_20260822
do
  if [[ -e "$old_prefix" ]]; then
    archived="${old_prefix}.paused_${stamp}"
    [[ ! -e "$archived" ]] || { echo "archive target exists: $archived" >&2; exit 4; }
    mv "$old_prefix" "$archived"
    printf '%s\t%s\n' "$old_prefix" "$archived" >> "$tool_audit/archived_partial_prefixes.tsv"
  fi
done

gmx_prefix=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v2
gorder_prefix=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v2
fatslim_prefix=/root/autodl-tmp/tools/fatslim_0_2_2_20260822_v2
for target in "$gmx_prefix" "$gorder_prefix" "$fatslim_prefix"; do
  [[ ! -e "$target" ]] || { echo "refusing existing fresh target: $target" >&2; exit 5; }
done

nohup env CONDA_EXE=/root/miniconda3/bin/conda \
  bash "$root/scripts/install_gmx_mmpbsa_1_6_5_cpu.sh" "$gmx_prefix" \
  > "$tool_audit/gmx_mmpbsa_install.log" 2>&1 < /dev/null &
gmx_pid=$!

nohup env CARGO_BUILD_JOBS=30 \
  bash "$root/scripts/install_gorder_1_5_0.sh" "$gorder_prefix" \
  > "$tool_audit/gorder_install.log" 2>&1 < /dev/null &
gorder_pid=$!

nohup bash "$root/scripts/install_fatslim_0_2_2.sh" "$fatslim_prefix" \
  > "$tool_audit/fatslim_install.log" 2>&1 < /dev/null &
fatslim_pid=$!

nohup /root/miniconda3/bin/python "$root/scripts/resume_primary_pbc_trajectories.py" \
  --release-root "$root" --replica rep01 --recovery-id "$stamp" \
  > "$pbc_audit/rep01_runner.log" 2>&1 < /dev/null &
pbc_pid=$!

cat > "$tool_audit/RESUME_LAUNCH.tsv" <<EOF
work	pid	path	log
gmx_mmpbsa	$gmx_pid	$gmx_prefix	$tool_audit/gmx_mmpbsa_install.log
gorder	$gorder_pid	$gorder_prefix	$tool_audit/gorder_install.log
fatslim	$fatslim_pid	$fatslim_prefix	$tool_audit/fatslim_install.log
rep01_pbc_resume	$pbc_pid	$pbc_audit	$pbc_audit/rep01_runner.log
EOF
sha256sum "$tool_audit/RESUME_LAUNCH.tsv" > "$tool_audit/RESUME_LAUNCH.tsv.sha256"

printf 'STAMP=%s\nGMX_MMPBSA_PID=%s\nGORDER_PID=%s\nFATSLIM_PID=%s\nREP01_PBC_PID=%s\n' \
  "$stamp" "$gmx_pid" "$gorder_pid" "$fatslim_pid" "$pbc_pid"
