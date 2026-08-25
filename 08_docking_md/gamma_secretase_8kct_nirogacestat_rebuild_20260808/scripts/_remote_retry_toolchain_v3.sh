#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
audit="$root/audit/toolchain_install/20260822_network_retry_v3/$stamp"
prior="$root/audit/toolchain_install/20260822_cpu_clone_resume/20260822T093916Z"
mkdir -p "$audit" /root/autodl-tmp/envs /root/autodl-tmp/tools

for pid in 1434 1435 1436; do
  [[ ! -r "/proc/$pid/stat" ]] || { echo "refusing active prior installer pid=$pid" >&2; exit 3; }
done

for log in "$prior/gmx_mmpbsa_install.log" "$prior/gorder_install.log" "$prior/fatslim_install.log"; do
  [[ -f "$log" ]] && sha256sum "$log" >> "$audit/PRIOR_FAILURE_LOGS.sha256"
done

for old_prefix in \
  /root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v2 \
  /root/autodl-tmp/tools/gorder_1_5_0_20260822_v2 \
  /root/autodl-tmp/tools/fatslim_0_2_2_20260822_v2
do
  if [[ -e "$old_prefix" ]]; then
    archived="${old_prefix}.failed_network_${stamp}"
    [[ ! -e "$archived" ]] || { echo "archive target exists: $archived" >&2; exit 4; }
    mv "$old_prefix" "$archived"
    printf '%s\t%s\n' "$old_prefix" "$archived" >> "$audit/ARCHIVED_FAILED_PREFIXES.tsv"
  fi
done

gmx_prefix=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v3
gorder_prefix=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v3
fatslim_prefix=/root/autodl-tmp/tools/fatslim_0_2_2_20260822_v3
for target in "$gmx_prefix" "$gorder_prefix" "$fatslim_prefix"; do
  [[ ! -e "$target" ]] || { echo "refusing existing v3 target: $target" >&2; exit 5; }
done

nohup env CONDA_EXE=/root/miniconda3/bin/conda \
  bash "$root/scripts/install_gmx_mmpbsa_1_6_5_cpu.sh" "$gmx_prefix" \
  > "$audit/gmx_mmpbsa_install.log" 2>&1 < /dev/null &
gmx_pid=$!

nohup env CARGO_BUILD_JOBS=30 \
  bash "$root/scripts/install_gorder_1_5_0.sh" "$gorder_prefix" \
  > "$audit/gorder_install.log" 2>&1 < /dev/null &
gorder_pid=$!

nohup bash "$root/scripts/install_fatslim_0_2_2.sh" "$fatslim_prefix" \
  > "$audit/fatslim_install.log" 2>&1 < /dev/null &
fatslim_pid=$!

cat > "$audit/RETRY_LAUNCH.tsv" <<EOF
work	pid	path	log
gmx_mmpbsa	$gmx_pid	$gmx_prefix	$audit/gmx_mmpbsa_install.log
gorder	$gorder_pid	$gorder_prefix	$audit/gorder_install.log
fatslim	$fatslim_pid	$fatslim_prefix	$audit/fatslim_install.log
EOF
sha256sum "$audit/RETRY_LAUNCH.tsv" > "$audit/RETRY_LAUNCH.tsv.sha256"
printf 'STAMP=%s\nGMX_MMPBSA_PID=%s\nGORDER_PID=%s\nFATSLIM_PID=%s\n' \
  "$stamp" "$gmx_pid" "$gorder_pid" "$fatslim_pid"
