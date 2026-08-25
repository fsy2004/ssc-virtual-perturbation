#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
audit="$root/audit/toolchain_install/20260822_gmx_mmpbsa_retry_v5/$stamp"
prior_audit="$root/audit/toolchain_install/20260822_gmx_mmpbsa_retry_v4/20260822T111420Z"
old=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v4
fresh=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v5
mkdir -p "$audit"

if [[ -r /proc/2441/cmdline ]] && tr '\0' ' ' < /proc/2441/cmdline | grep -q 'gmx_mmpbsa'; then
  echo 'refusing active prior gmx_MMPBSA installer pid=2441' >&2
  exit 3
fi
[[ ! -e "$fresh" ]] || { echo "refusing existing v5 target: $fresh" >&2; exit 4; }

if [[ -f "$prior_audit/gmx_mmpbsa_install.log" ]]; then
  sha256sum "$prior_audit/gmx_mmpbsa_install.log" > "$audit/PRIOR_MPI_SOLVER_FAILURE.sha256"
fi
if [[ -e "$old" ]]; then
  archived="${old}.failed_mpi_solver_${stamp}"
  [[ ! -e "$archived" ]]
  mv "$old" "$archived"
  printf '%s\t%s\n' "$old" "$archived" > "$audit/ARCHIVED_FAILED_PREFIX.tsv"
fi

nohup env CONDA_EXE=/root/miniconda3/bin/conda \
  bash "$root/scripts/install_gmx_mmpbsa_1_6_5_cpu.sh" "$fresh" \
  > "$audit/gmx_mmpbsa_install.log" 2>&1 < /dev/null &
pid=$!
printf 'work\tpid\tpath\tlog\ngmx_MMPBSA\t%s\t%s\t%s\n' \
  "$pid" "$fresh" "$audit/gmx_mmpbsa_install.log" > "$audit/RETRY_LAUNCH.tsv"
sha256sum "$audit/RETRY_LAUNCH.tsv" > "$audit/RETRY_LAUNCH.tsv.sha256"
printf 'STAMP=%s\nGMX_MMPBSA_PID=%s\n' "$stamp" "$pid"
