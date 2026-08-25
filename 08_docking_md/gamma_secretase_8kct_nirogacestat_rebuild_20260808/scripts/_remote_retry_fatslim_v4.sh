#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
audit="$root/audit/toolchain_install/20260822_fatslim_retry_v4/$stamp"
prior_audit="$root/audit/toolchain_install/20260822_network_retry_v3/20260822T110831Z"
old=/root/autodl-tmp/tools/fatslim_0_2_2_20260822_v3
fresh=/root/autodl-tmp/tools/fatslim_0_2_2_20260822_v4
mkdir -p "$audit"
[[ ! -r /proc/2007/stat ]] || { echo 'refusing active prior FATSLiM installer' >&2; exit 3; }
[[ ! -e "$fresh" ]] || { echo "refusing existing v4 target: $fresh" >&2; exit 4; }
if [[ -f "$prior_audit/fatslim_install.log" ]]; then
  sha256sum "$prior_audit/fatslim_install.log" > "$audit/PRIOR_SOLVER_FAILURE.sha256"
fi
if [[ -e "$old" ]]; then
  archived="${old}.failed_solver_${stamp}"
  [[ ! -e "$archived" ]]
  mv "$old" "$archived"
  printf '%s\t%s\n' "$old" "$archived" > "$audit/ARCHIVED_FAILED_PREFIX.tsv"
fi
nohup bash "$root/scripts/install_fatslim_0_2_2.sh" "$fresh" \
  > "$audit/fatslim_install.log" 2>&1 < /dev/null &
pid=$!
printf 'work\tpid\tpath\tlog\nFATSLiM\t%s\t%s\t%s\n' \
  "$pid" "$fresh" "$audit/fatslim_install.log" > "$audit/RETRY_LAUNCH.tsv"
sha256sum "$audit/RETRY_LAUNCH.tsv" > "$audit/RETRY_LAUNCH.tsv.sha256"
printf 'STAMP=%s\nFATSLIM_PID=%s\n' "$stamp" "$pid"
