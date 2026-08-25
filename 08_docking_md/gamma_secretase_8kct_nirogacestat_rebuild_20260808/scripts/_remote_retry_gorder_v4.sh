#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
audit="$root/audit/toolchain_install/20260822_gorder_retry_v4/$stamp"
prior_audit="$root/audit/toolchain_install/20260822_network_retry_v3/20260822T110831Z"
old=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v3
fresh=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v4
mkdir -p "$audit"

if [[ -r /proc/2006/cmdline ]] && tr '\0' ' ' < /proc/2006/cmdline | grep -q 'gorder'; then
  echo 'refusing active prior gorder installer pid=2006' >&2
  exit 3
fi
[[ ! -e "$fresh" ]] || { echo "refusing existing v4 target: $fresh" >&2; exit 4; }

if [[ -f "$prior_audit/gorder_install.log" ]]; then
  sha256sum "$prior_audit/gorder_install.log" > "$audit/PRIOR_CARGO182_FAILURE.sha256"
fi
if [[ -e "$old" ]]; then
  archived="${old}.failed_cargo182_${stamp}"
  [[ ! -e "$archived" ]]
  mv "$old" "$archived"
  printf '%s\t%s\n' "$old" "$archived" > "$audit/ARCHIVED_FAILED_PREFIX.tsv"
fi

nohup env CARGO_BUILD_JOBS=30 \
  bash "$root/scripts/install_gorder_1_5_0.sh" "$fresh" \
  > "$audit/gorder_install.log" 2>&1 < /dev/null &
pid=$!
printf 'work\tpid\tpath\tlog\ngorder\t%s\t%s\t%s\n' \
  "$pid" "$fresh" "$audit/gorder_install.log" > "$audit/RETRY_LAUNCH.tsv"
sha256sum "$audit/RETRY_LAUNCH.tsv" > "$audit/RETRY_LAUNCH.tsv.sha256"
printf 'STAMP=%s\nGORDER_PID=%s\n' "$stamp" "$pid"
