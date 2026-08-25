#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=$(date -u +%Y%m%dT%H%M%SZ)
audit="$root/audit/toolchain_install/20260822_gorder_retry_v5/$stamp"
prior_audit="$root/audit/toolchain_install/20260822_gorder_retry_v4/20260822T112115Z"
old=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v4
fresh=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v5
mkdir -p "$audit"

if [[ -r /proc/2986/cmdline ]] && tr '\0' ' ' < /proc/2986/cmdline | grep -q 'gorder'; then
  echo 'refusing active prior gorder installer pid=2986' >&2
  exit 3
fi
[[ ! -e "$fresh" ]] || { echo "refusing existing v5 target: $fresh" >&2; exit 4; }

if [[ -f "$prior_audit/gorder_install.log" ]]; then
  sha256sum "$prior_audit/gorder_install.log" > "$audit/PRIOR_RUST185_FAILURE.sha256"
fi
if [[ -e "$old" ]]; then
  archived="${old}.failed_rust185_${stamp}"
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
