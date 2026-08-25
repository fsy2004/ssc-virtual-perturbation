#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
prefix=/root/autodl-tmp/tools/gorder_1_5_0_20260822_v5
audit="$root/audit/toolchain_install/20260822_gorder_retry_v5/20260822T113333Z"
expected_commit=1beece37dc58a819be0a20b3ec691ef6cade365d

if [[ -r /proc/4064/cmdline ]] && tr '\0' ' ' < /proc/4064/cmdline | grep -q 'gorder'; then
  echo 'gorder installer is still active' >&2
  exit 3
fi
[[ -x "$prefix/install/bin/gorder" ]]
[[ "$(tr -d '\r\n' < "$prefix/GORDER_SOURCE_COMMIT.txt")" == "$expected_commit" ]]
(cd "$prefix" && sha256sum -c "$(basename "$prefix/gorder-v1.5.0-source.tar.gz.sha256")")
(cd "$prefix" && sha256sum -c gorder-binary.sha256)
[[ "$("$prefix/install/bin/gorder" --version)" == 'gorder 1.5.0' ]]
grep -q '^release: 1.87.0$' "$prefix/rustc-version.txt"
grep -q '^cargo 1.87.0 ' "$prefix/cargo-version.txt"

sha256sum \
  "$prefix/GORDER_SOURCE_COMMIT.txt" \
  "$prefix/gorder-v1.5.0-source.tar.gz" \
  "$prefix/gorder-binary.sha256" \
  "$prefix/gorder-version.txt" \
  "$prefix/rustc-version.txt" \
  "$prefix/cargo-version.txt" \
  "$prefix/source/Cargo.toml" \
  "$prefix/source/Cargo.lock" > "$audit/ACCEPTANCE_ARTIFACTS.sha256"
printf 'status\tpass\nprefix\t%s\nsource_commit\t%s\ngorder_version\t1.5.0\nrust_version\t1.87.0\n' \
  "$prefix" "$expected_commit" > "$audit/TOOLCHAIN_ACCEPTANCE.tsv"
sha256sum "$audit/TOOLCHAIN_ACCEPTANCE.tsv" > "$audit/TOOLCHAIN_ACCEPTANCE.tsv.sha256"
cat "$audit/TOOLCHAIN_ACCEPTANCE.tsv"
