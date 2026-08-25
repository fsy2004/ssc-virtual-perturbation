#!/bin/bash
set -euo pipefail

printf '%s\n' '--- build tools ---'
for tool in gcc g++ make cmake pkg-config; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '%s=%s\n' "$tool" "$(command -v "$tool")"
    "$tool" --version 2>/dev/null | head -n 1 || true
  else
    printf '%s=MISSING\n' "$tool"
  fi
done

printf '%s\n' '--- gorder crates.io ---'
curl -A 'o6u-md-reproducibility-audit/1.0' -fsSL https://crates.io/api/v1/crates/gorder \
  | /root/miniconda3/bin/python -c 'import json,sys; d=json.load(sys.stdin); c=d["crate"]; print("newest_version="+c["newest_version"]); print("max_version="+c["max_version"]); print("repository="+str(c.get("repository"))); print("updated_at="+str(c.get("updated_at")))'

printf '%s\n' '--- git refs ---'
git ls-remote https://github.com/VachaLab/gorder.git HEAD 'refs/tags/*' | tail -n 20
git ls-remote https://github.com/FATSLiM/fatslim.git HEAD 'refs/tags/*' | tail -n 20

printf '%s\n' '--- PyPI FATSLiM ---'
curl -A 'o6u-md-reproducibility-audit/1.0' -fsSL https://pypi.org/pypi/fatslim/json \
  | /root/miniconda3/bin/python -c 'import json,sys; d=json.load(sys.stdin); print("version="+d["info"]["version"]); print("project_url="+str(d["info"].get("project_url"))); print("release_versions="+",".join(sorted(d["releases"])))'
