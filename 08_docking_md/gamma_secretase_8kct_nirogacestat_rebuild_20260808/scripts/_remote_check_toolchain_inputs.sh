#!/bin/bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
for path in \
  /root/miniconda3/bin/conda \
  "$root/scripts/install_gmx_mmpbsa_1_6_5_cpu.sh" \
  "$root/scripts/capture_gmx_mmpbsa_toolchain.py"; do
  if [[ -f "$path" ]]; then
    sha256sum "$path"
  else
    printf 'MISSING %s\n' "$path"
  fi
done
/root/miniconda3/bin/conda --version
