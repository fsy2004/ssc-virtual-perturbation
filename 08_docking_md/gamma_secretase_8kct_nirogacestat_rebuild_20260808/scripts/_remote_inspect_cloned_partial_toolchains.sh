#!/usr/bin/env bash
set -euo pipefail

release=/root/autodl-tmp/o6u_md_release_3x500ns_v4
printf 'MINICONDA='
if [[ -x /root/miniconda3/bin/conda ]]; then
  echo yes
else
  echo no
fi

echo PARTIALS
for prefix in \
  /root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822 \
  /root/autodl-tmp/tools/gorder_1_5_0_20260822 \
  /root/autodl-tmp/tools/fatslim_0_2_2_20260822
do
  if [[ -e "$prefix" ]]; then
    du -sh "$prefix"
    find "$prefix" -maxdepth 2 -type f -printf '%p\t%s\n' | sort | tail -20
  else
    echo "MISSING:$prefix"
  fi
done

echo LOGTAILS
shopt -s nullglob
for log in "$release"/audit/toolchain_install/20260822/*.log; do
  echo "FILE:$log"
  tail -8 "$log"
done
