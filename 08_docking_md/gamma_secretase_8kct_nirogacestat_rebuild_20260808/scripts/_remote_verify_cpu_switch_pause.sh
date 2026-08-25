#!/bin/bash
set -euo pipefail

pause_dir=/root/autodl-tmp/o6u_md_release_3x500ns_v4/audit/cpu_instance_switch_pause/20260822T163915+0800
cd "$pause_dir"
cat PAUSE_METADATA.txt
printf '%s\n' '--- forced kill pids ---'
cat forced_kill_pids.txt
printf '%s\n' '--- evidence hashes ---'
sha256sum -c PAUSE_EVIDENCE_SHA256.txt
printf '%s\n' '--- matching live processes ---'
ps -eo pid,args \
  | grep -E 'prepare_primary_pbc_trajectories|install_gmx_mmpbsa_1_6_5_cpu|install_gorder_1_5_0|install_fatslim_0_2_2|gmx mindist|conda create|rustup-init' \
  | grep -v grep || true
