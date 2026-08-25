#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
cd "$root"
bash -n \
  scripts/_remote_resume_cpu_clone_work.sh \
  scripts/install_gmx_mmpbsa_1_6_5_cpu.sh \
  scripts/install_gorder_1_5_0.sh \
  scripts/install_fatslim_0_2_2.sh
/root/miniconda3/bin/python -c 'import sys; sys.path.insert(0, "scripts"); import resume_primary_pbc_trajectories; print("REMOTE_IMPORT_OK")'
[[ -x /root/miniconda3/bin/conda ]]
[[ -x /root/GROMACS-2025.2/bin/gmx ]]
echo REMOTE_PREFLIGHT_OK
