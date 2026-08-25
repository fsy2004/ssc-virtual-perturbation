#!/bin/bash
set -euo pipefail
root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
python=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
gmx=/root/GROMACS-2025.2/bin/gmx
test -d "$root"
test -x "$python"
test -x "$gmx"
printf '%s\n' '--- disk ---'
df -h /root/autodl-tmp
printf '%s\n' '--- active production or PBC processes ---'
pgrep -af 'mdrun|^/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python scripts/prepare_primary_pbc_trajectories.py |gmx check' || true
printf '%s\n' '--- output directories ---'
find "$root" -maxdepth 3 -type d -name 'analysis_primary*' -print
