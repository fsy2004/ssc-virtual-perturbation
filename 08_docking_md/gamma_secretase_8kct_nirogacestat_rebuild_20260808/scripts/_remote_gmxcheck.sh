#!/bin/bash
set -e
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
GMX=/root/GROMACS-2025.2/bin/gmx
for rep in rep01 rep02 rep03; do
  echo "==== $rep ===="
  $GMX check -f $rep/work/production.xtc -s1 $rep/work/production.tpr 2>&1 | tail -6
  echo "  xtc-exit=$?"
  $GMX check -e $rep/work/production.edr 2>&1 | tail -4
  echo "  edr-exit=$?"
  $GMX check -f $rep/work/production.gro 2>&1 | tail -3
  echo "  gro-exit=$?"
done
