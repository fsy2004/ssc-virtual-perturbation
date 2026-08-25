#!/bin/bash
set -euo pipefail
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
PYTHONPATH=scripts /root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python - <<'PY'
import MDAnalysis as mda
from analyze_primary_structure_mdanalysis import topology_identity_sha256

universe = mda.Universe('rep01/work/production.tpr')
print(len(universe.atoms), topology_identity_sha256(universe))
for index in (0, 93222, 123652, 330582):
    atom = universe.atoms[index]
    print(index, atom.name, atom.resname, atom.resid, repr(atom.segid), repr(getattr(atom, 'chainID', '')))
PY
