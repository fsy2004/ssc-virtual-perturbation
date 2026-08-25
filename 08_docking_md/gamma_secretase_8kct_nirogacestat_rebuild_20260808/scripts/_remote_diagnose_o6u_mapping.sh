set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
cd "$ROOT"
PYTHONPATH=scripts "$PY" - <<'PY'
import json
import numpy as np
import MDAnalysis as mda
from pathlib import Path

root=Path('/root/autodl-tmp/o6u_md_release_3x500ns_v4')
mapping=json.loads((root/'config/primary_atom_mapping_contacts.json').read_text())
u=mda.Universe(str(root/'rep01/work/production.tpr'))
o6u=np.flatnonzero(np.asarray(u.atoms.resnames)=='O6U')
heavy=set(int(i) for i in o6u[np.asarray(u.atoms[o6u].masses)>2.0])
mapped=set(int(row['trajectory']['index']) for row in mapping['atom_mappings']['o6u_heavy'])
print('trajectory_o6u_count',len(o6u),'heavy',len(heavy),'mapped',len(mapped))
for label, indices in [('missing',sorted(heavy-mapped)),('extra',sorted(mapped-heavy))]:
    print(label,len(indices))
    for i in indices:
        a=u.atoms[i]
        print(i,a.name,a.resname,a.resid,a.segid,a.mass)
print('mapped endpoints')
for row in mapping['atom_mappings']['o6u_heavy']:
    t=row['trajectory']; r=row['reference']
    print(r['name'],'=>',t['index'],t['name'],t['resname'],u.atoms[t['index']].mass)
PY
