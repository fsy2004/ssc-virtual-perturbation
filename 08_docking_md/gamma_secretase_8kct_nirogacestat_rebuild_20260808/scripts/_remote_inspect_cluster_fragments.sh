set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
cd "$ROOT"
"$PY" - <<'PY'
import MDAnalysis as mda
from pathlib import Path

root=Path('/root/autodl-tmp/o6u_md_release_3x500ns_v4')
u=mda.Universe(str(root/'rep01/work/production.tpr'))
seed=u.select_atoms('protein or resname O6U')
fragments=seed.fragments
indices=sorted({int(i) for fragment in fragments for i in fragment.indices})
print('atoms',len(u.atoms),'seed',len(seed),'fragments',len(fragments),'complete_group',len(indices))
for n,fragment in enumerate(fragments):
    resnames=sorted(set(str(x) for x in fragment.resnames))
    segids=sorted(set(str(x) for x in fragment.segids))
    print(n,'atoms',len(fragment),'range',int(fragment.indices.min()),int(fragment.indices.max()),'resnames',','.join(resnames),'segids',','.join(segids))
print('added_to_seed',len(set(indices)-set(int(i) for i in seed.indices)))
for i in sorted(set(indices)-set(int(i) for i in seed.indices))[:30]:
    a=u.atoms[i]
    print('added',i,a.name,a.resname,a.resid,a.segid)
PY
