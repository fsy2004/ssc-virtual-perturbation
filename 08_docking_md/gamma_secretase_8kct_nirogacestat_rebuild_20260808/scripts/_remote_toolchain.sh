#!/bin/bash
set -e
cd /root/autodl-tmp/ssc_md_work/gamma_secretase_8kct_nirogacestat_rebuild_20260808/server_records/toolchain/analysis_environment_v3
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
echo ====RECORD====
$PY - <<'PY'
import json
d = json.load(open('SERVER_TOOLCHAIN_RECORD.json'))
print('keys:', list(d.keys()))
for k in ('gromacs_version','gromacs','gmx','executable','toolchain'):
    if k in d:
        print(k, '=', json.dumps(d[k], default=str)[:300])
PY
echo ====GMX_IDENTITY====
$PY - <<'PY'
import json
d = json.load(open('gromacs_executable_identity.json'))
print('keys:', list(d.keys()))
print(json.dumps(d, default=str)[:400])
PY
