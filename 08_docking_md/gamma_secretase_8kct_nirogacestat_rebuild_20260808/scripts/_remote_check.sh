#!/bin/bash
set -e
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
echo ====REP02_COMPLETION====
/root/miniconda3/bin/python - <<'PYEOF'
import json
d = json.load(open('rep02/PRODUCTION_COMPLETION_500NS.json'))
print('status', d['status'])
print('final_step', d['final_step'])
print('final_time_ps', d['final_time_ps'])
print('tpr_sha', d['production_tpr_sha256'])
print('recovery_bound', d['checks']['cuda_recovery_audit_bound'])
PYEOF
echo ====REP03_CPT====
ls -la rep03/work/production_prev.cpt rep03/work/production.cpt
echo ====REP03_TPR_SHA====
sha256sum rep03/work/production.tpr
echo ====REP03_LOG_LAST====
tail -5 rep03/work/production.log
