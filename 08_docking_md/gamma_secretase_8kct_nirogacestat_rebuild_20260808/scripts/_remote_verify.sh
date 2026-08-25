#!/bin/bash
set -e
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
echo ====PY_SYNTAX====
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python -m py_compile scripts/make_analysis_trajectories.py scripts/validate_md_outputs.py scripts/validate_qc_stationarity_report.py scripts/md_contract.py scripts/primary_postprocessing_common.py scripts/analyze_primary_structure_mdanalysis.py scripts/analyze_membrane_qc_mdanalysis.py scripts/gmx_energy_qc.py scripts/validate_primary_postprocessing.py 2>&1 | head -5
echo "PY_COMPILE_EXIT=$?"
echo ====NDX_CHECK====
/root/GROMACS-2025.2/bin/gmx select -s common/minimized.gro -f common/minimized.gro -on /tmp/ndx_test.ndx -select 'none' 2>&1 | tail -3
echo ====MANIFEST_PARSE====
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python -c "
import json
m = json.load(open('config/study_manifest.json'))
print('manifest_status:', m['manifest_status'])
print('run dirs:', [r['run_directory'] for r in m['systems'][0]['realizations']])
p = json.load(open('config/analysis_plan.json'))
print('plan groups:', p['trajectory_processing']['groups']['fit'])
print('window:', p['trajectory_processing']['fixed_window_ns'])
"
echo ====NDX_GROUPS====
grep -E '^\[' builds/analysis.ndx
