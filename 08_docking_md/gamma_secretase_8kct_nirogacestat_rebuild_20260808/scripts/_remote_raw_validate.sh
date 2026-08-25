#!/bin/bash
set -e
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
GMX=/root/GROMACS-2025.2/bin/gmx
echo "=== validate_md_outputs production ==="
$PY scripts/validate_md_outputs.py --manifest config/study_manifest.json --phase production --strict --gmx $GMX --report reports/production_output_validation.json 2>&1 | tail -50
echo "EXIT=$?"
echo "=== report? ==="
ls -la reports/production_output_validation.json 2>/dev/null && echo "REPORT OK" || echo "NO REPORT"
