#!/bin/bash
set -e
echo ====MD_WORK====
ls /root/autodl-tmp/ssc_md_work/ 2>/dev/null | head -20
echo ====MD_WORK_PROJ====
ls /root/autodl-tmp/ssc_md_work/gamma_secretase_8kct_nirogacestat_rebuild_20260808/ 2>/dev/null | head -40
echo ====ANALYSIS_ENV====
ls /root/autodl-tmp/envs/ 2>/dev/null
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python -c 'import MDAnalysis; print("MDAnalysis", MDAnalysis.__version__)' 2>&1 | head -2
echo ====ENV_PKGS====
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python -m pip list 2>/dev/null | grep -iE 'mdanalysis|numpy|scipy|pandas|matplotlib|parmed' | head
