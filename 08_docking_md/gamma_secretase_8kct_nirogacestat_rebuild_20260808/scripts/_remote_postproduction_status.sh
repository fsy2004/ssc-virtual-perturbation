#!/bin/bash
set -euo pipefail
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
/root/miniconda3/bin/python3 - <<'PY'
import json

manifest = json.load(open("config/study_manifest.json", encoding="utf-8"))
print("manifest_status=", manifest["manifest_status"])
print("temperature_k=", manifest["global_model"]["temperature_k"])
print("time_step_ps=", manifest["simulation"]["time_step_ps"])
print("hmr=", manifest["simulation"]["hydrogen_mass_repartitioning"])
PY
ps -o pid=,ppid=,etime=,stat=,cmd= -p 249195,249196 || true
