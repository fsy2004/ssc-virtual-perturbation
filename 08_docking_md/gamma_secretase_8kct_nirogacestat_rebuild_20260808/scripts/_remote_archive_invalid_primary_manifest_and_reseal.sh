#!/usr/bin/env bash
set -euo pipefail
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
stamp=20260822T1430_old_draft_placeholders
faildir="audit/primary_manifest_seal_failed/${stamp}"
mkdir -p "$faildir"
if [[ -e config/primary_postprocessing_manifest.approved.json ]]; then
  mv config/primary_postprocessing_manifest.approved.json "$faildir/primary_postprocessing_manifest.approved.json"
fi
if [[ -e config/primary_postprocessing_manifest.approved.json.sha256 ]]; then
  mv config/primary_postprocessing_manifest.approved.json.sha256 "$faildir/primary_postprocessing_manifest.approved.json.sha256"
fi
find "$faildir" -maxdepth 1 -type f -printf 'ARCHIVED\t%f\t%s\n' | sort
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python scripts/seal_primary_postprocessing_manifest.py \
  --draft config/primary_postprocessing_manifest.json \
  --release-root /root/autodl-tmp/o6u_md_release_3x500ns_v4 \
  --protocol config/production_protocol_hmr4fs_303K_v1.json \
  --output config/primary_postprocessing_manifest.approved.json
