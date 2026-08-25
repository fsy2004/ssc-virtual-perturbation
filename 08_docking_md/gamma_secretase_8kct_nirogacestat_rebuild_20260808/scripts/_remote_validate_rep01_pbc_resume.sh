#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
out="$root/analysis/trajectories/8kct_nirogacestat_native/rep01"
recovery="$root/audit/pbc_resume/20260822T093916Z/rep01"
cd "$out"
sha256sum -c trajectory_provenance.pre_qc.json.sha256
cd "$recovery"
sha256sum -c RECOVERY_ARCHIVE.json.sha256

/root/miniconda3/bin/python - "$out" "$recovery" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
recovery = Path(sys.argv[2])
pbc = json.loads((out / "11_pbc_distance_invariance.json").read_text())
prov = json.loads((out / "trajectory_provenance.pre_qc.json").read_text())
archive = json.loads((recovery / "RECOVERY_ARCHIVE.json").read_text())
assert pbc["status"] == "pass", pbc
assert pbc["tolerance_nm"] == 0.01, pbc
assert pbc["maximum_absolute_difference_nm"] <= 0.01, pbc
assert prov["status"] == "pass_pending_scientific_qc_seal", prov
assert prov["realization_id"] == "rep01", prov
assert prov["production_tpr_sha256"] == "fd11c7287d5670c81ccb44fcb5b4215344726989f66f4f55db33643ba618678f"
assert archive["status"] == "archived_before_resume", archive
for name in (
    "check_centered_resume_preflight.json",
    "check_fitted_resume_preflight.json",
    "check_fixed_resume_preflight.json",
    "check_centered.json",
    "check_fitted.json",
    "check_fixed.json",
):
    check = json.loads((out / name).read_text())
    assert check["returncode"] == 0, (name, check["returncode"])
print(json.dumps({
    "status": "pass",
    "maximum_absolute_difference_nm": pbc["maximum_absolute_difference_nm"],
    "frame_count": pbc["frame_count"],
    "provenance_sha256": hashlib.sha256((out / "trajectory_provenance.pre_qc.json").read_bytes()).hexdigest(),
    "archived_partial_count": len(archive["archived_partial_outputs"]),
}, sort_keys=True))
PY
