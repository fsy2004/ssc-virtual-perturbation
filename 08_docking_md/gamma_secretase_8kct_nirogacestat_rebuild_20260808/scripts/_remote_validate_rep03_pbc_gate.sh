#!/usr/bin/env bash
set -euo pipefail
export REPLICA=rep03
cd /root/autodl-tmp/o6u_md_release_3x500ns_v4
/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python - <<'PY'
import hashlib
import json
import os
import pathlib
import sys

replica = os.environ.get("REPLICA", "rep03")
if replica not in {"rep01", "rep02", "rep03"}:
    raise SystemExit(f"invalid replica: {replica}")

root = pathlib.Path("/root/autodl-tmp/o6u_md_release_3x500ns_v4")
out = root / "analysis/trajectories/8kct_nirogacestat_native" / replica
errors: list[str] = []

pbc = json.loads((out / "11_pbc_distance_invariance.json").read_text())
if pbc.get("status") != "pass":
    errors.append(f"pbc status {pbc.get('status')}")
if pbc.get("tolerance_nm") != 0.01:
    errors.append(f"pbc tolerance {pbc.get('tolerance_nm')}")
if pbc.get("maximum_absolute_difference_nm", 999.0) > 0.01:
    errors.append(f"pbc maxdiff {pbc.get('maximum_absolute_difference_nm')}")
if int(pbc.get("frame_count", -1)) != 25001:
    errors.append(f"pbc frame_count {pbc.get('frame_count')}")

checks: dict[str, str | None] = {}
for name in ("check_centered.json", "check_fitted.json", "check_fixed.json"):
    data = json.loads((out / name).read_text())
    checks[name] = str(data.get("returncode"))
    if data.get("returncode") != 0:
        errors.append(f"{name} returncode {data.get('returncode')}")

prov = out / "trajectory_provenance.pre_qc.json"
sidecar = out / "trajectory_provenance.pre_qc.json.sha256"
provenance_sha = hashlib.sha256(prov.read_bytes()).hexdigest()
sidecar_sha = sidecar.read_text().split()[0]
if provenance_sha != sidecar_sha:
    errors.append("provenance sha sidecar mismatch")
prov_data = json.loads(prov.read_text())
if prov_data.get("status") != "pass_pending_scientific_qc_seal":
    errors.append(f"provenance status {prov_data.get('status')}")

report = {
    "replica": replica,
    "status": "pass" if not errors else "fail",
    "errors": errors,
    "pbc": {
        "maximum_absolute_difference_nm": pbc.get("maximum_absolute_difference_nm"),
        "frame_count": pbc.get("frame_count"),
        "tolerance_nm": pbc.get("tolerance_nm"),
    },
    "checks": checks,
    "provenance_sha256": provenance_sha,
}
print(json.dumps(report, indent=2, sort_keys=True))
if errors:
    sys.exit(2)
PY
