set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
TPR="$ROOT/rep01/work/production.tpr"
REFERENCE="$ROOT/docking_native_redock/plip_native/8KCT_protonated.pdb"
STAMP=$(date +%Y%m%dT%H%M%S%z)
ARCHIVE="$ROOT/audit/mapping_revisions/${STAMP}_pre_formal_tpr_mapping"

test -d "$ROOT"
test -x "$PY"
test -f "$TPR"
test -f "$REFERENCE"
mkdir -p "$ROOT/analysis_config_work" "$ARCHIVE"
if test ! -e "$ROOT/analysis_config_work/step5_input.pdb"; then
  ln -s ../common/step5_input.pdb "$ROOT/analysis_config_work/step5_input.pdb"
fi
if test ! -e "$ROOT/analysis_config_work/minimized.gro"; then
  ln -s ../common/minimized.gro "$ROOT/analysis_config_work/minimized.gro"
fi
for name in primary_atom_mapping_contacts.json primary_atom_mapping_contacts.json.sha256 membrane_qc_mapping.json membrane_qc_mapping.json.sha256; do
  if test -f "$ROOT/config/$name"; then
    cp -a "$ROOT/config/$name" "$ARCHIVE/$name"
  fi
done

cd "$ROOT"
PYTHONPATH=scripts "$PY" scripts/build_primary_mapping_records.py --trajectory-topology "$TPR"
PYTHONPATH=scripts "$PY" scripts/build_membrane_mapping.py --trajectory-topology "$TPR"

PYTHONPATH=scripts "$PY" - "$ROOT" "$TPR" "$REFERENCE" "$ARCHIVE" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import MDAnalysis as mda

import analyze_membrane_qc_mdanalysis as membrane
import analyze_primary_structure_mdanalysis as structure

root, tpr_path, reference_path, archive = map(Path, sys.argv[1:])

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def approve_and_validate(path, validator, *universes):
    draft = json.loads(path.read_text(encoding="utf-8"))
    draft_hash = sha256(path)
    candidate = dict(draft)
    candidate["approval_status"] = "approved"
    validator(candidate, *universes)
    path.write_text(json.dumps(candidate, indent=2) + "\n", encoding="utf-8")
    approved_hash = sha256(path)
    path.with_suffix(".json.sha256").write_text(f"{approved_hash}  {path.name}\n", encoding="ascii")
    validator(json.loads(path.read_text(encoding="utf-8")), *universes)
    return {"draft_sha256": draft_hash, "approved_sha256": approved_hash, "status": "pass"}

reference = mda.Universe(str(reference_path))
trajectory = mda.Universe(str(tpr_path))
structural_path = root / "config" / "primary_atom_mapping_contacts.json"
membrane_path = root / "config" / "membrane_qc_mapping.json"

structural = approve_and_validate(structural_path, structure._validate_mapping_record, reference, trajectory)
membrane_result = approve_and_validate(membrane_path, membrane._validate_mapping, trajectory)
report = {
    "schema_version": "1.0",
    "report_type": "formal_mapping_validation",
    "status": "pass",
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
    "trajectory_topology": {"path": str(tpr_path), "sha256": sha256(tpr_path), "atom_count": len(trajectory.atoms)},
    "reference": {"path": str(reference_path), "sha256": sha256(reference_path), "atom_count": len(reference.atoms)},
    "structural_mapping": structural,
    "membrane_mapping": membrane_result,
    "validation_functions": [
        "analyze_primary_structure_mdanalysis._validate_mapping_record",
        "analyze_membrane_qc_mdanalysis._validate_mapping",
    ],
    "interpretation": "Approval covers identity/mapping schema validation only; trajectory QC and scientific interpretation remain pending.",
}
report_path = archive / "FORMAL_MAPPING_VALIDATION.json"
report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
report_hash = sha256(report_path)
report_path.with_suffix(".json.sha256").write_text(f"{report_hash}  {report_path.name}\n", encoding="ascii")
print(json.dumps(report, indent=2))
print(f"FORMAL_MAPPING_VALIDATION_SHA256={report_hash}")
PY

echo "ARCHIVE=$ARCHIVE"
sha256sum "$ROOT/config/primary_atom_mapping_contacts.json" "$ROOT/config/membrane_qc_mapping.json"
