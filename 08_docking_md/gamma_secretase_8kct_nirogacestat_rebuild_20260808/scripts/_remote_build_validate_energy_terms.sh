set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
PY=/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python
OUT="$ROOT/audit/energy_term_inventory"
mkdir -p "$OUT"
cd "$ROOT"
PYTHONPATH=scripts "$PY" scripts/build_energy_terms_record.py
PYTHONPATH=scripts "$PY" - "$ROOT" "$OUT" <<'PY'
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from gmx_energy_qc import validate_terms_record

root, output = map(Path, sys.argv[1:])
record_path = root / "config" / "gromacs_energy_terms.json"
menu_path = output / "gmx_energy_menu.stderr"
version_path = output / "gmx_version.txt"
edr_path = root / "rep01" / "work" / "production.edr"

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

menu_text = menu_path.read_text(encoding="utf-8", errors="strict")
menu = {int(number): name for number, name in re.findall(r"(?<!\S)(\d+)\s+([^\s]+)", menu_text)}
if set(menu) != set(range(1, 50)):
    raise RuntimeError(f"energy menu indices differ: {sorted(menu)}")

record = json.loads(record_path.read_text(encoding="utf-8"))
required = [term["gmx_name"] for term in record["terms"]]
missing = [name for name in required if name not in set(menu.values())]
if missing:
    raise RuntimeError(f"required GROMACS energy terms are missing: {missing}")

record["approval_status"] = "approved"
record["formal_validation"] = {
    "status": "pass",
    "production_edr": {"path": str(edr_path), "sha256": sha256(edr_path)},
    "gromacs_version_capture": {"path": str(version_path), "sha256": sha256(version_path)},
    "gromacs_energy_menu": {"path": str(menu_path), "sha256": sha256(menu_path), "available_terms": menu},
    "required_terms_present_exactly": required,
    "validated_at_utc": datetime.now(timezone.utc).isoformat(),
}
validate_terms_record(record)
record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
approved_hash = sha256(record_path)
record_path.with_suffix(".json.sha256").write_text(f"{approved_hash}  {record_path.name}\n", encoding="ascii")
validate_terms_record(json.loads(record_path.read_text(encoding="utf-8")))

report = {
    "schema_version": "1.0",
    "report_type": "gromacs_energy_term_record_validation",
    "status": "pass",
    "approved_record": {"path": str(record_path), "sha256": approved_hash},
    "required_term_count": len(required),
    "available_term_count": len(menu),
    "production_edr_sha256": sha256(edr_path),
    "created_at_utc": datetime.now(timezone.utc).isoformat(),
}
report_path = output / "ENERGY_TERM_RECORD_VALIDATION.json"
report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
report_hash = sha256(report_path)
report_path.with_suffix(".json.sha256").write_text(f"{report_hash}  {report_path.name}\n", encoding="ascii")
print(json.dumps(report, indent=2))
print(f"ENERGY_TERM_RECORD_VALIDATION_SHA256={report_hash}")
PY
