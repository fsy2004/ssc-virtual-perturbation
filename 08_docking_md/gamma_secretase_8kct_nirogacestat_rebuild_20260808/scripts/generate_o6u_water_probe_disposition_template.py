#!/usr/bin/env python3
"""Freeze an auditable pending-review table for every O6U water probe."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_COUNT = 70
ALLOWED_PROBE_TYPES = {"A2", "A31", "AP", "APL", "D", "DOP"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orientation-da", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    source = args.orientation_da.resolve()
    policy = args.policy.resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise SystemExit(f"Missing or empty orientation file: {source}")
    if not policy.is_file() or policy.stat().st_size == 0:
        raise SystemExit(f"Missing or empty disposition policy: {policy}")
    policy_payload = json.loads(policy.read_text(encoding="utf-8"))
    if (
        policy_payload.get("schema_version") != "1.0"
        or policy_payload.get("report_type") != "o6u_water_probe_disposition_policy"
        or policy_payload.get("status") != "pass"
        or policy_payload.get("production_approved") is not False
    ):
        raise SystemExit("Disposition policy is not a valid fail-closed policy record")
    raw_lines = [line.strip() for line in source.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(raw_lines) != EXPECTED_COUNT:
        raise SystemExit(f"Expected {EXPECTED_COUNT} orientations, found {len(raw_lines)}")

    rows: list[dict[str, object]] = []
    for index, raw_line in enumerate(raw_lines, start=1):
        fields = raw_line.split()
        if len(fields) < 5 or fields[0] not in ALLOWED_PROBE_TYPES:
            raise SystemExit(f"Malformed orientation line {index}: {raw_line}")
        try:
            rotation_degrees = float(fields[1])
        except ValueError as exc:
            raise SystemExit(f"Non-numeric rotation on line {index}: {raw_line}") from exc
        rows.append(
            {
                "orientation_id": f"O6U_WP_{index:03d}",
                "source_line_number": index,
                "source_definition": raw_line,
                "probe_type": fields[0],
                "rotation_degrees": rotation_degrees,
                "target_atom": fields[2],
                "reference_atoms": fields[3:],
                "prospective_disposition": "pending_review",
                "allowed_final_dispositions": ["applicable", "excluded", "weak", "unfavourable"],
                "disposition_rationale": None,
                "selection_basis": None,
                "disposition_evidence_artifact": None,
                "reviewer": None,
                "reviewed_at_utc": None,
                "allowed_hf_631gd_status_by_final_disposition": {
                    "applicable": ["pass"],
                    "weak": ["pass"],
                    "unfavourable": ["pass"],
                    "excluded": ["not_required_prespecified_exclusion"],
                },
                "hf_631gd_distance_optimization_status": "not_started",
                "qm_input_artifact": None,
                "qm_output_artifact": None,
                "failed_attempt_artifacts": [],
            }
        )

    output = args.output.resolve()
    report = {
        "schema_version": "1.1",
        "report_type": "o6u_water_probe_prospective_disposition_table",
        "status": "frozen_pending_prospective_review",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "component_id": "O6U",
        "source_orientation_da": {
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "sha256": sha256(source),
        },
        "disposition_policy": {
            "path": str(policy),
            "size_bytes": policy.stat().st_size,
            "sha256": sha256(policy),
        },
        "orientation_count": len(rows),
        "pending_count": len(rows),
        "final_disposition_counts": {
            "applicable": 0,
            "excluded": 0,
            "weak": 0,
            "unfavourable": 0,
        },
        "decision_rule": (
            "Every generated orientation is registered before QM execution. Applicable, weak, and "
            "unfavourable rows require converged retained HF/6-31G(d) input/output. Prespecified "
            "exclusions require retained visual or chemical evidence and no QM artifacts. Failed "
            "attempts cannot be silently omitted or used as targets."
        ),
        "orientations": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "output": str(output), "sha256": sha256(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
