#!/usr/bin/env python3
"""Prepare a hash-bound pending post-QM review table for O6U water probes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


VALIDATION_PASS = "pass_raw_water_qm_independent_numerical_reconstruction"
AUTHORIZATION_PASS = "pass_frozen_preqm_orientation_authorization"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def verify_record(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    if not path.is_file() or path.stat().st_size != record.get("size_bytes") or sha256(path) != record.get("sha256"):
        raise RuntimeError(f"Artifact failed hash/size verification: {label}")
    return path


def atomic_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-report", required=True, type=Path)
    parser.add_argument("--authorization-report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", choices=("formal_postqm_template", "synthetic_canary"), required=True)
    args = parser.parse_args()
    validation_path = args.validation_report.resolve()
    authorization_path = args.authorization_report.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit("Output directory already exists; refusing reuse")

    validation = load_json(validation_path)
    authorization = load_json(authorization_path)
    if (
        validation.get("status") != VALIDATION_PASS
        or validation.get("parameter_fitting_authorized") is not False
        or validation.get("production_approved") is not False
        or validation.get("automatic_scientific_classification_applied") is not False
    ):
        raise SystemExit("Independent numerical validation does not pass its exact gate")
    if verify_record(validation.get("authorization_report"), "validation.authorization_report") != authorization_path:
        raise SystemExit("Validation is bound to a different authorization report")
    selected_ids = validation.get("selected_orientation_ids")
    if not isinstance(selected_ids, list) or not selected_ids or len(selected_ids) != len(set(selected_ids)):
        raise SystemExit("Validated selected orientation universe is missing or duplicated")
    if selected_ids != authorization.get("run_qm_orientation_ids"):
        raise SystemExit("Validated and authorized orientation orders differ")
    if args.role == "formal_postqm_template" and (
        authorization.get("status") != AUTHORIZATION_PASS
        or authorization.get("freeze_role") != "formal_preqm_authorization"
        or authorization.get("water_interaction_qm_authorized") is not True
        or authorization.get("production_approved") is not False
    ):
        raise SystemExit("Formal template requires the exact frozen formal pre-QM authorization")

    minimum_path = verify_record(validation.get("minimum_table"), "validation.minimum_table")
    curve_path = verify_record(validation.get("curve_table"), "validation.curve_table")
    with minimum_path.open("r", encoding="utf-8", newline="") as handle:
        minima = list(csv.DictReader(handle, delimiter="\t"))
    if [row.get("orientation_id") for row in minima] != selected_ids:
        raise SystemExit("Minimum table differs from the selected orientation order")
    authorization_rows = authorization.get("orientations")
    if not isinstance(authorization_rows, list):
        raise SystemExit("Authorization orientation records are missing")
    by_id = {str(row.get("orientation_id")): row for row in authorization_rows if isinstance(row, dict)}
    if len(by_id) != len(authorization_rows) or not set(selected_ids).issubset(by_id):
        raise SystemExit("Authorization orientation identities are incomplete or duplicated")

    rows: list[dict[str, object]] = []
    for minimum in minima:
        orientation_id = str(minimum["orientation_id"])
        distance = float(minimum["minimum_distance_angstrom"])
        if distance == 1.5:
            location = "lower_grid_boundary"
        elif distance == 3.0:
            location = "upper_grid_boundary"
        else:
            location = "grid_interior"
        source = by_id[orientation_id]
        rows.append(
            {
                "orientation_id": orientation_id,
                "probe_type": source.get("probe_type"),
                "target_atom": source.get("target_atom"),
                "source_definition": source.get("source_definition"),
                "minimum_distance_angstrom": minimum["minimum_distance_angstrom"],
                "minimum_interaction_energy_kcal_mol": minimum["minimum_interaction_energy_kcal_mol"],
                "minimum_grid_location_annotation": location,
                "postqm_disposition": "PENDING",
                "review_rationale": "",
                "fitting_target_selection_basis": "",
                "reviewed_at_utc": "",
                "reviewer_role": "",
            }
        )

    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "O6U_WATER_POSTQM_REVIEW_PENDING.tsv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_postqm_review_template",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "postqm_review_template_ready_pending_signed_dispositions",
        "template_role": args.role,
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "automatic_scientific_classification_applied": False,
        "allowed_dispositions": ["APPLICABLE", "WEAK", "UNFAVOURABLE"],
        "validation_report": artifact(validation_path),
        "authorization_report": artifact(authorization_path),
        "curve_table": artifact(curve_path),
        "minimum_table": artifact(minimum_path),
        "orientation_count": len(rows),
        "orientation_ids": selected_ids,
        "pending_table": artifact(table_path),
        "review_rule": (
            "Review every selected orientation using the complete raw curve, chemical role, minimum geometry, "
            "and failed-attempt provenance. Grid-boundary location is an annotation, not an automatic exclusion. "
            "No universal energy threshold is applied."
        ),
        "release_boundary": (
            "This pending template authorizes no fitting. Every row requires a signed APPLICABLE, WEAK, or "
            "UNFAVOURABLE disposition with rationale and within-target selection basis."
        ),
    }
    report_path = output_dir / "O6U_WATER_POSTQM_REVIEW_TEMPLATE.json"
    atomic_json(report_path, report)
    print(json.dumps({"status": report["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
