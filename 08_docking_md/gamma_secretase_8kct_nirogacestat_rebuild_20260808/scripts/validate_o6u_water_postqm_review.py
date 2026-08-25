#!/usr/bin/env python3
"""Validate a completed O6U post-QM review without judging its chemistry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ALLOWED = {"APPLICABLE", "WEAK", "UNFAVOURABLE"}
IMMUTABLE = (
    "orientation_id", "probe_type", "target_atom", "source_definition",
    "minimum_distance_angstrom", "minimum_interaction_energy_kcal_mol",
    "minimum_grid_location_annotation",
)


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-report", required=True, type=Path)
    parser.add_argument("--completed-tsv", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-template-role", required=True)
    args = parser.parse_args()
    template_path = args.template_report.resolve()
    completed_path = args.completed_tsv.resolve()
    template = load_json(template_path)
    if (
        template.get("status") != "postqm_review_template_ready_pending_signed_dispositions"
        or template.get("template_role") != args.expected_template_role
        or template.get("parameter_fitting_authorized") is not False
        or template.get("automatic_scientific_classification_applied") is not False
    ):
        raise SystemExit("Template report does not pass its exact gate")
    pending_path = verify_record(template.get("pending_table"), "template.pending_table")
    pending = read_tsv(pending_path)
    completed = read_tsv(completed_path)
    if len(pending) != template.get("orientation_count") or len(completed) != len(pending):
        raise SystemExit("Completed review row count differs from the frozen template")
    if [row.get("orientation_id") for row in completed] != template.get("orientation_ids"):
        raise SystemExit("Completed review order differs from the frozen orientation order")

    counts: Counter[str] = Counter()
    for index, (before, after) in enumerate(zip(pending, completed, strict=True), start=1):
        for field in IMMUTABLE:
            if before.get(field) != after.get(field):
                raise SystemExit(f"Immutable field changed on row {index}: {field}")
        disposition = str(after.get("postqm_disposition", "")).strip()
        if disposition not in ALLOWED:
            raise SystemExit(f"Invalid or pending disposition on row {index}")
        for field in ("review_rationale", "fitting_target_selection_basis", "reviewed_at_utc", "reviewer_role"):
            if not str(after.get(field, "")).strip():
                raise SystemExit(f"Missing signed-review field on row {index}: {field}")
        timestamp = datetime.fromisoformat(str(after["reviewed_at_utc"]).replace("Z", "+00:00"))
        if timestamp.tzinfo is None:
            raise SystemExit(f"Review timestamp lacks timezone on row {index}")
        counts[disposition] += 1

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_postqm_completed_review_integrity_validation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_completed_postqm_review_integrity_validation",
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "scientific_correctness_adjudicated": False,
        "automatic_scientific_classification_applied": False,
        "template_report": artifact(template_path),
        "pending_table": artifact(pending_path),
        "completed_table": artifact(completed_path),
        "orientation_count": len(completed),
        "disposition_counts": dict(sorted(counts.items())),
        "release_boundary": (
            "This report validates completeness, identity, immutability, vocabulary, and signatures only. It does "
            "not judge chemical correctness or authorize parameter fitting; a separate fitting-target freeze is required."
        ),
    }
    atomic_json(args.report.resolve(), report)
    print(json.dumps({"status": report["status"], "report": str(args.report.resolve()), "sha256": sha256(args.report.resolve())}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
