#!/usr/bin/env python3
"""Independently validate a completed O6U water-probe adjudication table.

This validator checks provenance and decision-record integrity. It does not
replace the geometry-specific chemical/visual review and cannot authorize QM.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ALLOWED_DECISIONS = {
    "RUN_QM",
    "EXCLUDE_PRESPECIFIED_STERIC_OR_COMPETING_INTERACTION",
}
REVIEW_COLUMNS = {
    "review_decision",
    "review_rationale",
    "selection_basis_within_target_atom",
    "reviewed_at_utc",
    "reviewer_role",
}
EXPECTED_TEMPLATE_STATUS = "pending_visual_adjudication_no_qm_authorized"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"TSV header is missing: {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def parse_review_time(value: str, orientation_id: str) -> datetime:
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RuntimeError(f"Invalid reviewed_at_utc for {orientation_id}: {value}") from exc
    if stamp.tzinfo is None:
        raise RuntimeError(f"reviewed_at_utc lacks timezone for {orientation_id}")
    return stamp.astimezone(timezone.utc)


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-report", required=True, type=Path)
    parser.add_argument("--completed-tsv", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--expected-template-role",
        required=True,
        choices=("canary_template", "formal_mp2_template"),
    )
    parser.add_argument(
        "--validation-role",
        required=True,
        choices=("synthetic_canary_fixture", "formal_completed_adjudication"),
    )
    args = parser.parse_args()

    template_report_path = args.template_report.resolve()
    completed_path = args.completed_tsv.resolve()
    if not template_report_path.is_file() or not completed_path.is_file():
        raise SystemExit("Template report and completed TSV must both exist")

    template = load_json(template_report_path)
    if (
        template.get("status") != EXPECTED_TEMPLATE_STATUS
        or template.get("role") != args.expected_template_role
        or template.get("production_approved") is not False
    ):
        raise SystemExit("Template report does not pass the exact pending-template gate")
    if sorted(template.get("allowed_review_decisions", [])) != sorted(ALLOWED_DECISIONS):
        raise SystemExit("Template allowed-decision set differs from the frozen set")

    source_record = template.get("adjudication_table")
    if not isinstance(source_record, dict):
        raise SystemExit("Template report lacks its source adjudication table record")
    source_path = Path(str(source_record.get("path", ""))).resolve()
    if (
        not source_path.is_file()
        or source_path.stat().st_size != source_record.get("size_bytes")
        or sha256(source_path) != source_record.get("sha256")
    ):
        raise SystemExit("Pending source adjudication table failed hash/size verification")

    source_header, source_rows = read_tsv(source_path)
    completed_header, completed_rows = read_tsv(completed_path)
    if source_header != completed_header:
        raise SystemExit("Completed TSV header or column order differs from the pending template")
    if not REVIEW_COLUMNS.issubset(source_header):
        raise SystemExit("Required review columns are missing")
    if len(source_rows) != template.get("orientation_count") or len(completed_rows) != len(source_rows):
        raise SystemExit("Completed TSV row count differs from the frozen template")

    immutable_columns = [column for column in source_header if column not in REVIEW_COLUMNS]
    source_by_id: dict[str, dict[str, str]] = {}
    for row in source_rows:
        orientation_id = row.get("orientation_id", "")
        if not orientation_id or orientation_id in source_by_id:
            raise SystemExit("Pending template contains blank or duplicate orientation_id")
        source_by_id[orientation_id] = row

    completed_by_id: dict[str, dict[str, str]] = {}
    decision_counts: Counter[str] = Counter()
    target_counts: dict[str, Counter[str]] = defaultdict(Counter)
    reviewers: set[str] = set()
    review_times: list[datetime] = []
    for row in completed_rows:
        orientation_id = row.get("orientation_id", "")
        if not orientation_id or orientation_id in completed_by_id:
            raise SystemExit("Completed table contains blank or duplicate orientation_id")
        if orientation_id not in source_by_id:
            raise SystemExit(f"Unexpected orientation_id in completed table: {orientation_id}")
        source = source_by_id[orientation_id]
        for column in immutable_columns:
            if row.get(column, "") != source.get(column, ""):
                raise SystemExit(f"Immutable field changed for {orientation_id}: {column}")

        decision = row.get("review_decision", "").strip()
        rationale = row.get("review_rationale", "").strip()
        selection_basis = row.get("selection_basis_within_target_atom", "").strip()
        reviewed_at = row.get("reviewed_at_utc", "").strip()
        reviewer_role = row.get("reviewer_role", "").strip()
        if decision not in ALLOWED_DECISIONS:
            raise SystemExit(f"Unresolved or invalid review decision for {orientation_id}: {decision}")
        if not rationale or not selection_basis or not reviewed_at or not reviewer_role:
            raise SystemExit(f"Incomplete review fields for {orientation_id}")
        if args.validation_role == "formal_completed_adjudication" and any(
            token in reviewer_role.lower() for token in ("synthetic", "fixture", "canary", "test")
        ):
            raise SystemExit(f"Formal adjudication has a non-formal reviewer role: {orientation_id}")
        review_time = parse_review_time(reviewed_at, orientation_id)
        if review_time > datetime.now(timezone.utc):
            raise SystemExit(f"Review timestamp is in the future: {orientation_id}")

        decision_counts[decision] += 1
        target_counts[row.get("target_atom", "")][decision] += 1
        reviewers.add(reviewer_role)
        review_times.append(review_time)
        completed_by_id[orientation_id] = row

    if set(completed_by_id) != set(source_by_id):
        raise SystemExit("Completed orientation_id set differs from the pending template")

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_completed_adjudication_independent_validation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_completed_adjudication_independent_integrity_validation",
        "validation_role": args.validation_role,
        "expected_template_role": args.expected_template_role,
        "production_approved": False,
        "template_report": artifact(template_report_path),
        "pending_template_tsv": artifact(source_path),
        "completed_adjudication_tsv": artifact(completed_path),
        "orientation_count": len(completed_rows),
        "decision_counts": dict(sorted(decision_counts.items())),
        "target_atom_decision_counts": {
            target: dict(sorted(counts.items())) for target, counts in sorted(target_counts.items())
        },
        "reviewer_roles": sorted(reviewers),
        "review_time_range_utc": {
            "minimum": min(review_times).isoformat(),
            "maximum": max(review_times).isoformat(),
        },
        "automatic_exclusion_applied": False,
        "scientific_review_boundary": (
            "This report independently validates source identity, immutable geometry fields, decision vocabulary, "
            "and review-record completeness. It cannot determine whether the chemical/visual decisions are correct."
        ),
        "release_boundary": (
            "This integrity validation alone does not authorize water-interaction QM, parameter fitting, "
            "CHARMM-GUI construction, or MD. A separately frozen formal disposition record must bind the completed "
            "table, this validation, and the reviewer-asserted geometry-specific decision basis."
        ),
    }
    atomic_json(args.report.resolve(), report)
    print(json.dumps({"status": report["status"], "report": str(args.report.resolve()), "sha256": sha256(args.report.resolve())}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
