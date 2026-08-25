#!/usr/bin/env python3
"""Prepare a hash-bound pending visual-adjudication table for O6U water probes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_RETAINED = 20
EXPECTED_TOTAL = 70


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def require_file(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SystemExit(f"Missing or empty input: {resolved}")
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pending-table", required=True, type=Path)
    parser.add_argument("--prescreen", required=True, type=Path)
    parser.add_argument("--generation-report", required=True, type=Path)
    parser.add_argument("--geometry-audit", required=True, type=Path)
    parser.add_argument("--independent-validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", choices=("canary_template", "formal_mp2_template"), required=True)
    args = parser.parse_args()

    pending_path = require_file(args.pending_table)
    prescreen_path = require_file(args.prescreen)
    generation_path = require_file(args.generation_report)
    audit_path = require_file(args.geometry_audit)
    validation_path = require_file(args.independent_validation)
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True)

    pending = json.loads(pending_path.read_text(encoding="utf-8"))
    prescreen = json.loads(prescreen_path.read_text(encoding="utf-8"))
    generation = json.loads(generation_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if (
        pending.get("status") != "frozen_pending_prospective_review"
        or pending.get("orientation_count") != EXPECTED_TOTAL
        or pending.get("pending_count") != EXPECTED_TOTAL
        or pending.get("production_approved") is not False
    ):
        raise SystemExit("Pending disposition registry differs from its frozen state")
    if (
        prescreen.get("status") != "pass_chemical_role_prescreen_visual_review_required"
        or prescreen.get("orientation_count") != EXPECTED_TOTAL
        or prescreen.get("retained_for_visual_review_count") != EXPECTED_RETAINED
        or prescreen.get("production_approved") is not False
    ):
        raise SystemExit("Chemical-role prescreen differs from its exact gate")
    if generation.get("status") != "pass_generation_only_visual_review_required" or generation.get("production_approved") is not False:
        raise SystemExit("Input-generation report differs from its generation-only gate")
    expected_generation_role = "generation_canary" if args.role == "canary_template" else "formal_mp2_target"
    if generation.get("role") != expected_generation_role:
        raise SystemExit("Input-generation role differs from requested adjudication role")
    if (
        audit.get("status") != "pass_geometry_integrity_visual_review_required"
        or audit.get("orientation_count") != EXPECTED_RETAINED
        or audit.get("production_approved") is not False
    ):
        raise SystemExit("Geometry audit differs from its integrity-only gate")
    if (
        validation.get("status") != "pass_geometry_audit_independently_reconstructed"
        or validation.get("orientation_count") != EXPECTED_RETAINED
        or validation.get("production_approved") is not False
    ):
        raise SystemExit("Independent geometry validation differs from its exact gate")
    if validation.get("audit_report", {}).get("sha256") != sha256(audit_path):
        raise SystemExit("Independent validation does not bind the supplied geometry audit")
    if validation.get("generation_report", {}).get("sha256") != sha256(generation_path):
        raise SystemExit("Independent validation does not bind the supplied generation report")
    if audit.get("inputs", {}).get("generation_report", {}).get("sha256") != sha256(generation_path):
        raise SystemExit("Geometry audit does not bind the supplied generation report")

    pending_by_id = {row["orientation_id"]: row for row in pending["orientations"]}
    prescreen_by_id = {row["orientation_id"]: row for row in prescreen["orientations"]}
    audit_by_id = {row["orientation_id"]: row for row in audit["orientations"]}
    validation_by_id = {row["orientation_id"]: row for row in validation["reconstructed_orientations"]}
    retained_ids = sorted(
        orientation_id
        for orientation_id, row in prescreen_by_id.items()
        if row.get("prescreen_suggestion") == "retain_for_visual_review"
    )
    if len(pending_by_id) != EXPECTED_TOTAL or len(prescreen_by_id) != EXPECTED_TOTAL:
        raise SystemExit("Pending registry or prescreen orientation identities are not unique")
    if len(retained_ids) != EXPECTED_RETAINED or set(retained_ids) != set(audit_by_id) or set(retained_ids) != set(validation_by_id):
        raise SystemExit("Retained orientation identities differ across evidence layers")

    rows: list[dict[str, object]] = []
    for orientation_id in retained_ids:
        registry = pending_by_id[orientation_id]
        screened = prescreen_by_id[orientation_id]
        geometry = audit_by_id[orientation_id]
        reconstructed = validation_by_id[orientation_id]
        if not (
            registry.get("source_definition") == screened.get("source_definition") == geometry.get("source_definition")
            and registry.get("target_atom") == screened.get("target_atom") == geometry.get("target_atom")
        ):
            raise SystemExit(f"Orientation identity fields differ for {orientation_id}")
        representative = geometry.get("representative_2p0A", {})
        pdb_record = representative.get("pdb", {})
        if pdb_record.get("sha256") != reconstructed.get("representative_pdb_sha256"):
            raise SystemExit(f"Independent PDB hash differs for {orientation_id}")
        rows.append(
            {
                "orientation_id": orientation_id,
                "source_line_number": registry["source_line_number"],
                "probe_type": registry["probe_type"],
                "target_atom": registry["target_atom"],
                "rotation_degrees": registry["rotation_degrees"],
                "source_definition": registry["source_definition"],
                "representative_distance_angstrom": representative["scan_distance_angstrom"],
                "nearest_non_target_distance_at_2p0A_angstrom": representative["nearest_non_target_distance_angstrom"],
                "nearest_non_target_ligand_atom_at_2p0A": representative["nearest_non_target_ligand_atom"],
                "nearest_non_target_water_atom_at_2p0A": representative["nearest_non_target_water_atom"],
                "minimum_non_target_distance_over_scan_angstrom": geometry["minimum_non_target_distance_over_scan_angstrom"],
                "sanity_collision_anywhere_in_scan": str(bool(geometry["sanity_collision_anywhere_in_scan"])).lower(),
                "representative_pdb_path": pdb_record["path"],
                "representative_pdb_sha256": pdb_record["sha256"],
                "review_decision": "PENDING",
                "review_rationale": "",
                "selection_basis_within_target_atom": "",
                "reviewed_at_utc": "",
                "reviewer_role": "",
            }
        )

    table_path = output_dir / "O6U_WATER_PROBE_VISUAL_ADJUDICATION.tsv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_visual_adjudication_template",
        "status": "pending_visual_adjudication_no_qm_authorized",
        "role": args.role,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "pending_table": artifact(pending_path),
            "prescreen": artifact(prescreen_path),
            "generation_report": artifact(generation_path),
            "geometry_audit": artifact(audit_path),
            "independent_validation": artifact(validation_path),
        },
        "orientation_count": len(rows),
        "pending_count": len(rows),
        "adjudication_table": artifact(table_path),
        "allowed_review_decisions": ["RUN_QM", "EXCLUDE_PRESPECIFIED_STERIC_OR_COMPETING_INTERACTION"],
        "decision_rule": (
            "Within each target atom, retain orientations that maximize the intended interaction and minimize "
            "competing parent-molecule interactions after direct PDB review. Geometry metrics are annotations, "
            "not automatic exclusions."
        ),
        "release_boundary": (
            "This is an unsigned pending template. No row is disposed, and no water-interaction QM, fitting, "
            "CHARMM-GUI build, or MD is authorized."
        ),
    }
    report_path = output_dir / "O6U_WATER_PROBE_VISUAL_ADJUDICATION_TEMPLATE.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "orientation_count": len(rows), "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
