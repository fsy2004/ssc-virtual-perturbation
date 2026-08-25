#!/usr/bin/env python3
"""Freeze the complete 70-orientation O6U pre-QM authorization record."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


VALIDATION_STATUS = "pass_completed_adjudication_independent_integrity_validation"
TEMPLATE_STATUS = "pending_visual_adjudication_no_qm_authorized"
RUN_QM = "RUN_QM"
VISUAL_EXCLUDE = "EXCLUDE_PRESPECIFIED_STERIC_OR_COMPETING_INTERACTION"


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Missing TSV header: {path}")
        return [dict(row) for row in reader]


def verify_recorded_artifact(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    if (
        not path.is_file()
        or path.stat().st_size != record.get("size_bytes")
        or sha256(path) != record.get("sha256")
    ):
        raise RuntimeError(f"Artifact failed size/hash verification: {label}")
    return path


def keyed(rows: list[dict[str, object]], label: str) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        orientation_id = str(row.get("orientation_id", ""))
        if not orientation_id or orientation_id in result:
            raise RuntimeError(f"Blank or duplicate orientation_id in {label}")
        result[orientation_id] = row
    return result


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pending-registry", required=True, type=Path)
    parser.add_argument("--prescreen", required=True, type=Path)
    parser.add_argument("--template-report", required=True, type=Path)
    parser.add_argument("--completed-tsv", required=True, type=Path)
    parser.add_argument("--validation-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--freeze-role",
        required=True,
        choices=("synthetic_canary_fixture", "formal_preqm_authorization"),
    )
    args = parser.parse_args()

    paths = {
        "pending_registry": args.pending_registry.resolve(),
        "prescreen": args.prescreen.resolve(),
        "template_report": args.template_report.resolve(),
        "completed_tsv": args.completed_tsv.resolve(),
        "validation_report": args.validation_report.resolve(),
    }
    if any(not path.is_file() for path in paths.values()):
        raise SystemExit("Every declared input must exist")

    registry = load_json(paths["pending_registry"])
    prescreen = load_json(paths["prescreen"])
    template = load_json(paths["template_report"])
    validation = load_json(paths["validation_report"])
    completed_rows = read_tsv(paths["completed_tsv"])

    expected_template_role = "canary_template" if args.freeze_role == "synthetic_canary_fixture" else "formal_mp2_template"
    expected_validation_role = "synthetic_canary_fixture" if args.freeze_role == "synthetic_canary_fixture" else "formal_completed_adjudication"
    if (
        template.get("status") != TEMPLATE_STATUS
        or template.get("role") != expected_template_role
        or template.get("production_approved") is not False
    ):
        raise SystemExit("Visual-adjudication template does not pass its exact role/status gate")
    if (
        validation.get("status") != VALIDATION_STATUS
        or validation.get("validation_role") != expected_validation_role
        or validation.get("production_approved") is not False
    ):
        raise SystemExit("Completed-adjudication validation does not pass its exact role/status gate")

    if verify_recorded_artifact(validation.get("template_report"), "validation.template_report") != paths["template_report"]:
        raise SystemExit("Validation is bound to a different template report")
    if verify_recorded_artifact(validation.get("completed_adjudication_tsv"), "validation.completed_adjudication_tsv") != paths["completed_tsv"]:
        raise SystemExit("Validation is bound to a different completed TSV")
    template_inputs = template.get("inputs")
    if not isinstance(template_inputs, dict):
        raise SystemExit("Template input map is missing")
    if verify_recorded_artifact(template_inputs.get("pending_table"), "template.pending_table") != paths["pending_registry"]:
        raise SystemExit("Template is bound to a different pending registry")
    if verify_recorded_artifact(template_inputs.get("prescreen"), "template.prescreen") != paths["prescreen"]:
        raise SystemExit("Template is bound to a different prescreen")

    registry_rows = registry.get("orientations")
    prescreen_rows = prescreen.get("orientations")
    if not isinstance(registry_rows, list) or not isinstance(prescreen_rows, list):
        raise SystemExit("Registry or prescreen orientations are missing")
    if registry.get("orientation_count") != 70 or prescreen.get("orientation_count") != 70:
        raise SystemExit("The prospective 70-orientation universe is not intact")

    registry_by_id = keyed(registry_rows, "pending registry")
    prescreen_by_id = keyed(prescreen_rows, "prescreen")
    completed_by_id = keyed(completed_rows, "completed visual adjudication")
    if set(registry_by_id) != set(prescreen_by_id):
        raise SystemExit("Registry and prescreen orientation universes differ")

    retained_ids = {
        orientation_id
        for orientation_id, row in prescreen_by_id.items()
        if row.get("prescreen_suggestion") == "retain_for_visual_review"
    }
    if len(retained_ids) != 20 or set(completed_by_id) != retained_ids:
        raise SystemExit("Completed visual adjudication does not cover the exact 20 retained orientations")

    frozen_rows: list[dict[str, object]] = []
    action_counts: Counter[str] = Counter()
    for orientation_id in sorted(registry_by_id, key=lambda x: int(x.rsplit("_", 1)[1])):
        registry_row = registry_by_id[orientation_id]
        prescreen_row = prescreen_by_id[orientation_id]
        base = {
            "orientation_id": orientation_id,
            "source_line_number": registry_row.get("source_line_number"),
            "source_definition": registry_row.get("source_definition"),
            "probe_type": registry_row.get("probe_type"),
            "target_atom": registry_row.get("target_atom"),
        }
        for field in ("source_line_number", "source_definition", "probe_type", "target_atom"):
            if base[field] != prescreen_row.get(field):
                raise SystemExit(f"Registry/prescreen identity mismatch for {orientation_id}: {field}")

        if orientation_id in retained_ids:
            completed = completed_by_id[orientation_id]
            decision = str(completed.get("review_decision", ""))
            if decision not in (RUN_QM, VISUAL_EXCLUDE):
                raise SystemExit(f"Invalid completed visual decision for {orientation_id}")
            action = "run_hf_631gd_water_interaction_qm" if decision == RUN_QM else "exclude_prespecified_steric_or_competing_interaction"
            evidence = {
                "decision_source": "geometry_specific_visual_adjudication",
                "review_decision": decision,
                "review_rationale": completed.get("review_rationale"),
                "selection_basis_within_target_atom": completed.get("selection_basis_within_target_atom"),
                "reviewed_at_utc": completed.get("reviewed_at_utc"),
                "reviewer_role": completed.get("reviewer_role"),
                "representative_pdb_path": completed.get("representative_pdb_path"),
                "representative_pdb_sha256": completed.get("representative_pdb_sha256"),
            }
        else:
            if prescreen_row.get("final_disposition") != "pending_signed_exclusion":
                raise SystemExit(f"Non-retained prescreen row lacks pending signed exclusion: {orientation_id}")
            rationale = str(prescreen_row.get("prescreen_rationale", "")).strip()
            suggestion = str(prescreen_row.get("prescreen_suggestion", "")).strip()
            if not rationale or not suggestion or suggestion == "retain_for_visual_review":
                raise SystemExit(f"Non-retained prescreen evidence is incomplete: {orientation_id}")
            action = "exclude_prespecified_chemical_role"
            evidence = {
                "decision_source": "prospective_chemical_role_prescreen",
                "interaction_role": prescreen_row.get("interaction_role"),
                "prescreen_suggestion": suggestion,
                "prescreen_rationale": rationale,
            }
        action_counts[action] += 1
        frozen_rows.append({**base, "pre_qm_action": action, "decision_evidence": evidence})

    run_ids = [row["orientation_id"] for row in frozen_rows if row["pre_qm_action"] == "run_hf_631gd_water_interaction_qm"]
    formal = args.freeze_role == "formal_preqm_authorization"
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_preqm_orientation_authorization",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_frozen_preqm_orientation_authorization" if formal and run_ids else (
            "scientific_no_go_no_water_qm_orientations_selected" if formal else "pass_synthetic_canary_structure_only"
        ),
        "freeze_role": args.freeze_role,
        "production_approved": False,
        "water_interaction_qm_authorized": bool(formal and run_ids),
        "inputs": {name: artifact(path) for name, path in paths.items()},
        "orientation_count": len(frozen_rows),
        "pre_qm_action_counts": dict(sorted(action_counts.items())),
        "run_qm_orientation_ids": run_ids,
        "orientations": frozen_rows,
        "automatic_distance_exclusion_applied": False,
        "release_boundary": (
            "Only the explicitly listed run_qm_orientation_ids may enter retained HF/6-31G(d) water-interaction "
            "calculations, and only when water_interaction_qm_authorized is true. This record does not authorize "
            "parameter fitting, CHARMM-GUI construction, or MD. Failed QM attempts must remain in provenance."
        ),
    }
    atomic_json(args.report.resolve(), report)
    print(json.dumps({"status": report["status"], "run_qm_count": len(run_ids), "report": str(args.report.resolve()), "sha256": sha256(args.report.resolve())}, sort_keys=True))
    return 0 if not formal or run_ids else 2


if __name__ == "__main__":
    raise SystemExit(main())
