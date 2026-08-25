#!/usr/bin/env python3
"""Wait for the frozen O6U MP2 ensemble, then generate and audit formal water inputs."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ENSEMBLE_STATUS = "pass_five_member_ensemble_independently_reconstructed"
GENERATION_STATUS = "pass_generation_only_visual_review_required"
AUDIT_STATUS = "pass_geometry_integrity_visual_review_required"
INDEPENDENT_VALIDATION_STATUS = "pass_geometry_audit_independently_reconstructed"
VISUAL_TEMPLATE_STATUS = "pending_visual_adjudication_no_qm_authorized"
RENDERING_STATUS = "pass_rendering_only_direct_review_required"
RENDERING_VALIDATION_STATUS = "pass_rendering_independent_geometry_reconstruction"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.poll_seconds < 10:
        raise SystemExit("Poll interval must be at least 10 seconds")

    scripts = root / "scripts"
    ensemble_report = root / "server_runs/o6u_parameterization/mp2_631gd_ensemble_freeze_20260811_v1/O6U_MP2_REPRESENTATIVE_ENSEMBLE_FREEZE.json"
    ensemble_controller = root / "server_runs/o6u_parameterization/mp2_631gd_ensemble_freeze_controller_20260811_v1/O6U_MP2_ENSEMBLE_FREEZE_CONTROLLER.json"
    formal_dir = root / "server_runs/o6u_parameterization/ffparam_water_input_generation_formal_mp2_20260811_v1"
    audit_dir = root / "server_runs/o6u_parameterization/ffparam_water_geometry_audit_formal_mp2_20260811_v1"
    validation_dir = root / "server_runs/o6u_parameterization/ffparam_water_geometry_audit_formal_mp2_independent_validation_20260811_v1"
    visual_template_dir = root / "server_runs/o6u_parameterization/ffparam_water_visual_adjudication_formal_mp2_20260811_v1"
    review_rendering_dir = root / "server_runs/o6u_parameterization/ffparam_water_review_rendering_formal_mp2_20260811_v1"
    controller_dir = root / "server_runs/o6u_parameterization/ffparam_water_formal_input_controller_20260811_v4"
    controller_report = controller_dir / "O6U_FFPARAM_WATER_FORMAL_INPUT_CONTROLLER.json"
    controller_dir.mkdir(parents=True, exist_ok=True)

    lock_path = controller_dir / ".controller.lock"
    lock_handle = lock_path.open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("Another formal-input controller owns the lock") from exc
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_ffparam_water_formal_input_controller",
        "status": "waiting_for_frozen_mp2_ensemble",
        "production_approved": False,
        "pid": os.getpid(),
        "started_at_utc": now(),
        "updated_at_utc": now(),
        "ensemble_report_expected": str(ensemble_report),
        "formal_generation_directory": str(formal_dir),
        "formal_geometry_audit_directory": str(audit_dir),
        "formal_independent_validation_directory": str(validation_dir),
        "formal_visual_adjudication_template_directory": str(visual_template_dir),
        "formal_review_rendering_directory": str(review_rendering_dir),
        "release_boundary": (
            "This controller may generate inputs, audit geometry, and independently reconstruct that audit only; "
            "it cannot execute water-interaction QM."
        ),
    }
    atomic_json(controller_report, state)

    while not ensemble_report.is_file():
        if ensemble_controller.is_file():
            upstream = load_json(ensemble_controller)
            upstream_status = str(upstream.get("status", ""))
            state["upstream_controller_status"] = upstream_status
            if any(token in upstream_status.lower() for token in ("fail", "error", "no_go", "blocked")):
                state.update({"status": "fail_closed_upstream_ensemble_controller", "updated_at_utc": now()})
                atomic_json(controller_report, state)
                return 2
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    try:
        ensemble = load_json(ensemble_report)
        if ensemble.get("status") != ENSEMBLE_STATUS or ensemble.get("production_approved") is not False:
            raise RuntimeError("Frozen MP2 ensemble report does not pass its exact release gate")
        geometry_record = ensemble.get("charge_water_target_optimized_xyz")
        if not isinstance(geometry_record, dict):
            raise RuntimeError("Frozen ensemble lacks its charge/water target geometry")
        geometry = Path(str(geometry_record.get("path", ""))).resolve()
        if not geometry.is_file() or sha256(geometry) != geometry_record.get("sha256"):
            raise RuntimeError("Frozen charge/water target geometry integrity differs")
        if (
            formal_dir.exists()
            or audit_dir.exists()
            or validation_dir.exists()
            or visual_template_dir.exists()
            or review_rendering_dir.exists()
        ):
            raise RuntimeError("Formal generation/audit/validation/adjudication/rendering output path already exists; refusing reuse")

        state.update(
            {
                "status": "generating_formal_inputs",
                "updated_at_utc": now(),
                "ensemble_report": artifact(ensemble_report),
                "charge_water_target_geometry": artifact(geometry),
            }
        )
        atomic_json(controller_report, state)
        generation_command = [
            sys.executable,
            str(scripts / "generate_o6u_ffparam_water_inputs.py"),
            "--source-sdf", str(root / "inputs/ligand_parameterization/O6U_neutral_hydrogen_complete_3D.sdf"),
            "--correspondence-tsv", str(root / "inputs/ligand_parameterization/O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.server.tsv"),
            "--rtf", str(root / "server_records/charmmgui_ligand_reader/official_outputs/o6u__o6u.rtf"),
            "--coordinate-template", str(root / "server_records/charmmgui_ligand_reader/official_outputs/ligandrm.crd"),
            "--geometry-xyz", str(geometry),
            "--orientation-da", str(root / "server_runs/o6u_parameterization/water_probe_orientation_plan_20260811_v1/O6U_WATER_PROBE_ORIENTATIONS.da"),
            "--prescreen", str(root / "server_runs/o6u_parameterization/water_probe_disposition_policy_20260811_v2/O6U_WATER_PROBE_CHEMICAL_ROLE_PRESCREEN.json"),
            "--policy", str(root / "O6U_WATER_PROBE_DISPOSITION_POLICY.json"),
            "--output-dir", str(formal_dir),
            "--role", "formal_mp2_target",
            "--ensemble-report", str(ensemble_report),
        ]
        generation = subprocess.run(generation_command, cwd=root, text=True, capture_output=True, check=False)
        state["generation_returncode"] = generation.returncode
        state["generation_stdout"] = generation.stdout[-4000:]
        state["generation_stderr"] = generation.stderr[-4000:]
        if generation.returncode != 0:
            raise RuntimeError("Formal FFParam input generation failed")
        generation_report = formal_dir / "O6U_FFPARAM_WATER_INPUT_GENERATION.json"
        generation_record = load_json(generation_report)
        if generation_record.get("status") != GENERATION_STATUS or generation_record.get("role") != "formal_mp2_target":
            raise RuntimeError("Formal generation report differs from its exact gate")

        state.update({"status": "auditing_formal_input_geometry", "updated_at_utc": now(), "generation_report": artifact(generation_report)})
        atomic_json(controller_report, state)
        audit_command = [
            sys.executable,
            str(scripts / "audit_o6u_ffparam_water_geometry.py"),
            "--generation-report", str(generation_report),
            "--rtf", str(root / "server_records/charmmgui_ligand_reader/official_outputs/o6u__o6u.rtf"),
            "--output-dir", str(audit_dir),
        ]
        audit = subprocess.run(audit_command, cwd=root, text=True, capture_output=True, check=False)
        state["audit_returncode"] = audit.returncode
        state["audit_stdout"] = audit.stdout[-4000:]
        state["audit_stderr"] = audit.stderr[-4000:]
        if audit.returncode != 0:
            raise RuntimeError("Formal water-input geometry audit failed")
        audit_report = audit_dir / "O6U_FFPARAM_WATER_GEOMETRY_AUDIT.json"
        audit_record = load_json(audit_report)
        if audit_record.get("status") != AUDIT_STATUS or audit_record.get("production_approved") is not False:
            raise RuntimeError("Formal geometry audit report differs from its exact gate")

        state.update(
            {
                "status": "independently_validating_formal_geometry_audit",
                "updated_at_utc": now(),
                "audit_report": artifact(audit_report),
            }
        )
        atomic_json(controller_report, state)
        validation_dir.mkdir(parents=True, exist_ok=False)
        validation_report = validation_dir / "O6U_FFPARAM_WATER_GEOMETRY_AUDIT_INDEPENDENT_VALIDATION.json"
        validation_command = [
            sys.executable,
            str(scripts / "validate_o6u_ffparam_water_geometry_audit.py"),
            "--audit-report", str(audit_report),
            "--generation-report", str(generation_report),
            "--rtf", str(root / "server_records/charmmgui_ligand_reader/official_outputs/o6u__o6u.rtf"),
            "--report", str(validation_report),
        ]
        validation = subprocess.run(validation_command, cwd=root, text=True, capture_output=True, check=False)
        state["independent_validation_returncode"] = validation.returncode
        state["independent_validation_stdout"] = validation.stdout[-4000:]
        state["independent_validation_stderr"] = validation.stderr[-4000:]
        if validation.returncode != 0:
            raise RuntimeError("Independent formal water-geometry validation failed")
        validation_record = load_json(validation_report)
        if (
            validation_record.get("status") != INDEPENDENT_VALIDATION_STATUS
            or validation_record.get("production_approved") is not False
        ):
            raise RuntimeError("Independent formal geometry validation report differs from its exact gate")

        state.update(
            {
                "status": "preparing_formal_visual_adjudication_template",
                "updated_at_utc": now(),
                "independent_validation_report": artifact(validation_report),
            }
        )
        atomic_json(controller_report, state)
        visual_command = [
            sys.executable,
            str(scripts / "prepare_o6u_water_probe_visual_adjudication.py"),
            "--pending-table", str(root / "server_runs/o6u_parameterization/water_probe_disposition_policy_20260811_v2/O6U_WATER_PROBE_PROSPECTIVE_DISPOSITIONS_V2.json"),
            "--prescreen", str(root / "server_runs/o6u_parameterization/water_probe_disposition_policy_20260811_v2/O6U_WATER_PROBE_CHEMICAL_ROLE_PRESCREEN.json"),
            "--generation-report", str(generation_report),
            "--geometry-audit", str(audit_report),
            "--independent-validation", str(validation_report),
            "--output-dir", str(visual_template_dir),
            "--role", "formal_mp2_template",
        ]
        visual = subprocess.run(visual_command, cwd=root, text=True, capture_output=True, check=False)
        state["visual_template_returncode"] = visual.returncode
        state["visual_template_stdout"] = visual.stdout[-4000:]
        state["visual_template_stderr"] = visual.stderr[-4000:]
        if visual.returncode != 0:
            raise RuntimeError("Formal visual-adjudication template preparation failed")
        visual_report = visual_template_dir / "O6U_WATER_PROBE_VISUAL_ADJUDICATION_TEMPLATE.json"
        visual_record = load_json(visual_report)
        if (
            visual_record.get("status") != VISUAL_TEMPLATE_STATUS
            or visual_record.get("role") != "formal_mp2_template"
            or visual_record.get("production_approved") is not False
        ):
            raise RuntimeError("Formal visual-adjudication template differs from its exact gate")

        state.update(
            {
                "status": "rendering_formal_geometry_review_bundle",
                "updated_at_utc": now(),
                "visual_adjudication_template_report": artifact(visual_report),
            }
        )
        atomic_json(controller_report, state)
        rendering_command = [
            sys.executable,
            str(scripts / "render_o6u_water_probe_review_panels.py"),
            "--template-report", str(visual_report),
            "--source-sdf", str(root / "inputs/ligand_parameterization/O6U_neutral_hydrogen_complete_3D.sdf"),
            "--output-dir", str(review_rendering_dir),
            "--role", "formal_geometry_review",
        ]
        rendering = subprocess.run(rendering_command, cwd=root, text=True, capture_output=True, check=False)
        state["review_rendering_returncode"] = rendering.returncode
        state["review_rendering_stdout"] = rendering.stdout[-4000:]
        state["review_rendering_stderr"] = rendering.stderr[-4000:]
        if rendering.returncode != 0:
            raise RuntimeError("Formal geometry-review rendering failed")
        rendering_report = review_rendering_dir / "O6U_WATER_PROBE_GEOMETRY_REVIEW_RENDERING.json"
        rendering_record = load_json(rendering_report)
        if (
            rendering_record.get("status") != RENDERING_STATUS
            or rendering_record.get("role") != "formal_geometry_review"
            or rendering_record.get("production_approved") is not False
        ):
            raise RuntimeError("Formal geometry-review rendering differs from its exact gate")

        state.update(
            {
                "status": "independently_validating_formal_review_rendering",
                "updated_at_utc": now(),
                "review_rendering_report": artifact(rendering_report),
            }
        )
        atomic_json(controller_report, state)
        rendering_validation_report = review_rendering_dir / "O6U_WATER_PROBE_GEOMETRY_REVIEW_RENDERING_INDEPENDENT_VALIDATION.json"
        rendering_validation_command = [
            sys.executable,
            str(scripts / "validate_o6u_water_probe_review_rendering.py"),
            "--rendering-report", str(rendering_report),
            "--report", str(rendering_validation_report),
            "--expected-role", "formal_geometry_review",
        ]
        rendering_validation = subprocess.run(
            rendering_validation_command, cwd=root, text=True, capture_output=True, check=False
        )
        state["review_rendering_validation_returncode"] = rendering_validation.returncode
        state["review_rendering_validation_stdout"] = rendering_validation.stdout[-4000:]
        state["review_rendering_validation_stderr"] = rendering_validation.stderr[-4000:]
        if rendering_validation.returncode != 0:
            raise RuntimeError("Independent formal geometry-review rendering validation failed")
        rendering_validation_record = load_json(rendering_validation_report)
        if (
            rendering_validation_record.get("status") != RENDERING_VALIDATION_STATUS
            or rendering_validation_record.get("production_approved") is not False
        ):
            raise RuntimeError("Independent formal rendering validation differs from its exact gate")
        state.update(
            {
                "status": "formal_visual_adjudication_review_bundle_ready_pending",
                "updated_at_utc": now(),
                "audit_report": artifact(audit_report),
                "independent_validation_report": artifact(validation_report),
                "visual_adjudication_template_report": artifact(visual_report),
                "review_rendering_report": artifact(rendering_report),
                "review_rendering_independent_validation_report": artifact(rendering_validation_report),
                "production_approved": False,
                "release_boundary": (
                    "Formal inputs are generated, geometrically audited, independently reconstructed, bound to a "
                    "pending visual-adjudication template, and rendered into an independently geometry-validated review "
                    "bundle. No water-interaction QM may run until geometry-specific visual/chemical dispositions are "
                    "completed, frozen, and independently validated."
                ),
            }
        )
        atomic_json(controller_report, state)
        return 0
    except Exception as exc:
        state.update({"status": "fail_closed", "updated_at_utc": now(), "error": str(exc), "production_approved": False})
        atomic_json(controller_report, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
