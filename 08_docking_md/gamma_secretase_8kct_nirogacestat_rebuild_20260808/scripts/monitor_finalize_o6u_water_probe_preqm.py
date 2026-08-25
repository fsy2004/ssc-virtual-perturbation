#!/usr/bin/env python3
"""Wait for formal visual review, then validate and freeze O6U pre-QM actions."""

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


UPSTREAM_READY = "formal_visual_adjudication_review_bundle_ready_pending"
VALIDATION_PASS = "pass_completed_adjudication_independent_integrity_validation"
FREEZE_PASS = "pass_frozen_preqm_orientation_authorization"


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


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def upstream_failed(status: str) -> bool:
    lowered = status.lower()
    return any(token in lowered for token in ("fail", "error", "no_go", "blocked"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.poll_seconds < 10:
        raise SystemExit("Poll interval must be at least 10 seconds")

    scripts = root / "scripts"
    runs = root / "server_runs/o6u_parameterization"
    upstream_controller = runs / "ffparam_water_formal_input_controller_20260811_v4/O6U_FFPARAM_WATER_FORMAL_INPUT_CONTROLLER.json"
    visual_dir = runs / "ffparam_water_visual_adjudication_formal_mp2_20260811_v1"
    template_report = visual_dir / "O6U_WATER_PROBE_VISUAL_ADJUDICATION_TEMPLATE.json"
    completed_tsv = visual_dir / "O6U_WATER_PROBE_VISUAL_ADJUDICATION_COMPLETED.tsv"
    validation_dir = runs / "ffparam_water_completed_adjudication_formal_mp2_20260811_v1"
    validation_report = validation_dir / "O6U_WATER_PROBE_COMPLETED_ADJUDICATION_VALIDATION.json"
    freeze_dir = runs / "ffparam_water_preqm_authorization_formal_mp2_20260811_v1"
    freeze_report = freeze_dir / "O6U_WATER_PROBE_PREQM_AUTHORIZATION.json"
    controller_dir = runs / "ffparam_water_preqm_release_controller_20260811_v1"
    controller_report = controller_dir / "O6U_WATER_PROBE_PREQM_RELEASE_CONTROLLER.json"
    controller_dir.mkdir(parents=True, exist_ok=True)

    lock_handle = (controller_dir / ".controller.lock").open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("Another pre-QM release controller owns the lock") from exc
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_preqm_release_controller",
        "status": "waiting_for_formal_review_bundle",
        "pid": os.getpid(),
        "started_at_utc": now(),
        "updated_at_utc": now(),
        "production_approved": False,
        "water_interaction_qm_authorized": False,
        "upstream_controller": str(upstream_controller),
        "completed_adjudication_tsv_expected": str(completed_tsv),
        "completed_adjudication_validation_report": str(validation_report),
        "preqm_authorization_report": str(freeze_report),
        "release_boundary": (
            "This controller may validate a completed formal visual adjudication and freeze the 70-row pre-QM "
            "orientation authorization. It cannot create review decisions or execute water-interaction QM."
        ),
    }
    atomic_json(controller_report, state)

    while True:
        if upstream_controller.is_file():
            upstream = load_json(upstream_controller)
            status = str(upstream.get("status", ""))
            state["upstream_status"] = status
            if upstream_failed(status):
                state.update({"status": "fail_closed_upstream_formal_review_controller", "updated_at_utc": now()})
                atomic_json(controller_report, state)
                return 2
            if status == UPSTREAM_READY:
                break
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    state.update(
        {
            "status": "waiting_for_completed_visual_adjudication",
            "updated_at_utc": now(),
            "upstream_controller_record": artifact(upstream_controller),
        }
    )
    atomic_json(controller_report, state)
    while not completed_tsv.is_file():
        upstream = load_json(upstream_controller)
        status = str(upstream.get("status", ""))
        if status != UPSTREAM_READY:
            state.update(
                {
                    "status": "fail_closed_upstream_status_changed_after_review_bundle",
                    "updated_at_utc": now(),
                    "upstream_status": status,
                }
            )
            atomic_json(controller_report, state)
            return 2
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    try:
        if validation_dir.exists() or freeze_dir.exists():
            raise RuntimeError("Formal adjudication validation/freeze output path already exists; refusing reuse")
        if not template_report.is_file():
            raise RuntimeError("Formal visual-adjudication template report is missing")

        state.update(
            {
                "status": "validating_completed_visual_adjudication",
                "updated_at_utc": now(),
                "completed_adjudication_tsv": artifact(completed_tsv),
                "visual_adjudication_template_report": artifact(template_report),
            }
        )
        atomic_json(controller_report, state)
        validation_dir.mkdir(parents=True, exist_ok=False)
        validation_command = [
            sys.executable,
            str(scripts / "validate_o6u_water_probe_completed_adjudication.py"),
            "--template-report", str(template_report),
            "--completed-tsv", str(completed_tsv),
            "--report", str(validation_report),
            "--expected-template-role", "formal_mp2_template",
            "--validation-role", "formal_completed_adjudication",
        ]
        validation = subprocess.run(validation_command, cwd=root, text=True, capture_output=True, check=False)
        state["validation_returncode"] = validation.returncode
        state["validation_stdout"] = validation.stdout[-4000:]
        state["validation_stderr"] = validation.stderr[-4000:]
        if validation.returncode != 0:
            raise RuntimeError("Completed formal visual adjudication failed independent integrity validation")
        validation_record = load_json(validation_report)
        if validation_record.get("status") != VALIDATION_PASS or validation_record.get("production_approved") is not False:
            raise RuntimeError("Completed-adjudication validation differs from its exact gate")

        state.update(
            {
                "status": "freezing_70_orientation_preqm_authorization",
                "updated_at_utc": now(),
                "completed_adjudication_validation": artifact(validation_report),
            }
        )
        atomic_json(controller_report, state)
        freeze_dir.mkdir(parents=True, exist_ok=False)
        freeze_command = [
            sys.executable,
            str(scripts / "freeze_o6u_water_probe_preqm_authorization.py"),
            "--pending-registry", str(runs / "water_probe_disposition_policy_20260811_v2/O6U_WATER_PROBE_PROSPECTIVE_DISPOSITIONS_V2.json"),
            "--prescreen", str(runs / "water_probe_disposition_policy_20260811_v2/O6U_WATER_PROBE_CHEMICAL_ROLE_PRESCREEN.json"),
            "--template-report", str(template_report),
            "--completed-tsv", str(completed_tsv),
            "--validation-report", str(validation_report),
            "--report", str(freeze_report),
            "--freeze-role", "formal_preqm_authorization",
        ]
        freeze = subprocess.run(freeze_command, cwd=root, text=True, capture_output=True, check=False)
        state["freeze_returncode"] = freeze.returncode
        state["freeze_stdout"] = freeze.stdout[-4000:]
        state["freeze_stderr"] = freeze.stderr[-4000:]
        if freeze.returncode != 0:
            if freeze_report.is_file() and load_json(freeze_report).get("status") == "scientific_no_go_no_water_qm_orientations_selected":
                state.update(
                    {
                        "status": "scientific_no_go_no_water_qm_orientations_selected",
                        "updated_at_utc": now(),
                        "preqm_authorization": artifact(freeze_report),
                    }
                )
                atomic_json(controller_report, state)
                return 2
            raise RuntimeError("Formal 70-orientation pre-QM freeze failed")
        freeze_record = load_json(freeze_report)
        if (
            freeze_record.get("status") != FREEZE_PASS
            or freeze_record.get("water_interaction_qm_authorized") is not True
            or freeze_record.get("production_approved") is not False
        ):
            raise RuntimeError("Formal pre-QM authorization differs from its exact gate")

        state.update(
            {
                "status": "preqm_orientation_authorization_frozen",
                "updated_at_utc": now(),
                "preqm_authorization": artifact(freeze_report),
                "water_interaction_qm_authorized": True,
                "production_approved": False,
                "release_boundary": (
                    "The 70-orientation pre-QM action universe is frozen. Only run_qm_orientation_ids in the bound "
                    "authorization may be executed. This controller does not run QM or authorize fitting/MD."
                ),
            }
        )
        atomic_json(controller_report, state)
        return 0
    except Exception as exc:
        state.update(
            {
                "status": "fail_closed",
                "updated_at_utc": now(),
                "error": str(exc),
                "water_interaction_qm_authorized": False,
                "production_approved": False,
            }
        )
        atomic_json(controller_report, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
