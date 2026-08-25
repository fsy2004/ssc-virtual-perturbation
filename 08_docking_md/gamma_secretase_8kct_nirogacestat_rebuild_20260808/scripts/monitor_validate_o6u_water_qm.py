#!/usr/bin/env python3
"""Wait for formal raw O6U water QM, then independently reconstruct curves."""

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


UPSTREAM_READY = "raw_water_qm_batch_complete_independent_validation_required"
VALIDATION_PASS = "pass_raw_water_qm_independent_numerical_reconstruction"


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
    execution_controller = runs / (
        "water_interaction_qm_execution_controller_20260811_v1/"
        "O6U_WATER_INTERACTION_QM_EXECUTION_CONTROLLER.json"
    )
    batch_report = runs / (
        "water_interaction_qm_batch_formal_mp2_20260811_v1/"
        "O6U_WATER_INTERACTION_QM_BATCH.json"
    )
    generation_report = runs / (
        "ffparam_water_input_generation_formal_mp2_20260811_v1/"
        "O6U_FFPARAM_WATER_INPUT_GENERATION.json"
    )
    authorization_report = runs / (
        "ffparam_water_preqm_authorization_formal_mp2_20260811_v1/"
        "O6U_WATER_PROBE_PREQM_AUTHORIZATION.json"
    )
    validation_dir = runs / "water_interaction_qm_independent_validation_formal_mp2_20260811_v1"
    validation_report = validation_dir / "O6U_WATER_INTERACTION_QM_INDEPENDENT_VALIDATION.json"
    controller_dir = runs / "water_interaction_qm_validation_controller_20260811_v1"
    controller_report = controller_dir / "O6U_WATER_INTERACTION_QM_VALIDATION_CONTROLLER.json"
    controller_dir.mkdir(parents=True, exist_ok=True)

    lock_handle = (controller_dir / ".controller.lock").open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("Another water-QM validation controller owns the lock") from exc
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_water_interaction_qm_validation_controller",
        "status": "waiting_for_formal_raw_water_qm_batch",
        "pid": os.getpid(),
        "started_at_utc": now(),
        "updated_at_utc": now(),
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "automatic_scientific_classification_applied": False,
        "execution_controller": str(execution_controller),
        "formal_batch_report": str(batch_report),
        "formal_generation_report": str(generation_report),
        "formal_authorization_report": str(authorization_report),
        "independent_validation_directory": str(validation_dir),
        "release_boundary": (
            "This controller may independently reconstruct raw no-BSSE HF/6-31G(d) curves only. It cannot "
            "classify orientations, authorize parameter fitting, construct a CHARMM-GUI system, or run MD."
        ),
    }
    atomic_json(controller_report, state)

    while True:
        if execution_controller.is_file():
            upstream = load_json(execution_controller)
            status = str(upstream.get("status", ""))
            state["upstream_status"] = status
            if upstream_failed(status):
                state.update({"status": "fail_closed_upstream_water_qm_execution", "updated_at_utc": now()})
                atomic_json(controller_report, state)
                return 2
            if status == UPSTREAM_READY:
                break
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    try:
        if validation_dir.exists():
            raise RuntimeError("Formal independent-validation directory already exists; refusing reuse")
        for required in (batch_report, generation_report, authorization_report):
            if not required.is_file():
                raise RuntimeError(f"Required formal input is missing: {required}")
        state.update(
            {
                "status": "reconstructing_formal_raw_water_qm_numerically",
                "updated_at_utc": now(),
                "execution_controller_record": artifact(execution_controller),
                "batch_report_record": artifact(batch_report),
                "generation_report_record": artifact(generation_report),
                "authorization_report_record": artifact(authorization_report),
            }
        )
        atomic_json(controller_report, state)
        command = [
            sys.executable,
            str(scripts / "validate_o6u_water_interaction_qm.py"),
            "--batch-report", str(batch_report),
            "--generation-report", str(generation_report),
            "--authorization-report", str(authorization_report),
            "--output-dir", str(validation_dir),
        ]
        run = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
        state["validator_returncode"] = run.returncode
        state["validator_stdout"] = run.stdout[-4000:]
        state["validator_stderr"] = run.stderr[-4000:]
        if run.returncode != 0 or not validation_report.is_file():
            raise RuntimeError("Independent numerical validation did not complete successfully")
        validation = load_json(validation_report)
        if (
            validation.get("status") != VALIDATION_PASS
            or validation.get("parameter_fitting_authorized") is not False
            or validation.get("production_approved") is not False
            or validation.get("automatic_scientific_classification_applied") is not False
        ):
            raise RuntimeError("Independent numerical validation differs from its exact fail-closed gate")
        state.update(
            {
                "status": "numerical_validation_complete_postqm_dispositions_required",
                "updated_at_utc": now(),
                "independent_validation": artifact(validation_report),
                "production_approved": False,
                "parameter_fitting_authorized": False,
                "automatic_scientific_classification_applied": False,
                "release_boundary": (
                    "Raw curves and interaction energies were independently reconstructed. Every selected "
                    "orientation still requires a signed applicable/weak/unfavourable disposition; parameter "
                    "fitting and all downstream structure/MD stages remain blocked."
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
                "production_approved": False,
                "parameter_fitting_authorized": False,
                "automatic_scientific_classification_applied": False,
            }
        )
        atomic_json(controller_report, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
