#!/usr/bin/env python3
"""Wait for the frozen O6U pre-QM authorization, then run raw water QM."""

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


UPSTREAM_READY = "preqm_orientation_authorization_frozen"
RUNNER_PASS = "pass_raw_water_qm_execution_outputs_present_validation_required"


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


def failed(status: str) -> bool:
    lowered = status.lower()
    return any(token in lowered for token in ("fail", "error", "no_go", "blocked"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.poll_seconds < 10 or not 1 <= args.workers <= 20:
        raise SystemExit("Poll interval must be >=10 seconds and workers between 1 and 20")

    runs = root / "server_runs/o6u_parameterization"
    scripts = root / "scripts"
    preqm_controller = runs / "ffparam_water_preqm_release_controller_20260811_v1/O6U_WATER_PROBE_PREQM_RELEASE_CONTROLLER.json"
    authorization = runs / "ffparam_water_preqm_authorization_formal_mp2_20260811_v1/O6U_WATER_PROBE_PREQM_AUTHORIZATION.json"
    generation = runs / "ffparam_water_input_generation_formal_mp2_20260811_v1/O6U_FFPARAM_WATER_INPUT_GENERATION.json"
    batch_dir = runs / "water_interaction_qm_batch_formal_mp2_20260811_v1"
    batch_report = batch_dir / "O6U_WATER_INTERACTION_QM_BATCH.json"
    controller_dir = runs / "water_interaction_qm_execution_controller_20260811_v1"
    controller_report = controller_dir / "O6U_WATER_INTERACTION_QM_EXECUTION_CONTROLLER.json"
    controller_dir.mkdir(parents=True, exist_ok=True)

    lock_handle = (controller_dir / ".controller.lock").open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("Another water-QM execution controller owns the lock") from exc
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_water_interaction_qm_execution_controller",
        "status": "waiting_for_frozen_preqm_authorization",
        "pid": os.getpid(),
        "started_at_utc": now(),
        "updated_at_utc": now(),
        "workers": args.workers,
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "preqm_controller": str(preqm_controller),
        "formal_generation_report": str(generation),
        "formal_authorization_report": str(authorization),
        "formal_batch_directory": str(batch_dir),
        "release_boundary": (
            "This controller may run only the orientation IDs authorized by the frozen 70-row record. It cannot "
            "classify post-QM dispositions, fit parameters, construct a CHARMM-GUI system, or run MD."
        ),
    }
    atomic_json(controller_report, state)

    while True:
        if preqm_controller.is_file():
            upstream = load_json(preqm_controller)
            status = str(upstream.get("status", ""))
            state["upstream_status"] = status
            if failed(status):
                state.update({"status": "fail_closed_upstream_preqm_controller", "updated_at_utc": now()})
                atomic_json(controller_report, state)
                return 2
            if status == UPSTREAM_READY:
                if upstream.get("water_interaction_qm_authorized") is not True:
                    state.update({"status": "fail_closed_upstream_ready_without_qm_authorization", "updated_at_utc": now()})
                    atomic_json(controller_report, state)
                    return 2
                break
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    try:
        if not authorization.is_file() or not generation.is_file():
            raise RuntimeError("Formal generation or authorization report is missing")
        if batch_dir.exists():
            raise RuntimeError("Formal water-QM batch directory already exists; refusing reuse")
        state.update(
            {
                "status": "executing_frozen_authorized_water_qm_batch",
                "updated_at_utc": now(),
                "preqm_controller_record": artifact(preqm_controller),
                "authorization_report": artifact(authorization),
                "generation_report": artifact(generation),
            }
        )
        atomic_json(controller_report, state)
        command = [
            sys.executable,
            str(scripts / "run_o6u_water_interaction_qm.py"),
            "--generation-report", str(generation),
            "--authorization-report", str(authorization),
            "--output-dir", str(batch_dir),
            "--role", "formal_execution",
            "--workers", str(args.workers),
        ]
        run = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
        state["runner_returncode"] = run.returncode
        state["runner_stdout"] = run.stdout[-4000:]
        state["runner_stderr"] = run.stderr[-4000:]
        if not batch_report.is_file():
            raise RuntimeError("Water-QM runner did not create its batch report")
        batch = load_json(batch_report)
        state["batch_report"] = artifact(batch_report)
        if run.returncode != 0 or batch.get("status") != RUNNER_PASS:
            state.update(
                {
                    "status": "fail_closed_raw_water_qm_batch",
                    "updated_at_utc": now(),
                    "batch_status": batch.get("status"),
                }
            )
            atomic_json(controller_report, state)
            return 2
        state.update(
            {
                "status": "raw_water_qm_batch_complete_independent_validation_required",
                "updated_at_utc": now(),
                "batch_status": batch.get("status"),
                "production_approved": False,
                "parameter_fitting_authorized": False,
                "release_boundary": (
                    "Authorized raw HF/6-31G(d) outputs are complete. Independent numerical reconstruction and "
                    "post-QM dispositions are required before any parameter fitting."
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
            }
        )
        atomic_json(controller_report, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
