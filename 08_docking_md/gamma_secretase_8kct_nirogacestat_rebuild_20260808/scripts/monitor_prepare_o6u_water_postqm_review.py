#!/usr/bin/env python3
"""Wait for independent water-QM reconstruction, then prepare post-QM review."""

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


UPSTREAM_READY = "numerical_validation_complete_postqm_dispositions_required"
PREPARATION_PASS = "postqm_review_template_ready_pending_signed_dispositions"


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
    args = parser.parse_args()
    root = args.root.resolve()
    if args.poll_seconds < 10:
        raise SystemExit("Poll interval must be at least 10 seconds")
    runs = root / "server_runs/o6u_parameterization"
    scripts = root / "scripts"
    upstream = runs / (
        "water_interaction_qm_validation_controller_20260811_v1/"
        "O6U_WATER_INTERACTION_QM_VALIDATION_CONTROLLER.json"
    )
    validation = runs / (
        "water_interaction_qm_independent_validation_formal_mp2_20260811_v1/"
        "O6U_WATER_INTERACTION_QM_INDEPENDENT_VALIDATION.json"
    )
    authorization = runs / (
        "ffparam_water_preqm_authorization_formal_mp2_20260811_v1/"
        "O6U_WATER_PROBE_PREQM_AUTHORIZATION.json"
    )
    review_dir = runs / "water_interaction_postqm_review_formal_mp2_20260811_v1"
    review_report = review_dir / "O6U_WATER_POSTQM_REVIEW_TEMPLATE.json"
    controller_dir = runs / "water_interaction_postqm_review_controller_20260811_v1"
    controller_report = controller_dir / "O6U_WATER_POSTQM_REVIEW_CONTROLLER.json"
    controller_dir.mkdir(parents=True, exist_ok=True)
    lock_handle = (controller_dir / ".controller.lock").open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise SystemExit("Another post-QM review controller owns the lock") from exc
    lock_handle.write(str(os.getpid()) + "\n")
    lock_handle.flush()

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_water_postqm_review_controller",
        "status": "waiting_for_independent_water_qm_numerical_validation",
        "pid": os.getpid(),
        "started_at_utc": now(),
        "updated_at_utc": now(),
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "automatic_scientific_classification_applied": False,
        "upstream_controller": str(upstream),
        "formal_review_directory": str(review_dir),
        "release_boundary": "This controller may prepare a pending review bundle only; it cannot decide dispositions or authorize fitting.",
    }
    atomic_json(controller_report, state)
    while True:
        if upstream.is_file():
            upstream_record = load_json(upstream)
            status = str(upstream_record.get("status", ""))
            state["upstream_status"] = status
            if failed(status):
                state.update({"status": "fail_closed_upstream_numerical_validation", "updated_at_utc": now()})
                atomic_json(controller_report, state)
                return 2
            if status == UPSTREAM_READY:
                break
        state["updated_at_utc"] = now()
        atomic_json(controller_report, state)
        time.sleep(args.poll_seconds)

    try:
        if review_dir.exists():
            raise RuntimeError("Formal post-QM review directory already exists; refusing reuse")
        if not validation.is_file() or not authorization.is_file():
            raise RuntimeError("Formal numerical validation or authorization report is missing")
        state.update({
            "status": "preparing_hash_bound_postqm_review_bundle",
            "updated_at_utc": now(),
            "upstream_controller_record": artifact(upstream),
            "numerical_validation_report": artifact(validation),
            "preqm_authorization_report": artifact(authorization),
        })
        atomic_json(controller_report, state)
        command = [
            sys.executable, str(scripts / "prepare_o6u_water_postqm_review.py"),
            "--validation-report", str(validation),
            "--authorization-report", str(authorization),
            "--output-dir", str(review_dir),
            "--role", "formal_postqm_template",
        ]
        run = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
        state["preparer_returncode"] = run.returncode
        state["preparer_stdout"] = run.stdout[-4000:]
        state["preparer_stderr"] = run.stderr[-4000:]
        if run.returncode != 0 or not review_report.is_file():
            raise RuntimeError("Formal post-QM review preparation failed")
        report = load_json(review_report)
        if (
            report.get("status") != PREPARATION_PASS
            or report.get("template_role") != "formal_postqm_template"
            or report.get("parameter_fitting_authorized") is not False
            or report.get("automatic_scientific_classification_applied") is not False
        ):
            raise RuntimeError("Post-QM review template differs from its exact gate")
        state.update({
            "status": "postqm_review_bundle_ready_pending_signed_dispositions",
            "updated_at_utc": now(),
            "review_template_report": artifact(review_report),
            "production_approved": False,
            "parameter_fitting_authorized": False,
            "automatic_scientific_classification_applied": False,
            "release_boundary": (
                "Every retained orientation must receive a signed APPLICABLE, WEAK, or UNFAVOURABLE disposition. "
                "A separate integrity validation and fitting-target freeze are required."
            ),
        })
        atomic_json(controller_report, state)
        return 0
    except Exception as exc:
        state.update({"status": "fail_closed", "updated_at_utc": now(), "error": str(exc)})
        atomic_json(controller_report, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
