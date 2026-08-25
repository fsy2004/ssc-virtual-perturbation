#!/usr/bin/env python3
"""Wait for the O6U QM batch, then independently freeze the five-member ensemble."""

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


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checked(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SystemExit(f"Missing or empty required input: {resolved}")
    return resolved


def load(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SystemExit(f"Expected JSON object: {path}")
    return value


def write(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-controller", required=True, type=Path)
    parser.add_argument("--selection-report", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--canary-record", required=True, type=Path)
    parser.add_argument("--canary-validation", required=True, type=Path)
    parser.add_argument("--batch-state", required=True, type=Path)
    parser.add_argument("--aggregator", required=True, type=Path)
    parser.add_argument("--ensemble-report", required=True, type=Path)
    parser.add_argument("--controller-report", required=True, type=Path)
    parser.add_argument("--interval-seconds", type=int, default=300)
    args = parser.parse_args()

    if not 30 <= args.interval_seconds <= 1800:
        raise SystemExit("Polling interval must be 30-1800 seconds")
    upstream = checked(args.upstream_controller)
    selection = checked(args.selection_report)
    source = checked(args.source_sdf)
    canary = checked(args.canary_record)
    aggregator = checked(args.aggregator)
    canary_validation = args.canary_validation.resolve()
    batch_state = args.batch_state.resolve()
    ensemble_report = args.ensemble_report.resolve()
    controller_report = args.controller_report.resolve()
    if ensemble_report.exists():
        raise SystemExit(f"Refusing pre-existing ensemble report: {ensemble_report}")

    lock_path = controller_report.with_suffix(controller_report.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_handle = lock_path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another MP2 ensemble-freeze controller holds the lock")

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_ensemble_freeze_controller",
        "status": "waiting_for_qm_batch",
        "production_approved": False,
        "created_at_utc": now(),
        "updated_at_utc": now(),
        "interval_seconds": args.interval_seconds,
        "inputs": {
            "upstream_controller": {"path": str(upstream), "sha256_at_start": sha256(upstream)},
            "selection_report": {"path": str(selection), "sha256": sha256(selection)},
            "source_sdf": {"path": str(source), "sha256": sha256(source)},
            "canary_record": {"path": str(canary), "sha256_at_start": sha256(canary)},
            "aggregator": {"path": str(aggregator), "sha256": sha256(aggregator)},
        },
        "future_inputs": {
            "canary_validation": str(canary_validation),
            "batch_state": str(batch_state),
        },
        "ensemble_report": str(ensemble_report),
        "release_boundary": "Geometry ensemble only; force-field fitting, CHARMM-GUI, and MD remain blocked.",
    }
    write(controller_report, state)

    terminal_upstream_failures = {
        "fail_closed_incomplete_canary_record",
        "fail_closed_canary_not_passed",
        "fail_closed_independent_canary_validation",
        "fail_closed_invalid_validation_record",
        "fail_closed_nonempty_batch_output",
        "fail_closed_representative_batch",
        "fail_closed_invalid_batch_state",
    }
    while True:
        upstream_state = load(upstream)
        upstream_status = str(upstream_state.get("status", ""))
        state.update(
            {
                "updated_at_utc": now(),
                "observed_upstream_status": upstream_status,
                "observed_upstream_sha256": sha256(upstream),
            }
        )
        if upstream_status in terminal_upstream_failures:
            state.update({"status": "fail_closed_upstream_qm", "failed_at_utc": now()})
            write(controller_report, state)
            return 1
        if upstream_status != "pass_qm_representative_batch":
            write(controller_report, state)
            time.sleep(args.interval_seconds)
            continue

        checked(canary_validation)
        checked(batch_state)
        upstream_canary_validation = upstream_state.get("canary_validation", {})
        upstream_batch_state = upstream_state.get("batch_state", {})
        if upstream_canary_validation.get("path") != str(canary_validation) or upstream_canary_validation.get("sha256") != sha256(canary_validation):
            state.update({"status": "fail_closed_canary_validation_binding", "failed_at_utc": now()})
            write(controller_report, state)
            return 1
        if upstream_batch_state.get("path") != str(batch_state) or upstream_batch_state.get("sha256") != sha256(batch_state):
            state.update({"status": "fail_closed_batch_binding", "failed_at_utc": now()})
            write(controller_report, state)
            return 1
        if ensemble_report.exists():
            state.update({"status": "fail_closed_preexisting_ensemble_report", "failed_at_utc": now()})
            write(controller_report, state)
            return 1

        command = [
            sys.executable,
            str(aggregator),
            "--selection-report", str(selection),
            "--source-sdf", str(source),
            "--canary-record", str(canary),
            "--canary-validation", str(canary_validation),
            "--batch-state", str(batch_state),
            "--report", str(ensemble_report),
        ]
        state.update({"status": "running_independent_ensemble_freeze", "command": command, "started_at_utc": now()})
        write(controller_report, state)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0 or not ensemble_report.is_file():
            state.update({"status": "fail_closed_ensemble_freeze", "returncode": completed.returncode, "failed_at_utc": now()})
            write(controller_report, state)
            return 1
        ensemble = load(ensemble_report)
        if ensemble.get("status") != "pass_five_member_ensemble_independently_reconstructed" or ensemble.get("production_approved") is not False:
            state.update({"status": "fail_closed_invalid_ensemble_report", "failed_at_utc": now()})
            write(controller_report, state)
            return 1
        state.update(
            {
                "status": "pass_mp2_ensemble_frozen",
                "completed_at_utc": now(),
                "returncode": completed.returncode,
                "ensemble": {"path": str(ensemble_report), "sha256": sha256(ensemble_report)},
            }
        )
        write(controller_report, state)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
