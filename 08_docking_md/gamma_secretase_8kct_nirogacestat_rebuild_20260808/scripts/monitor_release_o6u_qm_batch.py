#!/usr/bin/env python3
"""Wait for the frozen O6U canary, validate it, then release the QM batch."""

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


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    os.replace(temporary, path)


def pid_is_live(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def checked_file(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SystemExit(f"Missing or empty required input: {resolved}")
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canary-record", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--selection-report", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--batch-runner", required=True, type=Path)
    parser.add_argument("--batch-output-dir", required=True, type=Path)
    parser.add_argument("--controller-report", required=True, type=Path)
    parser.add_argument("--interval-seconds", type=int, default=300)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    args = parser.parse_args()

    if not 30 <= args.interval_seconds <= 1800:
        raise SystemExit("Polling interval must be between 30 and 1800 seconds")
    if not 1 <= args.threads <= 24 or not 16 <= args.memory_gib <= 128:
        raise SystemExit("Resources fall outside the frozen project limits")

    canary_record = checked_file(args.canary_record)
    source_sdf = checked_file(args.source_sdf)
    validator = checked_file(args.validator)
    selection_report = checked_file(args.selection_report)
    runner = checked_file(args.runner)
    batch_runner = checked_file(args.batch_runner)
    report_path = args.controller_report.resolve()
    batch_output = args.batch_output_dir.resolve()
    lock_path = report_path.with_suffix(report_path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    lock_handle = lock_path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another O6U QM release controller already holds the lock")

    inputs = {
        "canary_record": {"path": str(canary_record), "sha256_at_start": sha256(canary_record)},
        "source_sdf": {"path": str(source_sdf), "sha256": sha256(source_sdf)},
        "validator": {"path": str(validator), "sha256": sha256(validator)},
        "selection_report": {"path": str(selection_report), "sha256": sha256(selection_report)},
        "runner": {"path": str(runner), "sha256": sha256(runner)},
        "batch_runner": {"path": str(batch_runner), "sha256": sha256(batch_runner)},
    }
    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_qm_fail_closed_release_controller",
        "status": "waiting_for_canary",
        "production_approved": False,
        "created_at_utc": now(),
        "updated_at_utc": now(),
        "inputs": inputs,
        "batch_output_dir": str(batch_output),
        "interval_seconds": args.interval_seconds,
        "release_boundary": (
            "Only the four frozen non-canary QM representatives may be released after "
            "an independent canary PASS; force-field fitting and MD remain blocked."
        ),
    }
    write_json(report_path, state)

    while True:
        record = read_json(canary_record)
        status = record.get("status")
        state.update(
            {
                "updated_at_utc": now(),
                "observed_canary_status": status,
                "observed_canary_record_sha256": sha256(canary_record),
            }
        )
        if status == "running":
            pid = int(record.get("pid", -1))
            state["observed_canary_pid"] = pid
            state["observed_canary_pid_live"] = pid_is_live(pid)
            if not pid_is_live(pid):
                time.sleep(15)
                record = read_json(canary_record)
                if record.get("status") == "running":
                    state.update({"status": "fail_closed_incomplete_canary_record", "failed_at_utc": now()})
                    write_json(report_path, state)
                    return 1
                continue
            write_json(report_path, state)
            time.sleep(args.interval_seconds)
            continue

        if status != "pass_optimization_canary":
            state.update({"status": "fail_closed_canary_not_passed", "failed_at_utc": now()})
            write_json(report_path, state)
            return 1

        validation_path = report_path.parent / "O6U_MP2_CANARY_INDEPENDENT_VALIDATION.json"
        validation_command = [
            sys.executable,
            str(validator),
            "--record",
            str(canary_record),
            "--source-sdf",
            str(source_sdf),
            "--report",
            str(validation_path),
            "--expected-role",
            "canary",
        ]
        state.update({"status": "validating_canary", "validation_command": validation_command})
        write_json(report_path, state)
        validation_run = subprocess.run(validation_command, check=False)
        if validation_run.returncode != 0 or not validation_path.is_file():
            state.update(
                {
                    "status": "fail_closed_independent_canary_validation",
                    "validation_returncode": validation_run.returncode,
                    "failed_at_utc": now(),
                }
            )
            write_json(report_path, state)
            return 1
        validation = read_json(validation_path)
        if (
            validation.get("status") != "pass_canary_independently_reconstructed"
            or validation.get("schema_version") != "1.1"
            or validation.get("optimizer_strategy_version") != "cartesian_rfo_trust020_v1"
            or validation.get("production_approved") is not False
        ):
            state.update({"status": "fail_closed_invalid_validation_record", "failed_at_utc": now()})
            write_json(report_path, state)
            return 1

        if batch_output.exists() and any(batch_output.iterdir()):
            state.update({"status": "fail_closed_nonempty_batch_output", "failed_at_utc": now()})
            write_json(report_path, state)
            return 1

        batch_command = [
            sys.executable,
            str(batch_runner),
            "--selection-report",
            str(selection_report),
            "--canary-validation",
            str(validation_path),
            "--source-sdf",
            str(source_sdf),
            "--runner",
            str(runner),
            "--validator",
            str(validator),
            "--output-dir",
            str(batch_output),
            "--threads",
            str(args.threads),
            "--memory-gib",
            str(args.memory_gib),
            "--scratch-root",
            str(args.scratch_root.resolve()),
        ]
        state.update(
            {
                "status": "running_representative_batch",
                "canary_validation": {"path": str(validation_path), "sha256": sha256(validation_path)},
                "batch_command": batch_command,
                "batch_started_at_utc": now(),
            }
        )
        write_json(report_path, state)
        batch_run = subprocess.run(batch_command, check=False)
        batch_state_path = batch_output / "O6U_MP2_REPRESENTATIVE_BATCH.json"
        if batch_run.returncode != 0 or not batch_state_path.is_file():
            state.update(
                {
                    "status": "fail_closed_representative_batch",
                    "batch_returncode": batch_run.returncode,
                    "failed_at_utc": now(),
                }
            )
            write_json(report_path, state)
            return 1
        batch_state = read_json(batch_state_path)
        if batch_state.get("status") != "pass_all_representatives":
            state.update({"status": "fail_closed_invalid_batch_state", "failed_at_utc": now()})
            write_json(report_path, state)
            return 1

        state.update(
            {
                "status": "pass_qm_representative_batch",
                "completed_at_utc": now(),
                "batch_returncode": batch_run.returncode,
                "batch_state": {"path": str(batch_state_path), "sha256": sha256(batch_state_path)},
            }
        )
        write_json(report_path, state)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
