#!/usr/bin/env python3
"""Run the four non-canary O6U QM representatives sequentially after release."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-report", required=True, type=Path)
    parser.add_argument("--canary-validation", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    args = parser.parse_args()

    required = [args.selection_report, args.canary_validation, args.source_sdf, args.runner, args.validator]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    if not 1 <= args.threads <= 24 or not 16 <= args.memory_gib <= 128:
        raise SystemExit("Batch resources fall outside the frozen project limits")

    selection = json.loads(args.selection_report.read_text(encoding="utf-8"))
    validation = json.loads(args.canary_validation.read_text(encoding="utf-8"))
    if selection.get("status") != "pass" or selection.get("selected_count") != 5:
        raise SystemExit("Frozen representative selection report is not a five-structure pass")
    if selection.get("method", {}).get("tfd_cutoff") != 0.2:
        raise SystemExit("Frozen TFD cutoff differs from 0.20")
    if validation.get("status") != "pass_canary_independently_reconstructed":
        raise SystemExit("Independent canary validation has not released the batch")
    if validation.get("schema_version") != "1.1":
        raise SystemExit("Independent canary validation is not the repaired v1.1 gate")
    if validation.get("optimizer_strategy_version") != "cartesian_rfo_trust020_v1":
        raise SystemExit("Independent canary validation did not approve the frozen Cartesian/RFO strategy")
    if validation.get("production_approved") is not False:
        raise SystemExit("Canary validation improperly claims production approval")

    global_minimum = int(selection["crest_global_minimum_frame_1based"])
    targets = [item for item in selection["selected"] if int(item["crest_frame_1based"]) != global_minimum]
    if len(targets) != 4:
        raise SystemExit(f"Expected four non-canary representatives, found {len(targets)}")
    if len({int(item["crest_frame_1based"]) for item in targets}) != 4:
        raise SystemExit("Non-canary representative frame identities are not unique")

    output_dir = args.output_dir.resolve()
    state_path = output_dir / "O6U_MP2_REPRESENTATIVE_BATCH.json"
    if state_path.exists():
        previous = json.loads(state_path.read_text(encoding="utf-8"))
        if previous.get("status") == "pass_all_representatives":
            print(json.dumps({"status": previous["status"], "resume": "already_complete"}, sort_keys=True))
            return 0
        raise SystemExit("Refusing to overwrite or silently resume an incomplete prior batch")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"Refusing non-empty batch output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_representative_batch",
        "status": "running",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_report": {"path": str(args.selection_report.resolve()), "sha256": sha256(args.selection_report)},
        "canary_validation": {"path": str(args.canary_validation.resolve()), "sha256": sha256(args.canary_validation)},
        "source_sdf": {"path": str(args.source_sdf.resolve()), "sha256": sha256(args.source_sdf)},
        "runner": {"path": str(args.runner.resolve()), "sha256": sha256(args.runner)},
        "validator": {"path": str(args.validator.resolve()), "sha256": sha256(args.validator)},
        "threads": args.threads,
        "memory_gib": args.memory_gib,
        "execution": "strictly sequential",
        "targets": [],
        "release_boundary": "QM targets only; ligand parameters and production MD remain blocked.",
    }
    write_json(state_path, state)

    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": str(args.threads),
            "MKL_NUM_THREADS": str(args.threads),
            "OPENBLAS_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    for item in targets:
        frame = int(item["crest_frame_1based"])
        start = Path(item["output_xyz"])
        if not start.is_file() or sha256(start) != item["output_xyz_sha256"]:
            state["status"] = "fail_input_integrity"
            state["failed_frame_1based"] = frame
            write_json(state_path, state)
            raise SystemExit(f"Representative frame {frame} input hash mismatch")
        target_dir = output_dir / f"frame_{frame:04d}"
        command = [
            sys.executable,
            str(args.runner.resolve()),
            "--source-sdf",
            str(args.source_sdf.resolve()),
            "--start-xyz",
            str(start.resolve()),
            "--output-dir",
            str(target_dir),
            "--threads",
            str(args.threads),
            "--memory-gib",
            str(args.memory_gib),
            "--scratch-root",
            str(args.scratch_root.resolve()),
            "--role",
            "representative_target",
        ]
        target_state = {
            "crest_frame_1based": frame,
            "start_xyz": str(start.resolve()),
            "start_xyz_sha256": sha256(start),
            "status": "running",
            "command": command,
        }
        state["targets"].append(target_state)
        write_json(state_path, state)
        completed = subprocess.run(command, env=environment, check=False)
        record_path = target_dir / "o6u_mp2_631gd_optimization_representative_target.json"
        if completed.returncode != 0 or not record_path.is_file():
            target_state.update({"status": "fail", "returncode": completed.returncode})
            state.update({"status": "fail_representative", "failed_frame_1based": frame})
            write_json(state_path, state)
            return 1
        record = json.loads(record_path.read_text(encoding="utf-8"))
        if record.get("status") != "pass_optimization_representative_target":
            target_state.update({"status": "fail_record", "returncode": completed.returncode})
            state.update({"status": "fail_representative", "failed_frame_1based": frame})
            write_json(state_path, state)
            return 1
        validation_path = target_dir / "O6U_MP2_REPRESENTATIVE_TARGET_VALIDATION.json"
        validation_command = [
            sys.executable,
            str(args.validator.resolve()),
            "--record",
            str(record_path.resolve()),
            "--source-sdf",
            str(args.source_sdf.resolve()),
            "--report",
            str(validation_path.resolve()),
            "--expected-role",
            "representative_target",
        ]
        validation_completed = subprocess.run(validation_command, env=environment, check=False)
        if validation_completed.returncode != 0 or not validation_path.is_file():
            target_state.update(
                {
                    "status": "fail_independent_validation",
                    "validation_returncode": validation_completed.returncode,
                    "validation_command": validation_command,
                }
            )
            state.update({"status": "fail_representative_validation", "failed_frame_1based": frame})
            write_json(state_path, state)
            return 1
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if (
            validation.get("status") != "pass_representative_independently_reconstructed"
            or validation.get("production_approved") is not False
            or validation.get("role") != "representative_target"
        ):
            target_state.update({"status": "fail_independent_validation_record"})
            state.update({"status": "fail_representative_validation", "failed_frame_1based": frame})
            write_json(state_path, state)
            return 1
        target_state.update(
            {
                "status": "pass",
                "returncode": completed.returncode,
                "record": str(record_path),
                "record_sha256": sha256(record_path),
                "validation": str(validation_path),
                "validation_sha256": sha256(validation_path),
            }
        )
        write_json(state_path, state)

    state["status"] = "pass_all_representatives"
    state["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(state_path, state)
    print(json.dumps({"status": state["status"], "state_sha256": sha256(state_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
