#!/usr/bin/env python3
"""Fail-closed integrity audit for the three GROMACS realizations."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from md_contract import (
    REALIZATION_IDS,
    SYSTEM_ID,
    artifact_path,
    build_contract_sha256,
    canonical_json_sha256,
    release_manifest_contract_sha256,
    report_payload_sha256,
    validate_common_minimization_record,
    validate_bound_gmx_identity,
)


NEVER_ALLOW = (
    re.compile(r"fatal error", re.IGNORECASE),
    re.compile(r"segmentation fault", re.IGNORECASE),
    re.compile(r"\b(?:nan|inf)\b", re.IGNORECASE),
    re.compile(r"lincs warning", re.IGNORECASE),
    re.compile(r"settle.*(?:error|warning)", re.IGNORECASE),
    re.compile(r"constraint.*(?:error|warning)", re.IGNORECASE),
    re.compile(r"pressure scaling more than", re.IGNORECASE),
    re.compile(r"domain decomposition.*(?:error|failed|fatal)", re.IGNORECASE),
)
GROMPP_WARNING = re.compile(r"^\s*WARNING\s+[0-9]+", re.IGNORECASE | re.MULTILINE)
LAST_FRAME = re.compile(r"Last frame\s+\d+\s+time\s+([0-9.+\-Ee]+)", re.IGNORECASE)
LAST_ENERGY_FRAME = re.compile(r"Last energy frame read\s+\d+\s+time\s+([0-9.+\-Ee]+)", re.IGNORECASE)
LOG_STEP_TIME = re.compile(
    r"^\s*Step\s+Time\s*$\s*^\s*([0-9]+)\s+([0-9.+\-Ee]+)\s*$",
    re.IGNORECASE | re.MULTILINE,
)
CHECKPOINT_STEP = re.compile(r"^\s*step\s*=\s*([0-9]+)", re.IGNORECASE | re.MULTILINE)
CHECKPOINT_TIME = re.compile(r"^\s*t\s*=\s*([0-9.+\-Ee]+)", re.IGNORECASE | re.MULTILINE)
GRO_TIME = re.compile(r"(?:^|\s)t\s*=\s*([0-9.+\-Ee]+)(?:\s|$)", re.IGNORECASE)
TRAJECTORY_ATOMS = re.compile(r"^\s*#\s*Atoms\s+([0-9]+)\s*$", re.IGNORECASE | re.MULTILINE)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def projected_storage_bytes(canary_bytes: int, production_ns: float, canary_ns: float) -> tuple[int, float]:
    if canary_bytes <= 0:
        raise ValueError("Canary output byte count must be positive")
    if not all(math.isfinite(value) and value > 0.0 for value in (production_ns, canary_ns)):
        raise ValueError("Production and canary durations must be positive and finite")
    multiplier = production_ns / canary_ns
    return math.ceil(canary_bytes * multiplier), multiplier


def inside(root: Path, value: str) -> Path:
    path = (root / value).resolve()
    if path != root and root not in path.parents:
        raise ValueError(f"Path escapes package root: {value}")
    return path


def run_check(argv: list[str], cwd: Path) -> dict[str, Any]:
    result = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, errors="replace", check=False)
    return {
        "argv": argv,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def last_float(pattern: re.Pattern[str], text: str) -> float | None:
    matches = pattern.findall(text)
    return float(matches[-1]) if matches else None


def last_log_step_time(path: Path) -> tuple[int | None, float | None]:
    matches = LOG_STEP_TIME.findall(path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        return None, None
    step, time_ps = matches[-1]
    return int(step), float(time_ps)


def gro_time_ps(path: Path) -> float | None:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        title = handle.readline()
    return last_float(GRO_TIME, title)


def gro_atom_count_and_box(path: Path) -> tuple[int, list[float]]:
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    if len(lines) < 3:
        raise ValueError(f"GRO file is truncated: {path}")
    atom_count = int(lines[1].strip())
    if atom_count <= 0 or len(lines) != atom_count + 3:
        raise ValueError(f"GRO atom count/line count mismatch: {path}")
    box = [float(value) for value in lines[-1].split()]
    if len(box) not in {3, 9} or not all(math.isfinite(value) for value in box):
        raise ValueError(f"GRO box vector is invalid: {path}")
    return atom_count, box


def audit_append_only_command_chain(run_dir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    errors: list[str] = []
    records: list[dict[str, Any]] = []
    log_dir = run_dir / "command_logs"
    finished_paths = sorted(log_dir.glob("[0-9][0-9][0-9][0-9]_*.finished.json")) if log_dir.is_dir() else []
    previous_finished_hash: str | None = None
    previous_sequence = 0
    for finished_path in finished_paths:
        try:
            finished = load_json(finished_path)
            sequence = int(finished.get("sequence"))
            if sequence <= previous_sequence:
                errors.append("append-only command sequence is not strictly increasing")
            started_path = Path(str(finished.get("started_record", ""))).resolve()
            if started_path.parent != log_dir.resolve() or not started_path.is_file():
                errors.append(f"finished record points to missing/outside started record: {finished_path}")
                continue
            if finished.get("started_record_sha256") != sha256(started_path):
                errors.append(f"started-record hash mismatch: {finished_path}")
            started = load_json(started_path)
            if started.get("sequence") != sequence:
                errors.append(f"started/finished sequence mismatch: {finished_path}")
            if started.get("previous_finished_record_sha256") != previous_finished_hash:
                errors.append(f"append-only command hash chain mismatch: {finished_path}")
            for stream in ("stdout", "stderr"):
                stream_artifact = finished.get(stream)
                if not isinstance(stream_artifact, dict):
                    errors.append(f"{finished_path} lacks immutable {stream} provenance")
                    continue
                stream_path = Path(str(stream_artifact.get("path", ""))).resolve()
                if stream_path.parent != log_dir.resolve() or not stream_path.is_file() or stream_artifact.get("sha256") != sha256(stream_path):
                    errors.append(f"{finished_path} {stream} file is missing, outside the log directory, or changed")
            argv = started.get("argv")
            if not isinstance(argv, list) or not all(isinstance(item, str) for item in argv):
                errors.append(f"invalid argv in {started_path}")
            elif "-maxwarn" in argv:
                errors.append(f"-maxwarn is prohibited: {started_path}")
            records.append({
                "sequence": sequence,
                "started": {"path": str(started_path), "sha256": sha256(started_path)},
                "finished": {"path": str(finished_path), "sha256": sha256(finished_path)},
                "argv": argv,
                "returncode": finished.get("returncode"),
                "runtime_outputs": finished.get("runtime_outputs"),
            })
            if finished.get("returncode") != 0:
                errors.append(f"retained command failed: {finished_path}")
            previous_finished_hash = sha256(finished_path)
            previous_sequence = sequence
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"cannot validate append-only command record {finished_path}: {exc}")
    return records, errors


def audit_production_tpr_record(manifest: dict[str, Any], run_dir: Path, rid: str, tpr: Path) -> tuple[dict[str, Any] | None, list[str]]:
    errors: list[str] = []
    record_path = run_dir / "production_tpr_record.json"
    if not record_path.is_file():
        return None, [f"missing production TPR provenance record: {record_path}"]
    try:
        record = load_json(record_path)
        expected = {
            "system_id": SYSTEM_ID,
            "construction_id": "build01",
            "realization_id": rid,
            "build_contract_sha256": build_contract_sha256(manifest),
            "production_tpr_sha256": sha256(tpr),
        }
        for key, value in expected.items():
            if record.get(key) != value:
                errors.append(f"{rid} production TPR record mismatch at {key}")
        grompp_record = Path(str(record.get("grompp_finished_record", ""))).resolve()
        if grompp_record.parent != (run_dir / "command_logs").resolve() or not grompp_record.is_file() or record.get("grompp_finished_record_sha256") != sha256(grompp_record):
            errors.append(f"{rid} production TPR grompp provenance is missing or changed")
        return record, errors
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return None, [f"cannot validate {rid} production TPR record: {exc}"]


def audit_equilibration_one(
    package_root: Path,
    manifest: dict[str, Any],
    realization: dict[str, Any],
    construction: dict[str, Any],
    common_minimization_sha256: str,
    gmx: str,
    use_gmx: bool,
    strict: bool,
) -> dict[str, Any]:
    rid = str(realization.get("id"))
    run_dir = inside(package_root, str(realization.get("run_directory", "")))
    work = run_dir / "work"
    errors: list[str] = []
    warnings: list[str] = []
    stages: list[dict[str, Any]] = []
    expected_hashes = construction.get("equilibration_mdp_sha256", {})
    for stage_index in range(1, 7):
        stem = f"step6.{stage_index}_equilibration"
        mdp_name = "step6.1_equilibration.frozen.mdp" if stage_index == 1 else f"{stem}.mdp"
        required = {label: work / f"{stem}.{suffix}" for label, suffix in (
            ("tpr", "tpr"), ("energy", "edr"), ("log", "log"),
            ("final_coordinates", "gro"), ("checkpoint", "cpt"),
        )}
        required["mdp"] = work / mdp_name
        stage_errors: list[str] = []
        artifacts: dict[str, Any] = {}
        for label, path in required.items():
            if not path.is_file() or path.stat().st_size == 0:
                stage_errors.append(f"missing or empty {label}: {path}")
            else:
                artifacts[label] = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        source_mdp = work / f"{stem}.mdp"
        expected_source_hash = expected_hashes.get(source_mdp.name)
        if not source_mdp.is_file() or sha256(source_mdp) != expected_source_hash:
            stage_errors.append(f"{rid} {source_mdp.name} differs from frozen CHARMM-GUI source hash")
        if required["log"].is_file():
            text = required["log"].read_text(encoding="utf-8", errors="replace")
            if "Finished mdrun" not in text:
                stage_errors.append(f"{rid} {stem} lacks the Finished mdrun marker")
            for pattern in NEVER_ALLOW:
                if pattern.search(text):
                    stage_errors.append(f"{rid} {stem} contains forbidden runtime pattern {pattern.pattern!r}")
        checks: list[dict[str, Any]] = []
        if use_gmx and all(path.is_file() for key, path in required.items() if key != "mdp"):
            checks = [
                run_check([gmx, "check", "-e", required["energy"].name], work),
                run_check([gmx, "check", "-c", required["final_coordinates"].name], work),
            ]
            if any(check["returncode"] != 0 for check in checks):
                stage_errors.append(f"gmx check failed for {rid} {stem}")
        elif strict:
            stage_errors.append(f"strict equilibration audit requires GROMACS checks for {rid} {stem}")
        stages.append({
            "stage": stem, "status": "pass" if not stage_errors else "fail",
            "errors": stage_errors, "artifacts": artifacts, "gmx_checks": checks,
        })
        errors.extend(stage_errors)

    staged_source = run_dir / "staged_source.json"
    if not staged_source.is_file():
        errors.append(f"missing staged-source provenance: {staged_source}")
    else:
        record = load_json(staged_source)
        if record.get("construction_id") != "build01" or record.get("common_minimization_sha256") != common_minimization_sha256:
            errors.append(f"{rid} equilibration does not derive from build01/common minimization")
    seed_record_path = run_dir / "velocity_seed_record.json"
    frozen_first_nvt = work / "step6.1_equilibration.frozen.mdp"
    if not seed_record_path.is_file() or not frozen_first_nvt.is_file():
        errors.append(f"{rid} lacks frozen first-NVT seed provenance")
    else:
        try:
            seed_record = load_json(seed_record_path)
            expected_seed = int(realization.get("velocity_seed"))
            if seed_record.get("velocity_seed") != expected_seed:
                errors.append(f"{rid} velocity-seed record differs from manifest")
            if seed_record.get("source_sha256") != expected_hashes.get("step6.1_equilibration.mdp"):
                errors.append(f"{rid} velocity-seed source hash differs from the frozen first-NVT MDP")
            if seed_record.get("derived_sha256") != sha256(frozen_first_nvt):
                errors.append(f"{rid} frozen first-NVT MDP hash differs from its seed record")
            match = re.search(r"^\s*gen[_-]seed\s*=\s*([0-9]+)", frozen_first_nvt.read_text(encoding="utf-8", errors="replace"), re.IGNORECASE | re.MULTILINE)
            if match is None or int(match.group(1)) != expected_seed:
                errors.append(f"{rid} frozen first-NVT MDP does not contain the manifest seed")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"cannot validate {rid} first-NVT seed provenance: {exc}")
    command_records, command_errors = audit_append_only_command_chain(run_dir)
    errors.extend(command_errors)
    grompp_argvs = [item["argv"] for item in command_records if isinstance(item.get("argv"), list) and len(item["argv"]) > 1 and item["argv"][1] == "grompp"]
    first_matching = [argv for argv in grompp_argvs if "step6.1_equilibration.frozen.mdp" in argv]
    if len(first_matching) != 1 or "-t" in first_matching[0]:
        errors.append(f"{rid} first NVT must be generated exactly once without a checkpoint")
    elif "-c" not in first_matching[0] or first_matching[0][first_matching[0].index("-c") + 1] != "step6.0_minimization.gro":
        errors.append(f"{rid} first NVT did not consume the frozen common-minimization coordinates")
    for stage_index in range(2, 7):
        mdp_name = f"step6.{stage_index}_equilibration.mdp"
        matching = [argv for argv in grompp_argvs if mdp_name in argv]
        expected_cpt = f"step6.{stage_index - 1}_equilibration.cpt"
        if len(matching) != 1 or "-t" not in matching[0] or matching[0][matching[0].index("-t") + 1] != expected_cpt:
            errors.append(f"{rid} {mdp_name} did not consume exactly {expected_cpt}")
    command_log_dir = run_dir / "command_logs"
    command_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace") for path in sorted(command_log_dir.glob("*.log"))
    ) if command_log_dir.is_dir() else ""
    warning_count = len(GROMPP_WARNING.findall(command_text))
    if warning_count:
        errors.append(f"{rid} equilibration grompp emitted {warning_count} numbered warning(s)")
    for pattern in NEVER_ALLOW:
        if pattern.search(command_text):
            errors.append(f"{rid} equilibration command logs contain forbidden runtime pattern {pattern.pattern!r}")
    return {
        "system_id": SYSTEM_ID, "construction_id": "build01", "realization_id": rid,
        "phase": "equilibration", "status": "pass" if not errors else "fail",
        "errors": errors, "warnings": warnings, "stages": stages,
        "append_only_command_records": command_records,
    }


def audit_one(
    package_root: Path,
    manifest: dict[str, Any],
    realization: dict[str, Any],
    phase: str,
    expected_ns: float,
    common_minimization_sha256: str,
    gmx: str,
    use_gmx: bool,
    strict: bool,
) -> dict[str, Any]:
    rid = str(realization.get("id"))
    run_dir = inside(package_root, str(realization.get("run_directory", "")))
    work = run_dir / "work"
    stem = "production"
    required = {
        "tpr": work / f"{stem}.tpr",
        "trajectory": work / f"{stem}.xtc",
        "energy": work / f"{stem}.edr",
        "log": work / f"{stem}.log",
        "final_coordinates": work / f"{stem}.gro",
        "checkpoint": work / f"{stem}.cpt",
    }
    errors: list[str] = []
    warnings: list[str] = []
    artifacts: dict[str, Any] = {}
    for label, path in required.items():
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty {label}: {path}")
        else:
            artifacts[label] = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}

    seed_record_path = run_dir / "velocity_seed_record.json"
    staged_source_path = run_dir / "staged_source.json"
    frozen_first_nvt = work / "step6.1_equilibration.frozen.mdp"
    for label, path in (
        ("velocity seed record", seed_record_path),
        ("staged source record", staged_source_path),
        ("frozen first-NVT MDP", frozen_first_nvt),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty {label}: {path}")
    if seed_record_path.is_file() and frozen_first_nvt.is_file():
        try:
            seed_record = load_json(seed_record_path)
            expected_seed = int(realization.get("velocity_seed"))
            if seed_record.get("velocity_seed") != expected_seed:
                errors.append(f"velocity seed record differs from manifest for {rid}")
            if seed_record.get("derived_sha256") != sha256(frozen_first_nvt):
                errors.append(f"frozen first-NVT MDP SHA-256 differs for {rid}")
            match = re.search(r"^\s*gen[_-]seed\s*=\s*([0-9]+)", frozen_first_nvt.read_text(encoding="utf-8", errors="replace"), re.IGNORECASE | re.MULTILINE)
            if match is None or int(match.group(1)) != expected_seed:
                errors.append(f"frozen first-NVT MDP does not contain the manifest seed for {rid}")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"cannot validate velocity seed record for {rid}: {exc}")
    if staged_source_path.is_file():
        try:
            staged = load_json(staged_source_path)
            if staged.get("construction_id") != "build01":
                errors.append(f"{rid} staged source does not reference build01")
            if staged.get("common_minimization_sha256") != common_minimization_sha256:
                errors.append(f"{rid} does not use the frozen common minimization")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"cannot validate staged-source record for {rid}: {exc}")

    tpr_record: dict[str, Any] | None = None
    if required["tpr"].is_file():
        tpr_record, tpr_errors = audit_production_tpr_record(manifest, run_dir, rid, required["tpr"])
        errors.extend(tpr_errors)
    command_records, command_errors = audit_append_only_command_chain(run_dir)
    errors.extend(command_errors)
    production_mdruns = [
        item["argv"] for item in command_records
        if isinstance(item.get("argv"), list) and len(item["argv"]) > 1 and
        item["argv"][1] == "mdrun" and "production" in item["argv"]
    ]
    if phase == "canary":
        canary_commands = [argv for argv in production_mdruns if "-nsteps" in argv and argv[argv.index("-nsteps") + 1] == "2500000"]
        if not canary_commands or len(canary_commands) != len(production_mdruns):
            errors.append(f"{rid} requires one or more restart-safe canary segments, all retaining the 2,500,000-step absolute target")
    elif phase == "production":
        continuation_commands = [argv for argv in production_mdruns if "-cpi" in argv and "-append" in argv]
        if not continuation_commands:
            errors.append(f"{rid} has no retained same-TPR checkpoint-continuation command")

    text_files = [path for path in (work / f"{stem}.log",) if path.is_file()]
    command_logs = run_dir / "command_logs"
    if command_logs.is_dir():
        text_files.extend(sorted(command_logs.glob("*.log")))
    combined = "\n".join(path.read_text(encoding="utf-8", errors="replace") for path in text_files)
    for pattern in NEVER_ALLOW:
        hits = pattern.findall(combined)
        if hits:
            errors.append(f"forbidden runtime pattern {pattern.pattern!r}: {len(hits)} occurrence(s)")
    warning_count = len(GROMPP_WARNING.findall(combined))
    if warning_count:
        errors.append(f"grompp reported {warning_count} numbered warning(s); -maxwarn is prohibited")
    if combined and "Finished mdrun" not in combined:
        errors.append("no 'Finished mdrun' completion marker in retained logs")

    production_command_records = [
        item for item in command_records
        if isinstance(item.get("argv"), list) and len(item["argv"]) > 1
        and item["argv"][1] == "mdrun" and "production" in item["argv"]
    ]
    final_coordinates_endpoint_binding = False
    if production_command_records and required["final_coordinates"].is_file():
        latest_outputs = production_command_records[-1].get("runtime_outputs")
        gro_record = latest_outputs.get("gro", {}) if isinstance(latest_outputs, dict) else {}
        final_coordinates_endpoint_binding = (
            isinstance(gro_record, dict)
            and gro_record.get("sha256") == sha256(required["final_coordinates"])
            and production_command_records[-1].get("returncode") == 0
        )
        if not final_coordinates_endpoint_binding and strict:
            errors.append("final GRO is not bound to the latest successful production mdrun record")
        elif not final_coordinates_endpoint_binding:
            warnings.append("final GRO was not bound to the latest mdrun because strict provenance was not requested")
    elif strict:
        errors.append("strict audit requires a latest successful production mdrun record for the final GRO")

    checks: list[dict[str, Any]] = []
    endpoint_times_ps: dict[str, float | None] = {
        "trajectory": None,
        "energy": None,
        "log": None,
        "final_coordinates": None,
        "checkpoint": None,
    }
    endpoint_steps: dict[str, int | None] = {"log": None, "checkpoint": None}
    if use_gmx and all(path.is_file() for path in required.values()):
        checks = [
            run_check([gmx, "check", "-f", required["trajectory"].name, "-s1", required["tpr"].name], work),
            run_check([gmx, "check", "-e", required["energy"].name], work),
            run_check([gmx, "dump", "-cp", required["checkpoint"].name], work),
        ]
        for check in checks:
            if check["returncode"] != 0:
                errors.append(f"gmx check failed: {check['argv']}")
            output = check["stdout"] + "\n" + check["stderr"]
            for pattern in NEVER_ALLOW:
                if pattern.search(output):
                    errors.append(f"gmx check emitted forbidden pattern {pattern.pattern!r}")
        trajectory_output = checks[0]["stdout"] + "\n" + checks[0]["stderr"]
        energy_output = checks[1]["stdout"] + "\n" + checks[1]["stderr"]
        checkpoint_output = checks[2]["stdout"] + "\n" + checks[2]["stderr"]
        endpoint_times_ps["trajectory"] = last_float(LAST_FRAME, trajectory_output)
        endpoint_times_ps["energy"] = last_float(LAST_ENERGY_FRAME, energy_output)
        endpoint_times_ps["final_coordinates"] = gro_time_ps(required["final_coordinates"])
        endpoint_times_ps["checkpoint"] = last_float(CHECKPOINT_TIME, checkpoint_output)
        checkpoint_steps = CHECKPOINT_STEP.findall(checkpoint_output)
        endpoint_steps["checkpoint"] = int(checkpoint_steps[-1]) if checkpoint_steps else None
        endpoint_steps["log"], endpoint_times_ps["log"] = last_log_step_time(required["log"])
        trajectory_atom_matches = TRAJECTORY_ATOMS.findall(trajectory_output)
        try:
            gro_atoms, gro_box = gro_atom_count_and_box(required["final_coordinates"])
            if not trajectory_atom_matches or gro_atoms != int(trajectory_atom_matches[-1]):
                errors.append("final GRO atom count differs from the validated trajectory")
            if any(value <= 0.0 for value in gro_box[:3]):
                errors.append("final GRO has a nonpositive primary box vector")
        except (OSError, UnicodeError, ValueError) as exc:
            errors.append(f"final GRO structural parse failed: {exc}")
        for label, value in endpoint_times_ps.items():
            if label == "final_coordinates" and value is None and final_coordinates_endpoint_binding:
                continue
            if value is None:
                errors.append(f"could not parse the final {label} time")
        for label, value in endpoint_steps.items():
            if value is None:
                errors.append(f"could not parse the final {label} step")
    elif strict:
        errors.append("strict audit requires the GROMACS executable and all raw outputs for gmx check")
    else:
        warnings.append("gmx check was not executed")

    expected_ps = expected_ns * 1000.0
    expected_step = int(round(expected_ps / float(manifest.get("simulation", {}).get("time_step_ps", 0.002))))
    for label, final_time_ps in endpoint_times_ps.items():
        if final_time_ps is not None and (
            not math.isfinite(final_time_ps) or abs(final_time_ps - expected_ps) > 1e-3
        ):
            errors.append(f"{label} ended at {final_time_ps} ps, expected exactly {expected_ps} ps")
    for label, final_step in endpoint_steps.items():
        if final_step is not None and final_step != expected_step:
            errors.append(f"{label} ended at step {final_step}, expected exactly {expected_step}")

    return {
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "realization_id": rid,
        "phase": phase,
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "warnings": warnings,
        "final_time_ps": endpoint_times_ps["trajectory"],
        "endpoint_times_ps": endpoint_times_ps,
        "endpoint_steps": endpoint_steps,
        "final_coordinates_endpoint_binding": final_coordinates_endpoint_binding,
        "artifacts": artifacts,
        "production_tpr_record": tpr_record,
        "append_only_command_records": command_records,
        "gmx_checks": checks,
    }


def synthetic_self_test() -> None:
    projected, multiplier = projected_storage_bytes(123, 500.0, 5.0)
    if projected != 12_300 or multiplier != 100.0:
        raise RuntimeError("The 5-to-500 ns storage projection is not exactly 100-fold")
    if last_float(LAST_FRAME, "Last frame 5000 time 500000.000") != 500000.0:
        raise RuntimeError("Trajectory endpoint parser failed")
    if last_float(LAST_ENERGY_FRAME, "Last energy frame read 50000 time 500000.000") != 500000.0:
        raise RuntimeError("Energy endpoint parser failed")
    if last_float(CHECKPOINT_TIME, "step = 250000000\nt = 500000.000\n") != 500000.0:
        raise RuntimeError("Checkpoint endpoint parser failed")
    if TRAJECTORY_ATOMS.findall("# Atoms  10839\n") != ["10839"]:
        raise RuntimeError("Trajectory atom-count parser failed")
    with tempfile.TemporaryDirectory(prefix="single_system_md_output_audit_") as temporary:
        root = Path(temporary)
        config = root / "config"
        config.mkdir()
        common_run = root / "runs" / SYSTEM_ID / "common_minimization"
        common_work = common_run / "work"
        common_work.mkdir(parents=True)
        common_coordinates = common_work / "step6.0_minimization.gro"
        common_coordinates.write_bytes(b"synthetic common minimized coordinates")
        common_hash = sha256(common_coordinates)
        for suffix in ("mdp", "tpr", "edr"):
            (common_work / f"step6.0_minimization.{suffix}").write_bytes(f"synthetic common {suffix}".encode())
        (common_work / "step6.0_minimization.log").write_text("Finished mdrun\n", encoding="utf-8")
        realizations: list[dict[str, Any]] = []
        for realization_index, rid in enumerate(REALIZATION_IDS):
            run_directory = root / "runs" / SYSTEM_ID / rid
            work = run_directory / "work"
            work.mkdir(parents=True)
            for suffix in ("tpr", "xtc", "edr", "gro", "cpt"):
                (work / f"production.{suffix}").write_bytes(f"synthetic {rid} {suffix}".encode())
            (work / "production.log").write_text("Finished mdrun\n", encoding="utf-8")
            seed = 26080801 + realization_index
            first_nvt = work / "step6.1_equilibration.frozen.mdp"
            first_nvt.write_text(f"gen_seed = {seed}\n", encoding="utf-8")
            (run_directory / "velocity_seed_record.json").write_text(json.dumps({
                "velocity_seed": seed,
                "derived_sha256": sha256(first_nvt),
            }, indent=2) + "\n", encoding="utf-8")
            (run_directory / "staged_source.json").write_text(json.dumps({
                "construction_id": "build01",
                "common_minimization_sha256": common_hash,
            }, indent=2) + "\n", encoding="utf-8")
            realizations.append({"id": rid, "velocity_seed": seed, "run_directory": str(run_directory.relative_to(root))})
        manifest = {
            "simulation": {"production_ns": 500.0},
            "systems": [{
                "id": SYSTEM_ID,
                "construction": {
                    "id": "build01",
                    "charmm_gui_archive": {"sha256": "a" * 64},
                    "gromacs_input_tree_manifest": {"sha256": "b" * 64},
                    "minimization_mdp": {
                        "path": str((common_work / "step6.0_minimization.mdp").relative_to(root)),
                        "sha256": sha256(common_work / "step6.0_minimization.mdp"),
                    },
                    "common_minimization_run_directory": str(common_run.relative_to(root)),
                },
                "realizations": realizations,
            }],
        }
        manifest_path = config / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        build_hash = build_contract_sha256(manifest)
        common_logs = common_run / "command_logs"; common_logs.mkdir()
        common_finished: list[Path] = []
        previous_common_hash: str | None = None
        for sequence, argv in enumerate((
            ["gmx", "grompp", "-f", "step6.0_minimization.mdp", "-o", "step6.0_minimization.tpr"],
            ["gmx", "mdrun", "-deffnm", "step6.0_minimization"],
        ), start=1):
            started = common_logs / f"{sequence:04d}_common.started.json"
            started.write_text(json.dumps({"sequence": sequence, "argv": argv, "previous_finished_record_sha256": previous_common_hash}), encoding="utf-8")
            stdout = common_logs / f"{sequence:04d}_common.stdout.log"; stdout.write_text("Finished mdrun\n", encoding="utf-8")
            stderr = common_logs / f"{sequence:04d}_common.stderr.log"; stderr.write_text("", encoding="utf-8")
            finished = common_logs / f"{sequence:04d}_common.finished.json"
            finished.write_text(json.dumps({
                "sequence": sequence, "started_record": str(started), "started_record_sha256": sha256(started),
                "returncode": 0, "stdout": {"path": str(stdout), "sha256": sha256(stdout)},
                "stderr": {"path": str(stderr), "sha256": sha256(stderr)},
            }), encoding="utf-8")
            previous_common_hash = sha256(finished); common_finished.append(finished)
        common_artifacts = {
            label: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for label, path in {
                "mdp": common_work / "step6.0_minimization.mdp", "tpr": common_work / "step6.0_minimization.tpr",
                "energy": common_work / "step6.0_minimization.edr", "log": common_work / "step6.0_minimization.log",
                "coordinates": common_coordinates,
            }.items()
        }
        common_record = {
            "schema_version": "2.0", "report_type": "common_minimization_validation", "status": "pass",
            "technical_integrity_pass": True, "system_id": SYSTEM_ID, "construction_id": "build01",
            "build_contract_sha256": build_hash, "archive_sha256": "a" * 64,
            "gromacs_input_tree_manifest_sha256": "b" * 64,
            "minimization_mdp_sha256": manifest["systems"][0]["construction"]["minimization_mdp"]["sha256"],
            "coordinates_sha256": common_hash, "artifacts": common_artifacts,
            "command_records": [{"path": str(path), "sha256": sha256(path)} for path in common_finished],
            "grompp_warning_count": 0, "forbidden_runtime_pattern_count": 0, "errors": [],
            "integrity": {"payload_sha256": "UNSEALED"},
        }
        common_record["integrity"]["payload_sha256"] = report_payload_sha256(common_record, ("integrity", "payload_sha256"))
        (common_run / "common_minimization_record.json").write_text(json.dumps(common_record, indent=2) + "\n", encoding="utf-8")
        for rid in REALIZATION_IDS:
            run_directory = root / "runs" / SYSTEM_ID / rid
            work = run_directory / "work"
            command_logs = run_directory / "command_logs"; command_logs.mkdir()
            previous_hash: str | None = None
            finished_paths: list[Path] = []
            for sequence, argv in enumerate((
                ["gmx", "grompp", "-f", "production.frozen.mdp", "-o", "production.tpr"],
                ["gmx", "mdrun", "-s", "production.tpr", "-deffnm", "production", "-cpi", "production.cpt", "-append"],
            ), start=1):
                started = command_logs / f"{sequence:04d}_synthetic.started.json"
                started.write_text(json.dumps({
                    "sequence": sequence, "argv": argv, "previous_finished_record_sha256": previous_hash,
                }), encoding="utf-8")
                stdout = command_logs / f"{sequence:04d}_synthetic.stdout.log"; stdout.write_text("synthetic\n", encoding="utf-8")
                stderr = command_logs / f"{sequence:04d}_synthetic.stderr.log"; stderr.write_text("", encoding="utf-8")
                finished = command_logs / f"{sequence:04d}_synthetic.finished.json"
                finished_payload = {
                    "sequence": sequence, "started_record": str(started),
                    "started_record_sha256": sha256(started), "returncode": 0,
                    "stdout": {"path": str(stdout), "sha256": sha256(stdout)},
                    "stderr": {"path": str(stderr), "sha256": sha256(stderr)},
                }
                if argv[1] == "mdrun":
                    finished_payload["runtime_outputs"] = {
                        suffix: {
                            "path": str(work / f"production.{suffix}"),
                            "bytes": (work / f"production.{suffix}").stat().st_size,
                            "sha256": sha256(work / f"production.{suffix}"),
                        }
                        for suffix in ("xtc", "edr", "log", "gro", "cpt")
                    }
                finished.write_text(json.dumps(finished_payload), encoding="utf-8")
                previous_hash = sha256(finished); finished_paths.append(finished)
            tpr = work / "production.tpr"
            (run_directory / "production_tpr_record.json").write_text(json.dumps({
                "system_id": SYSTEM_ID, "construction_id": "build01", "realization_id": rid,
                "build_contract_sha256": build_hash, "production_tpr_sha256": sha256(tpr),
                "grompp_finished_record": str(finished_paths[0]),
                "grompp_finished_record_sha256": sha256(finished_paths[0]),
            }), encoding="utf-8")
        report = root / "pass_report.json"
        base = [sys.executable, str(Path(__file__).resolve()), "--manifest", str(manifest_path), "--phase", "production", "--skip-gmx-check"]
        passed = subprocess.run(base + ["--report", str(report)], capture_output=True, text=True, check=False)
        if passed.returncode != 0 or load_json(report).get("status") != "pass":
            raise RuntimeError(f"Synthetic complete-output audit failed:\n{passed.stdout}\n{passed.stderr}")
        (root / "runs" / SYSTEM_ID / "rep03" / "work" / "production.xtc").unlink()
        failed_report = root / "fail_report.json"
        failed = subprocess.run(base + ["--report", str(failed_report)], capture_output=True, text=True, check=False)
        if failed.returncode == 0 or load_json(failed_report).get("status") != "fail":
            raise RuntimeError("Synthetic missing-output audit did not fail closed")
    print("SELF-TEST PASS: complete rep01-rep03 outputs pass non-GROMACS audit, endpoint parsers are exact, storage scaling is 100-fold, and a missing trajectory fails closed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--phase", choices=("equilibration", "canary", "production"), default="production")
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--skip-gmx-check", action="store_true", help="Allowed only outside strict mode")
    parser.add_argument("--self-test", action="store_true", help="Run a synthetic fail-closed output audit")
    args = parser.parse_args()
    if args.self_test:
        if args.manifest is not None or args.strict:
            raise SystemExit("--self-test cannot be combined with manifest or strict mode")
        synthetic_self_test()
        return 0
    if args.manifest is None:
        parser.error("--manifest is required unless --self-test is used")
    if args.strict and args.skip_gmx_check:
        raise SystemExit("--strict cannot be combined with --skip-gmx-check")

    manifest_path = args.manifest.resolve()
    package_root = manifest_path.parent.parent.resolve()
    manifest = load_json(manifest_path)
    gmx_available = shutil.which(args.gmx) is not None
    if args.strict and not gmx_available:
        raise SystemExit(f"GROMACS executable not found: {args.gmx}")
    use_gmx = gmx_available and not args.skip_gmx_check
    gmx_identity: dict[str, Any] | None = None
    if args.strict:
        try:
            gmx_identity = validate_bound_gmx_identity(package_root, manifest, args.gmx)
        except (OSError, ValueError, RuntimeError) as exc:
            raise SystemExit(f"Bound GROMACS executable validation failed: {exc}")
    systems = manifest.get("systems", [])
    if not isinstance(systems, list) or len(systems) != 1 or systems[0].get("id") != SYSTEM_ID:
        raise SystemExit(f"Manifest must contain exactly one system: {SYSTEM_ID}")
    realizations = systems[0].get("realizations", [])
    if not isinstance(realizations, list) or [item.get("id") for item in realizations] != list(REALIZATION_IDS):
        raise SystemExit(f"Manifest must contain exactly {list(REALIZATION_IDS)}")
    seeds = [item.get("velocity_seed") for item in realizations]
    if not all(isinstance(seed, int) and seed > 0 for seed in seeds) or len(set(seeds)) != 3:
        raise SystemExit("rep01-rep03 must have three distinct positive velocity seeds")
    simulation = manifest.get("simulation", {})
    production_ns = float(simulation.get("production_ns", 0.0))
    if abs(production_ns - 500.0) > 1e-12:
        raise SystemExit("Production duration must be exactly 500 ns")
    construction = systems[0].get("construction", {})
    common_run = inside(package_root, str(construction.get("common_minimization_run_directory", "")))
    common_coordinates = common_run / "work" / "step6.0_minimization.gro"
    common_record_path = common_run / "common_minimization_record.json"
    if not common_coordinates.is_file() or not common_record_path.is_file():
        raise SystemExit("Validated common-minimization coordinates and record are required")
    common_record = load_json(common_record_path)
    common_hash = sha256(common_coordinates)
    common_errors = validate_common_minimization_record(package_root, manifest, construction, common_record_path)
    if common_errors or common_record.get("coordinates_sha256") != common_hash:
        raise SystemExit("Common-minimization record is invalid: " + "; ".join(common_errors))
    if args.phase == "equilibration":
        results = [
            audit_equilibration_one(
                package_root, manifest, realization, construction, common_hash, args.gmx, use_gmx, args.strict
            )
            for realization in realizations
        ]
        expected_ns = None
    else:
        expected_ns = 5.0 if args.phase == "canary" else production_ns
        if args.phase == "canary":
            release = simulation.get("release_contract", {})
            if release.get("canary_target_ns_per_realization") != 5.0 or release.get("canary_total_ns") != 15.0:
                raise SystemExit("Canary contract must be exactly 5 ns per realization / 15 ns total")
        results = [
            audit_one(package_root, manifest, realization, args.phase, expected_ns, common_hash, args.gmx, use_gmx, args.strict)
            for realization in realizations
        ]
    storage: dict[str, Any] | None = None
    if args.phase == "canary":
        canary_dynamic_bytes = sum(
            int(artifact.get("bytes", 0))
            for result in results
            for label, artifact in result.get("artifacts", {}).items()
            if label in {"trajectory", "energy", "log", "final_coordinates", "checkpoint"}
        )
        canary_ns_per_realization = float(
            simulation.get("release_contract", {}).get("canary_target_ns_per_realization", 0.0)
        )
        if canary_ns_per_realization <= 0.0 or not math.isfinite(canary_ns_per_realization):
            raise SystemExit("Canary duration is invalid for storage projection")
        projected_dynamic_bytes, projection_multiplier = projected_storage_bytes(
            canary_dynamic_bytes, production_ns, canary_ns_per_realization
        )
        if abs(projection_multiplier - 100.0) > 1e-12:
            raise SystemExit("The 5-to-500 ns storage projection multiplier must be exactly 100")
        required_with_headroom_bytes = math.ceil(projected_dynamic_bytes * 1.30)
        budget_value = simulation.get("release_contract", {}).get("storage_budget_bytes")
        try:
            storage_budget_bytes = int(budget_value)
        except (TypeError, ValueError):
            storage_budget_bytes = -1
        storage_probe_path = inside(package_root, str(realizations[0].get("run_directory", "")))
        live_usage = shutil.disk_usage(storage_probe_path)
        live_device_id = int(os.stat(storage_probe_path).st_dev)
        live_free_bytes = int(live_usage.free)
        storage = {
            "canary_dynamic_output_bytes_all_three": canary_dynamic_bytes,
            "canary_ns_per_realization": canary_ns_per_realization,
            "production_ns_per_realization": production_ns,
            "projection_multiplier": projection_multiplier,
            "linear_projection_to_3x500ns_bytes": projected_dynamic_bytes,
            "required_with_30_percent_headroom_bytes": required_with_headroom_bytes,
            "predeclared_storage_budget_bytes": storage_budget_bytes,
            "live_storage_probe": {
                "path": str(storage_probe_path),
                "device_id": live_device_id,
                "total_bytes": int(live_usage.total),
                "used_bytes": int(live_usage.used),
                "free_bytes": live_free_bytes,
                "measured_at_utc": utc_now(),
            },
            "status": "pass" if (
                storage_budget_bytes >= required_with_headroom_bytes
                and live_free_bytes >= required_with_headroom_bytes
                and projected_dynamic_bytes > 0
            ) else "fail",
            "note": "Projection is a technical release estimate from retained 3x5 ns output; it is not a scientific threshold.",
        }
        if storage["status"] != "pass":
            for result in results:
                result["errors"].append(
                    "predeclared budget and live filesystem free space must both cover the measured canary projection plus 30% headroom"
                )
                result["status"] = "fail"
    report = {
        "schema_version": "2.0",
        "report_type": "md_stage_output_validation",
        "created_at_utc": utc_now(),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256(manifest_path),
        "manifest_release_contract_sha256": release_manifest_contract_sha256(manifest),
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "construction_archive_sha256": systems[0]["construction"]["charmm_gui_archive"]["sha256"],
        "build_contract_sha256": build_contract_sha256(manifest),
        "build_validation_report_sha256": construction.get("build_validation_report", {}).get("sha256"),
        "phase": args.phase,
        "strict": args.strict,
        "gmx_executable": shutil.which(args.gmx),
        "gmx_identity": gmx_identity,
        "status": "pass" if all(item["status"] == "pass" for item in results) else "fail",
        "analysis_disposition": (
            ({"equilibration": "eligible_for_5ns_canary", "canary": "eligible_for_checkpoint_continuation", "production": "eligible_for_downstream_qc"}[args.phase])
            if all(item["status"] == "pass" for item in results) else "inconclusive_fail_closed"
        ),
        "qc_and_stationarity_gate": {
            "status": "not_evaluated" if all(item["status"] == "pass" for item in results) else "fail",
            "failure_policy": "inconclusive_if_any_realization_fails_qc_or_stationarity",
            "results": [
                {
                    "realization_id": item["realization_id"],
                    "qc_status": "not_evaluated" if item["status"] == "pass" else "fail",
                    "stationarity_status": "not_evaluated" if item["status"] == "pass" else "fail",
                }
                for item in results
            ],
        },
        "runs": results,
        "storage_release_check": storage,
        "integrity": {"payload_sha256": "UNSEALED"},
    }
    report["integrity"]["payload_sha256"] = report_payload_sha256(report, ("integrity", "payload_sha256"))
    default_name = f"{args.phase}_output_validation.json"
    report_path = args.report.resolve() if args.report else package_root / "reports" / default_name
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": str(report_path)}, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
