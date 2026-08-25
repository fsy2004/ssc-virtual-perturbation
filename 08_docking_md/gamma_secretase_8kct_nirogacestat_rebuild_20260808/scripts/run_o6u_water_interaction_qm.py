#!/usr/bin/env python3
"""Preflight or execute the authorized O6U HF/6-31G(d) water-QM batch."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


GENERATION_STATUS = "pass_generation_only_visual_review_required"
AUTHORIZATION_STATUS = "pass_frozen_preqm_orientation_authorization"


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


def verify_record(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    if not path.is_file() or path.stat().st_size != record.get("size_bytes") or sha256(path) != record.get("sha256"):
        raise RuntimeError(f"Artifact failed hash/size verification: {label}")
    return path


def output_for_run_file(path: Path) -> Path:
    return path.with_suffix(".out")


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def run_one(job: dict[str, object], python: str, log_dir: Path) -> dict[str, object]:
    label = str(job["label"])
    run_file = Path(str(job["run_file"])).resolve()
    output_file = output_for_run_file(run_file)
    stdout_path = log_dir / f"{label}.stdout.txt"
    stderr_path = log_dir / f"{label}.stderr.txt"
    started = now()
    start = time.monotonic()
    with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, stderr_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as stderr:
        completed = subprocess.run(
            [python, str(run_file)],
            cwd=run_file.parent,
            stdout=stdout,
            stderr=stderr,
            check=False,
            text=True,
        )
    result: dict[str, object] = {
        "label": label,
        "orientation_id": job.get("orientation_id"),
        "started_at_utc": started,
        "finished_at_utc": now(),
        "elapsed_seconds": time.monotonic() - start,
        "returncode": completed.returncode,
        "run_file": artifact(run_file),
        "stdout": artifact(stdout_path),
        "stderr": artifact(stderr_path),
        "expected_output": str(output_file),
    }
    if output_file.is_file():
        result["raw_output"] = artifact(output_file)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generation-report", required=True, type=Path)
    parser.add_argument("--authorization-report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", required=True, choices=("canary_preflight", "formal_execution"))
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()
    if not 1 <= args.workers <= 20:
        raise SystemExit("Workers must be between 1 and 20")

    generation_path = args.generation_report.resolve()
    authorization_path = args.authorization_report.resolve()
    output_dir = args.output_dir.resolve()
    if not generation_path.is_file() or not authorization_path.is_file():
        raise SystemExit("Generation and authorization reports must exist")
    if output_dir.exists():
        raise SystemExit("Output directory already exists; refusing reuse")
    output_dir.mkdir(parents=True, exist_ok=False)
    report_path = output_dir / "O6U_WATER_INTERACTION_QM_BATCH.json"

    state: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_water_interaction_qm_batch",
        "created_at_utc": now(),
        "updated_at_utc": now(),
        "status": "preflighting",
        "role": args.role,
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "water_interaction_qm_executed": False,
    }
    atomic_json(report_path, state)

    try:
        generation = load_json(generation_path)
        authorization = load_json(authorization_path)
        expected_generation_role = "generation_canary" if args.role == "canary_preflight" else "formal_mp2_target"
        expected_freeze_role = "synthetic_canary_fixture" if args.role == "canary_preflight" else "formal_preqm_authorization"
        if (
            generation.get("status") != GENERATION_STATUS
            or generation.get("role") != expected_generation_role
            or generation.get("production_approved") is not False
        ):
            raise RuntimeError("Generation report differs from the exact role/status gate")
        if authorization.get("freeze_role") != expected_freeze_role or authorization.get("production_approved") is not False:
            raise RuntimeError("Authorization report differs from the exact role/production gate")
        if args.role == "formal_execution":
            if authorization.get("status") != AUTHORIZATION_STATUS or authorization.get("water_interaction_qm_authorized") is not True:
                raise RuntimeError("Formal water-interaction QM is not authorized by the frozen 70-row record")
        else:
            if authorization.get("status") != "pass_synthetic_canary_structure_only" or authorization.get("water_interaction_qm_authorized") is not False:
                raise RuntimeError("Canary preflight requires a non-authorizing synthetic freeze record")

        auth_inputs = authorization.get("inputs")
        if not isinstance(auth_inputs, dict):
            raise RuntimeError("Authorization input map is missing")
        template_path = verify_record(auth_inputs.get("template_report"), "authorization.template_report")
        template = load_json(template_path)
        template_inputs = template.get("inputs")
        if not isinstance(template_inputs, dict):
            raise RuntimeError("Template input map is missing")
        if verify_record(template_inputs.get("generation_report"), "template.generation_report") != generation_path:
            raise RuntimeError("Authorization chain is bound to a different generation report")

        selected_ids = authorization.get("run_qm_orientation_ids")
        if not isinstance(selected_ids, list) or not selected_ids or len(selected_ids) != len(set(selected_ids)):
            raise RuntimeError("Selected water-QM orientation list is blank or contains duplicates")
        pairs = generation.get("generated_pairs")
        monomers = generation.get("monomer_files")
        if not isinstance(pairs, list) or not isinstance(monomers, list) or len(monomers) != 2:
            raise RuntimeError("Generated pair or monomer file registry is incomplete")
        pair_by_id: dict[str, dict[str, object]] = {}
        for pair in pairs:
            if not isinstance(pair, dict):
                raise RuntimeError("Generated pair record is not an object")
            orientation_id = str(pair.get("orientation_id", ""))
            if not orientation_id or orientation_id in pair_by_id:
                raise RuntimeError("Generated pair orientation identity is blank or duplicate")
            pair_by_id[orientation_id] = pair
        if any(orientation_id not in pair_by_id for orientation_id in selected_ids):
            raise RuntimeError("Authorization selects an orientation absent from generated inputs")

        jobs: list[dict[str, object]] = []
        for index, monomer in enumerate(monomers, start=1):
            run_file = verify_record(monomer, f"generation.monomer{index}")
            jobs.append({"label": f"monomer{index}", "orientation_id": None, "run_file": str(run_file)})
        planned_pairs: list[dict[str, object]] = []
        for orientation_id in selected_ids:
            pair = pair_by_id[orientation_id]
            run_file = verify_record(pair.get("run_file"), f"{orientation_id}.run_file")
            coordinate_file = verify_record(pair.get("coordinate_file"), f"{orientation_id}.coordinate_file")
            grid = pair.get("distance_grid_angstrom")
            expected_grid = [round(1.5 + 0.05 * index, 2) for index in range(31)]
            if grid != expected_grid:
                raise RuntimeError(f"Distance grid differs for {orientation_id}")
            planned_pairs.append(
                {
                    "orientation_id": orientation_id,
                    "run_file": artifact(run_file),
                    "coordinate_file": artifact(coordinate_file),
                    "distance_grid_angstrom": grid,
                }
            )
            jobs.append({"label": orientation_id, "orientation_id": orientation_id, "run_file": str(run_file)})

        state.update(
            {
                "status": "pass_preflight_no_qm_executed" if args.role == "canary_preflight" else "executing_authorized_raw_water_qm",
                "updated_at_utc": now(),
                "generation_report": artifact(generation_path),
                "authorization_report": artifact(authorization_path),
                "selected_orientation_count": len(selected_ids),
                "selected_orientation_ids": selected_ids,
                "planned_pairs": planned_pairs,
                "monomer_run_files": [artifact(Path(str(job["run_file"]))) for job in jobs[:2]],
                "workers": args.workers,
            }
        )
        atomic_json(report_path, state)
        if args.role == "canary_preflight":
            state["release_boundary"] = (
                "This preflight validated the complete authorization chain and exact input hashes. It executed no QM "
                "and cannot authorize formal execution, parameter fitting, CHARMM-GUI construction, or MD."
            )
            atomic_json(report_path, state)
            print(json.dumps({"status": state["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
            return 0

        for job in jobs:
            raw_output = output_for_run_file(Path(str(job["run_file"])))
            if raw_output.exists():
                raise RuntimeError(f"Raw output already exists; refusing reuse: {raw_output}")
        log_dir = output_dir / "process_logs"
        log_dir.mkdir()
        results: list[dict[str, object]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(args.workers, len(jobs))) as executor:
            futures = {executor.submit(run_one, job, sys.executable, log_dir): job for job in jobs}
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
                state.update({"updated_at_utc": now(), "completed_job_count": len(results)})
                atomic_json(report_path, state)
        results.sort(key=lambda item: str(item["label"]))

        failures: list[str] = []
        nan_pattern = re.compile(r"(?<![A-Za-z])NaN(?![A-Za-z])", re.IGNORECASE)
        for result in results:
            label = str(result["label"])
            raw_record = result.get("raw_output")
            if result.get("returncode") != 0 or not isinstance(raw_record, dict):
                failures.append(f"{label}:process_or_output")
                continue
            raw_path = verify_record(raw_record, f"{label}.raw_output")
            text = raw_path.read_text(encoding="utf-8", errors="replace")
            if nan_pattern.search(text) or "PsiException" in text or "Traceback (most recent call last)" in text:
                failures.append(f"{label}:fatal_marker")
            if label.startswith("monomer"):
                if "INTERACTION MONOMER ENERGY is" not in text:
                    failures.append(f"{label}:missing_monomer_energy")
            elif not all(token in text for token in ("INTERACTION TABLE NOBSSE START", "INTERACTION TABLE NOBSSE END", "INTERACTION DISTANCE and ENERGY are:")):
                failures.append(f"{label}:missing_interaction_table")
        state.update(
            {
                "updated_at_utc": now(),
                "jobs": results,
                "completed_job_count": len(results),
                "water_interaction_qm_executed": True,
                "technical_failures": failures,
            }
        )
        if failures:
            state["status"] = "fail_closed_raw_water_qm_batch"
            state["release_boundary"] = "Every failed attempt remains recorded. No failed orientation may be dropped or used for fitting."
            atomic_json(report_path, state)
            return 2
        state["status"] = "pass_raw_water_qm_execution_outputs_present_validation_required"
        state["release_boundary"] = (
            "All authorized raw tasks produced technically complete outputs. Independent numerical reconstruction and "
            "post-QM applicable/weak/unfavourable/excluded dispositions are still required before parameter fitting."
        )
        atomic_json(report_path, state)
        print(json.dumps({"status": state["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
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
        atomic_json(report_path, state)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
