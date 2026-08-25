#!/usr/bin/env python3
"""Fail-closed validation for the gamma-secretase rebuild manifest.

The validator deliberately rejects the template. It becomes executable only after
all scientific decisions and build artifacts have been frozen and evidenced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Iterable

from hash_tree_manifest import validate_report as validate_tree_manifest_report
from md_contract import (
    analysis_plan_contract_sha256,
    build_contract_sha256,
    release_manifest_contract_sha256,
    report_payload_sha256,
    validate_common_minimization_record,
    validate_production_mdp,
)


PRIMARY_ID = "8kct_nirogacestat_native"
REALIZATION_IDS = ("rep01", "rep02", "rep03")
O6U_CCD_SHA256 = "c5c3c1ac73c9cb512612a855e39caddd9d840eaf823822f677c2675256659acd"
O6U_OFFICIAL_SDF_SHA256 = "985ac8899eb97efff73a57682da4794fbe55acf06be25717f9d885d64cb5962a"
O6U_PARAMETER_INPUT_SDF_SHA256 = "2cb9d769cde4157181a6199b83294cad56cade14ab34a5e86a6deb6790fc28d5"
O6U_PREMAPPING_SHA256 = "a773708cf030ba03cdee924a51359b84810742fbce63b99cc95b8d6309e2eca3"
O6U_CGENFF_CORRESPONDENCE_TSV_SHA256 = "62b5a9500a0c5e0c2d85eb3fa51fb4e4cb82881dafd3d0b2b38f02231d6935f5"
O6U_WATER_ORIENTATION_DA_SHA256 = "5ea7e12c464750b8f35d9fa3feed875bb6b69a4091867b8fb6ddf2ad3c6272a3"
O6U_PARAMETERIZATION_TOOLCHAIN_RECORD_SHA256 = "7cc4888eb93d4c805cbdc6d5d769ec20a6af479415daf65a631e413b6af6afc7"
O6U_PARAMETERIZATION_TOOLCHAIN_VALIDATION_SHA256 = "c3abc9cc6a8f1dbbe326fa28fc3c677bf42cffb874ae91f863541549454a6a5a"
O6U_WATER_PRECISION_AUDIT_SHA256 = "4726b3623f650a4d2ccb6dab56559aa10df4b882f15ed6d3e6ad51383b1e3886"
O6U_WATER_DISPOSITION_POLICY_SHA256 = "4e65a9ecb9d90ec9d0c9ee849c9c05f1fd0aae09c7e16a74bb59d8fb5043fccb"
O6U_WATER_CHEMICAL_PRESCREEN_SHA256 = "976211319f88c7c6ef19779d2e0e4ac011eed443ece525d7fb2c030556deb96a"
KCT8_CIF_SHA256 = "2ed75442ca2c503a014b4e5e8bac67e107201c31776238d0ae94069b39013da9"
LIGAND_TOOLCHAIN_ROLES = {
    "cgenff_initial_assignment", "ffparam_optimization", "qm_target_generation",
    "charmm_reference_energy", "gromacs_regression_and_canary",
}
LIGAND_COMMAND_ROLES = {
    "initial_assignment", "qm_target_generation", "ffparam_fit_and_validation",
    "charmm_energy_regression", "gromacs_energy_regression", "cheap_pre_membrane_canaries",
}
TODO_RE = re.compile(r"\b(?:TODO|TBD|PENDING|PLACEHOLDER)\b", re.IGNORECASE)
SECRET_KEY_RE = re.compile(r"(?:password|passwd|bearer|api[_-]?key|session[_-]?token|secret)", re.IGNORECASE)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
ALLOWED_DYAD_STATES = {"deprotonated", "protonated_od1", "protonated_od2"}
RELEASE_STATUS_ORDER = (
    "design_frozen",
    "chemistry_approved",
    "build_approved",
    "minimization_approved",
    "equilibration_approved",
    "ready_for_production",
)
STAGE_MINIMUM_STATUS = {
    "design": "design_frozen",
    "chemistry": "chemistry_approved",
    "builds": "build_approved",
    "equilibration": "minimization_approved",
    "canary": "equilibration_approved",
    "production": "ready_for_production",
}
FROZEN_PRODUCTION_PROTOCOL = (
    Path(__file__).resolve().parent.parent
    / "config"
    / "production_protocol_hmr4fs_303K_v1.json"
)


class Audit:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        if not condition:
            self.errors.append(message)

    def warn(self, condition: bool, message: str) -> None:
        if not condition:
            self.warnings.append(message)


def validate_frozen_production_protocol_fields(manifest: dict[str, Any], audit: Audit) -> None:
    """Reject drift from the protocol frozen before the completed production runs."""
    protocol = load_json(FROZEN_PRODUCTION_PROTOCOL)
    model = manifest.get("global_model", {})
    simulation = manifest.get("simulation", {})
    expected_temperature = float(protocol["ensemble"]["temperature_K"])
    expected_timestep = float(protocol["ensemble"]["production_time_step_ps"])
    expected_duration = float(protocol["production"]["duration_ns_per_realization"])
    expected_hmr = bool(protocol["force_field"]["hydrogen_mass_repartitioning"])

    temperature = finite_number(model.get("temperature_k"), "temperature_k", audit)
    if temperature is not None:
        audit.require(
            abs(temperature - expected_temperature) < 1e-12,
            f"temperature must match frozen protocol value {expected_temperature} K",
        )
    timestep = finite_number(simulation.get("time_step_ps"), "time_step_ps", audit)
    if timestep is not None:
        audit.require(
            abs(timestep - expected_timestep) < 1e-12,
            f"time step must match frozen protocol value {expected_timestep} ps",
        )
    duration = finite_number(simulation.get("production_ns"), "production_ns", audit)
    if duration is not None:
        audit.require(
            abs(duration - expected_duration) < 1e-12,
            f"production must match frozen protocol duration {expected_duration} ns",
        )
    audit.require(
        simulation.get("hydrogen_mass_repartitioning") is expected_hmr,
        f"hydrogen_mass_repartitioning must be {expected_hmr} for the frozen protocol",
    )


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_number(value: Any, label: str, audit: Audit) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        audit.errors.append(f"{label} must be numeric, got {value!r}")
        return None
    if not math.isfinite(result):
        audit.errors.append(f"{label} must be finite")
        return None
    return result


def contains_todo(value: Any) -> bool:
    if isinstance(value, str):
        return bool(TODO_RE.search(value))
    if isinstance(value, dict):
        return any(contains_todo(item) for item in value.values())
    if isinstance(value, list):
        return any(contains_todo(item) for item in value)
    return False


def scan_secret_keys(value: Any, prefix: str = "manifest") -> Iterable[str]:
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{prefix}.{key}"
            if SECRET_KEY_RE.search(str(key)):
                yield child
            yield from scan_secret_keys(item, child)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from scan_secret_keys(item, f"{prefix}[{index}]")
    elif isinstance(value, str) and value.lower().startswith("bearer "):
        yield prefix


def resolve_inside(package_root: Path, relative: str, label: str, audit: Audit) -> Path | None:
    if not isinstance(relative, str) or not relative.strip() or contains_todo(relative):
        audit.errors.append(f"{label} is missing or unresolved")
        return None
    candidate = (package_root / relative).resolve()
    if candidate != package_root and package_root not in candidate.parents:
        audit.errors.append(f"{label} escapes package root: {relative}")
        return None
    if not candidate.is_file() or candidate.stat().st_size == 0:
        audit.errors.append(f"{label} is missing or empty: {candidate}")
        return None
    return candidate


def validate_artifact(package_root: Path, value: Any, label: str, audit: Audit) -> Path | None:
    if not isinstance(value, dict):
        audit.errors.append(f"{label} must be an object containing path and sha256")
        return None
    path = resolve_inside(package_root, value.get("path", ""), f"{label}.path", audit)
    expected = str(value.get("sha256", "")).lower()
    audit.require(bool(SHA256_RE.fullmatch(expected)), f"{label}.sha256 must be a lowercase SHA-256 digest")
    if path is not None and SHA256_RE.fullmatch(expected):
        audit.require(sha256(path) == expected, f"{label} SHA-256 mismatch: {path}")
    return path


def validate_orientation_record(
    package_root: Path, artifact: Any, construction: dict[str, Any], study_id: str, audit: Audit
) -> None:
    path = validate_artifact(package_root, artifact, "build01 membrane orientation record", audit)
    if path is None:
        return
    try:
        record = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"cannot read membrane orientation record: {exc}")
        return
    audit.require(record.get("schema_version") == "1.0", "orientation record schema must be 1.0")
    audit.require(not contains_todo(record), "orientation record contains unresolved TODO/TBD values")
    audit.require(record.get("record_type") == "single_ppm_membrane_orientation_attestation", "orientation record type is invalid")
    audit.require(record.get("approval_status") == "approved", "orientation record is not approved")
    audit.require(record.get("study_id") == study_id, "orientation record study_id differs from manifest")
    audit.require(record.get("system_id") == PRIMARY_ID and record.get("construction_id") == "build01", "orientation record system/build ID is wrong")
    audit.require(record.get("pdb_reader_jobid") == construction.get("pdb_reader_jobid"), "orientation record PDB Reader job differs")
    audit.require(record.get("quick_bilayer_jobid") == construction.get("quick_bilayer_jobid"), "orientation record Quick Bilayer job differs")
    audit.require(record.get("ppm_application_count") == 1, "orientation record must document exactly one PPM application")
    audit.require(record.get("ppm_applied_by") == "CHARMM-GUI Quick Bilayer ppm=true", "orientation record PPM route is wrong")
    validate_artifact(package_root, record.get("pre_ppm_coordinates"), "orientation pre-PPM coordinates", audit)
    validate_artifact(package_root, record.get("post_ppm_coordinates"), "orientation post-PPM coordinates", audit)
    matrix = record.get("homogeneous_transform_4x4_row_major")
    parsed: list[float] = []
    if not isinstance(matrix, list) or len(matrix) != 16:
        audit.errors.append("orientation transform must contain exactly 16 row-major values")
    else:
        for index, value in enumerate(matrix):
            number = finite_number(value, f"orientation transform[{index}]", audit)
            if number is not None:
                parsed.append(number)
        if len(parsed) == 16:
            audit.require(all(abs(parsed[12 + index] - expected) <= 1e-12 for index, expected in enumerate((0.0, 0.0, 0.0, 1.0))), "orientation transform last row must be [0,0,0,1]")
            rotation = [parsed[0:3], parsed[4:7], parsed[8:11]]
            for row_index, row in enumerate(rotation):
                audit.require(abs(sum(value * value for value in row) - 1.0) <= 1e-6, f"orientation rotation row {row_index} is not unit length")
            for left in range(3):
                for right in range(left + 1, 3):
                    audit.require(abs(sum(rotation[left][i] * rotation[right][i] for i in range(3))) <= 1e-6, "orientation rotation rows are not orthogonal")
            determinant = (
                rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
                - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
                + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
            )
            audit.require(abs(determinant - 1.0) <= 1e-6, "orientation rotation determinant must be +1")
    boundaries = record.get("membrane_boundaries_angstrom", {})
    lower = finite_number(boundaries.get("lower"), "orientation lower membrane boundary", audit) if isinstance(boundaries, dict) else None
    upper = finite_number(boundaries.get("upper"), "orientation upper membrane boundary", audit) if isinstance(boundaries, dict) else None
    if lower is not None and upper is not None:
        audit.require(lower < upper, "orientation membrane boundaries are not ordered")
    audit.require(isinstance(boundaries, dict) and isinstance(boundaries.get("source"), str) and bool(boundaries["source"].strip()), "orientation boundary source is unresolved")
    visual = record.get("visual_evidence", {})
    for view in ("orthogonal_xy", "orthogonal_xz", "orthogonal_yz"):
        validate_artifact(package_root, visual.get(view) if isinstance(visual, dict) else None, f"orientation {view} screenshot", audit)
    request_path = validate_artifact(package_root, record.get("ppm_request_record"), "sanitized PPM request/response record", audit)
    if request_path is not None:
        try:
            request_record = load_json(request_path)
            expected_request = {
                "jobid": construction.get("pdb_reader_jobid"),
                "upper": "POPC=1", "lower": "POPC=1", "margin": "20.0", "wdist": "22.5",
                "Ion_conc": "0.15", "Ion_type": "NaCl", "clone_job": "false", "ppm": "true",
                "topologyIn": "true", "heteroatoms": "true",
            }
            audit.require(request_record.get("endpoint") == "https://www.charmm-gui.org/api/quick_bilayer", "PPM request did not use the pinned official CHARMM-GUI endpoint")
            audit.require(request_record.get("system_id") == PRIMARY_ID and request_record.get("build_id") == "build01", "PPM request system/build identity differs")
            audit.require(request_record.get("request") == expected_request, "PPM request payload differs from the frozen Quick Bilayer contract")
            audit.require(request_record.get("dry_run") is False, "a dry-run PPM request cannot validate the build")
            audit.require(request_record.get("quick_bilayer_jobid") == construction.get("quick_bilayer_jobid"), "PPM response job ID differs from the build")
            response = request_record.get("response", {})
            audit.require(isinstance(response, dict) and str(response.get("jobid")) == construction.get("quick_bilayer_jobid"), "PPM response payload does not bind the Quick Bilayer job")
            client_path = package_root / "scripts" / "charmmgui_api.py"
            audit.require(request_record.get("client_sha256") == sha256(client_path), "PPM request record was created by another or changed API client")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot validate PPM request provenance: {exc}")
    for key in ("curator", "independent_reviewer", "approved_at_utc"):
        audit.require(isinstance(record.get(key), str) and bool(record[key].strip()) and not contains_todo(record[key]), f"orientation record {key} is unresolved")


def validate_environment_report(
    package_root: Path, artifact: Any, manifest: dict[str, Any], audit: Audit
) -> None:
    path = validate_artifact(package_root, artifact, "server environment validation report", audit)
    if path is None:
        return
    try:
        report = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"cannot read environment validation report: {exc}")
        return
    audit.require(report.get("schema_version") == "1.0" and report.get("report_type") == "server_environment_validation", "environment report schema/type is invalid")
    audit.require(not contains_todo(report), "environment report contains unresolved TODO/TBD values")
    audit.require(report.get("status") == "pass" and report.get("approval_status") == "approved", "environment report is not approved PASS")
    audit.require(report.get("study_id") == manifest.get("study_id"), "environment report study_id differs")
    audit.require(report.get("system_id") == PRIMARY_ID and report.get("construction_id") == "build01", "environment report system/build identity differs")
    audit.require(report.get("manifest_release_contract_sha256") == release_manifest_contract_sha256(manifest), "environment report is stale relative to release manifest")
    software = manifest.get("software", {})
    expected_gromacs = manifest.get("simulation", {}).get("required_version")
    audit.require(report.get("gromacs_version") == expected_gromacs, "environment GROMACS version differs from manifest")
    audit.require(report.get("container_digest") == software.get("container_digest"), "environment container digest differs from manifest")
    gmx_identity_path = validate_artifact(
        package_root, report.get("gmx_executable"), "environment GROMACS executable identity", audit
    )
    if gmx_identity_path is not None:
        try:
            identity = load_json(gmx_identity_path)
            audit.require(
                identity.get("schema_version") == "1.0"
                and identity.get("record_type") == "gromacs_executable_identity",
                "GROMACS executable identity schema/type is invalid",
            )
            audit.require(identity.get("gromacs_version") == expected_gromacs, "GROMACS identity version differs")
            audit.require(
                isinstance(identity.get("resolved_path"), str) and Path(identity["resolved_path"]).is_absolute(),
                "GROMACS identity path is not absolute",
            )
            audit.require(
                isinstance(identity.get("bytes"), int) and identity.get("bytes") > 0
                and isinstance(identity.get("sha256"), str) and bool(SHA256_RE.fullmatch(identity["sha256"])),
                "GROMACS identity binary hash/size is invalid",
            )
            linked = identity.get("linked_libraries")
            audit.require(isinstance(linked, list) and bool(linked), "GROMACS identity lacks linked libraries")
            for library in linked if isinstance(linked, list) else []:
                audit.require(
                    isinstance(library, dict)
                    and isinstance(library.get("path"), str) and Path(library["path"]).is_absolute()
                    and isinstance(library.get("bytes"), int) and library.get("bytes") > 0
                    and isinstance(library.get("sha256"), str) and bool(SHA256_RE.fullmatch(library["sha256"])),
                    "GROMACS linked-library identity is invalid",
                )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read GROMACS executable identity: {exc}")
    for key in ("gmx_dump_test", "gpu_or_cpu_runtime_test", "checkpoint_restart_test"):
        validate_artifact(package_root, report.get(key), f"environment {key}", audit)
    integrity = report.get("integrity", {})
    audit.require(
        isinstance(integrity, dict) and integrity.get("payload_sha256") == report_payload_sha256(report, ("integrity", "payload_sha256")),
        "environment report payload checksum is invalid",
    )


def validate_md_stage_report(
    package_root: Path, artifact: Any, expected_phase: str, manifest: dict[str, Any], audit: Audit
) -> None:
    path = validate_artifact(package_root, artifact, f"strict {expected_phase} validation report", audit)
    if path is None:
        return
    try:
        report = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"cannot read {expected_phase} validation report: {exc}")
        return
    audit.require(report.get("schema_version") == "2.0" and report.get("report_type") == "md_stage_output_validation", f"{expected_phase} report schema/type is invalid")
    audit.require(not contains_todo(report), f"{expected_phase} report contains unresolved TODO/TBD values")
    audit.require(report.get("phase") == expected_phase and report.get("strict") is True and report.get("status") == "pass", f"{expected_phase} report is not a strict PASS")
    audit.require(report.get("system_id") == PRIMARY_ID and report.get("construction_id") == "build01", f"{expected_phase} report system/build differs")
    audit.require(report.get("build_contract_sha256") == build_contract_sha256(manifest), f"{expected_phase} report build contract is stale")
    audit.require(report.get("manifest_release_contract_sha256") == release_manifest_contract_sha256(manifest), f"{expected_phase} report release-manifest contract is stale")
    construction = manifest["systems"][0]["construction"]
    audit.require(report.get("construction_archive_sha256") == construction["charmm_gui_archive"]["sha256"], f"{expected_phase} report archive hash differs")
    runs = report.get("runs")
    audit.require(isinstance(runs, list) and [item.get("realization_id") for item in runs] == list(REALIZATION_IDS), f"{expected_phase} report must contain rep01-rep03 in order")
    if isinstance(runs, list):
        audit.require(all(item.get("status") == "pass" for item in runs), f"{expected_phase} report contains a failed realization")
    if expected_phase == "canary":
        audit.require(report.get("analysis_disposition") == "eligible_for_checkpoint_continuation", "canary report does not authorize continuation")
        storage = report.get("storage_release_check", {})
        audit.require(isinstance(storage, dict) and storage.get("status") == "pass", "canary storage/output-rate release check did not pass")
        for item in runs if isinstance(runs, list) else []:
            rid = item.get("realization_id")
            audit.require(item.get("final_time_ps") == 5000.0, f"{rid} canary did not end exactly at 5 ns")
            tpr = item.get("artifacts", {}).get("tpr", {})
            tpr_record = item.get("production_tpr_record", {}) or {}
            audit.require(tpr.get("sha256") == tpr_record.get("production_tpr_sha256"), f"{rid} canary TPR hash is not frozen")
            matching_realizations = [
                value for value in manifest["systems"][0].get("realizations", [])
                if isinstance(value, dict) and value.get("id") == rid
            ]
            audit.require(len(matching_realizations) == 1, f"{rid} canary has no unique manifest realization")
            if len(matching_realizations) == 1:
                run_dir = (package_root / str(matching_realizations[0].get("run_directory", ""))).resolve()
                checkpoint = run_dir / "work" / "production.cpt"
                checkpoint_record = item.get("artifacts", {}).get("checkpoint", {})
                expected_checkpoint_sha256 = (
                    checkpoint_record.get("sha256") if isinstance(checkpoint_record, dict) else None
                )
                previous_finished_sha256: str | None = None
                command_logs = run_dir / "command_logs"
                for finished_path in sorted(command_logs.glob("[0-9][0-9][0-9][0-9]_*.finished.json")):
                    try:
                        finished = load_json(finished_path)
                        started_path = Path(str(finished.get("started_record", ""))).resolve()
                        audit.require(
                            started_path.is_file() and finished.get("started_record_sha256") == sha256(started_path),
                            f"{rid} command provenance is invalid before continuation: {finished_path}",
                        )
                        if not started_path.is_file():
                            continue
                        started = load_json(started_path)
                        audit.require(
                            started.get("previous_finished_record_sha256") == previous_finished_sha256,
                            f"{rid} command hash chain is broken before continuation: {finished_path}",
                        )
                        previous_finished_sha256 = sha256(finished_path)
                        argv = started.get("argv")
                        is_continuation = (
                            isinstance(argv, list) and len(argv) > 1 and argv[1] == "mdrun"
                            and "production.tpr" in argv and "-cpi" in argv and "-append" in argv
                            and "-nsteps" not in argv
                        )
                        if not is_continuation:
                            continue
                        audit.require(finished.get("returncode") == 0, f"{rid} has a failed production continuation")
                        checkpoint_input = started.get("checkpoint_input", {})
                        audit.require(
                            isinstance(checkpoint_input, dict)
                            and checkpoint_input.get("sha256") == expected_checkpoint_sha256,
                            f"{rid} production checkpoint input chain is broken",
                        )
                        checkpoint_output = finished.get("checkpoint_output", {})
                        next_sha256 = checkpoint_output.get("sha256") if isinstance(checkpoint_output, dict) else None
                        audit.require(
                            isinstance(next_sha256, str) and bool(SHA256_RE.fullmatch(next_sha256)),
                            f"{rid} production continuation lacks a sealed checkpoint output",
                        )
                        if isinstance(next_sha256, str) and SHA256_RE.fullmatch(next_sha256):
                            expected_checkpoint_sha256 = next_sha256
                    except (OSError, ValueError, json.JSONDecodeError) as exc:
                        audit.errors.append(f"{rid} cannot validate continuation checkpoint chain: {exc}")
                audit.require(
                    checkpoint.is_file() and checkpoint.stat().st_size > 0,
                    f"{rid} current production checkpoint is missing or empty",
                )
                if checkpoint.is_file() and checkpoint.stat().st_size > 0:
                    audit.require(
                        expected_checkpoint_sha256 == sha256(checkpoint),
                        f"{rid} current production checkpoint differs from its latest approved chain state",
                    )
    integrity = report.get("integrity", {})
    audit.require(
        isinstance(integrity, dict) and integrity.get("payload_sha256") == report_payload_sha256(report, ("integrity", "payload_sha256")),
        f"{expected_phase} report payload checksum is invalid",
    )


def validate_build_report_binding(
    report: dict[str, Any], construction: dict[str, Any], manifest: dict[str, Any], audit: Audit
) -> None:
    """Validate every immutable build-report identity and payload binding."""
    audit.require(report.get("schema_version") == "2.0" and report.get("report_type") == "charmm_gui_build_validation", "build01 validation report schema/type is invalid")
    audit.require(not contains_todo(report), "build01 validation report contains unresolved TODO/TBD values")
    audit.require(report.get("status") == "pass" and report.get("strict") is True, "build01 validation did not strict-pass")
    audit.require(report.get("study_id") == manifest.get("study_id"), "build report study_id differs")
    audit.require(report.get("system_id") == PRIMARY_ID and report.get("construction_id") == "build01", "build report system/build differs")
    audit.require(report.get("pdb_reader_jobid") == construction.get("pdb_reader_jobid"), "build report PDB Reader job differs")
    audit.require(report.get("quick_bilayer_jobid") == construction.get("quick_bilayer_jobid"), "build report Quick Bilayer job differs")
    audit.require(report.get("archive_sha256") == construction.get("charmm_gui_archive", {}).get("sha256"), "build report archive SHA differs")
    audit.require(report.get("build_contract_sha256") == build_contract_sha256(manifest), "build report is stale relative to build contract")
    expected_hmr = manifest.get("simulation", {}).get("hydrogen_mass_repartitioning")
    audit.require(
        report.get("hydrogen_mass_repartitioning_detected") is expected_hmr,
        "build report HMR state differs from the frozen production protocol",
    )
    audit.require(report.get("equilibration_mdp_sha256") == construction.get("equilibration_mdp_sha256"), "build report staged-equilibration MDP hashes differ")
    audit.require(report.get("gromacs_input_tree_manifest_sha256") == construction.get("gromacs_input_tree_manifest", {}).get("sha256"), "build report GROMACS tree hash differs")
    tree_binding = report.get("archive_extracted_gromacs_tree_binding")
    audit.require(isinstance(tree_binding, list) and bool(tree_binding), "build report lacks archive-to-extracted GROMACS tree binding")
    if isinstance(tree_binding, list):
        audit.require(
            all(
                isinstance(item, dict) and item.get("match") is True
                and item.get("archive_sha256") == item.get("extracted_sha256")
                and item.get("archive_bytes") == item.get("extracted_bytes")
                for item in tree_binding
            ),
            "build report contains an archive/extracted GROMACS file mismatch",
        )
        relative_paths = [item.get("relative_path") for item in tree_binding if isinstance(item, dict)]
        audit.require(len(relative_paths) == len(set(relative_paths)), "archive/extracted GROMACS binding has duplicate paths")
    stage_binding = report.get("staged_mdp_archive_binding_and_physics")
    expected_stage_names = ["step6.0_minimization.mdp"] + [
        f"step6.{index}_equilibration.mdp" for index in range(1, 7)
    ]
    audit.require(
        isinstance(stage_binding, list)
        and [item.get("name") for item in stage_binding if isinstance(item, dict)] == expected_stage_names,
        "build report lacks the exact seven ordered CHARMM-GUI staged MDP bindings",
    )
    if isinstance(stage_binding, list):
        for item in stage_binding:
            if not isinstance(item, dict):
                continue
            name = item.get("name")
            expected_hash = (
                construction.get("minimization_mdp", {}).get("sha256")
                if name == "step6.0_minimization.mdp"
                else construction.get("equilibration_mdp_sha256", {}).get(name)
            )
            audit.require(
                item.get("archive_sha256") == expected_hash
                and item.get("manifest_sha256") == expected_hash
                and isinstance(item.get("physics"), dict),
                f"build report staged MDP archive/manifest/physics binding failed for {name}",
            )
    for report_key, construction_key in (
        ("membrane_orientation_record_sha256", "membrane_orientation_record"),
        ("starting_coordinates_sha256", "starting_coordinates"),
        ("topology_sha256", "topology"),
        ("index_sha256", "index"),
        ("analysis_index_sha256", "analysis_index"),
        ("production_mdp_sha256", "production_mdp"),
        ("minimization_mdp_sha256", "minimization_mdp"),
    ):
        audit.require(report.get(report_key) == construction.get(construction_key, {}).get("sha256"), f"build report {report_key} differs")
    integrity = report.get("integrity", {})
    audit.require(isinstance(integrity, dict) and integrity.get("payload_sha256") == report_payload_sha256(report, ("integrity", "payload_sha256")), "build report payload checksum is invalid")


def validate_structure_record(
    package_root: Path, record_path: str, role: str, audit: Audit
) -> dict[str, Any] | None:
    path = resolve_inside(package_root, record_path, f"{role} structure record", audit)
    if path is None:
        return None
    try:
        record = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"Cannot read structure record {path}: {exc}")
        return None
    audit.require(record.get("schema_version") == "2.0", f"{role} structure record schema must be 2.0")
    audit.require(not contains_todo(record), f"{role} structure record contains unresolved TODO/TBD values")
    audit.require(record.get("approval_status") == "approved", f"{role} structure record is not approved")
    audit.require(record.get("system_id") == PRIMARY_ID, f"{role} structure record system ID is wrong")
    audit.require(str(record.get("pdb_id", "")).upper() == "8KCT", f"{role} must derive from 8KCT")
    audit.require(record.get("biological_assembly") == 1, f"{role} must use biological assembly 1")
    audit.require(record.get("chain_mapping_verified") is True, f"{role} chain mapping is not verified")
    audit.require(record.get("ligand_stereochemistry_verified") is True, f"{role} O6U stereochemistry is not verified")
    audit.require(record.get("psen1_292_376_not_modelled") is True, f"{role} must not invent unresolved PSEN1 residues 292-376")
    audit.require(record.get("psen1_ntf_ctf_separate_topology_segments") is True, f"{role} PSEN1 NTF/CTF topology policy is not verified")
    audit.require(record.get("glycosylation_sites_verified") is True, f"{role} glycosylation sites are not verified")
    audit.require(record.get("glycan_covalent_links_verified") is True, f"{role} glycan covalent links are not verified")
    audit.require(record.get("covalent_link_counts") == {
        "disulfides": 4, "protein_to_nag": 12, "glycan_internal": 9, "total_struct_conn": 25,
    }, f"{role} deposited covalent-link counts are incomplete or changed")
    counts = record.get("component_counts")
    audit.require(isinstance(counts, dict), f"{role} component_counts must be an object")
    if isinstance(counts, dict):
        audit.require(counts.get("O6U") == 1, f"{role} must contain one O6U")
        audit.require(counts.get("CLR") == 3, f"{role} must retain three resolved CLR molecules")
        audit.require(counts.get("PC1") == 2, f"{role} must retain two resolved PC1 molecules")
        audit.require(counts.get("NAG") == 18, f"{role} must retain exactly 18 resolved NAG residues")
        audit.require(counts.get("BMA") == 3, f"{role} must retain exactly 3 resolved BMA residues")
    audit.require(record.get("resolved_segments") == {
        "NCSTN_A": "34-700", "PSEN1_B_NTF": "76-291", "PSEN1_B_CTF": "377-467",
        "APH1A_C": "2-244", "PEN2_D": "6-101",
    }, f"{role} resolved segment mapping is incomplete or changed")
    audit.require(record.get("disulfides") == [
        "NCSTN Cys50-Cys62", "NCSTN Cys140-Cys159", "NCSTN Cys230-Cys248", "NCSTN Cys586-Cys620",
    ], f"{role} disulfide list is incomplete or changed")
    audit.require(record.get("glycosylation_sites") == [45, 55, 187, 264, 387, 435, 464, 506, 530, 562, 573, 580], f"{role} glycosylation-site list is incomplete or changed")
    audit.require(isinstance(record.get("termini_and_patches"), list) and bool(record["termini_and_patches"]) and not contains_todo(record["termini_and_patches"]), f"{role} termini/patch decisions are unresolved")
    audit.require(isinstance(record.get("native_component_dispositions"), list) and bool(record["native_component_dispositions"]) and not contains_todo(record["native_component_dispositions"]), f"{role} native-component dispositions are unresolved")
    audit.require(record.get("ligand_formal_charge_e") == 0, f"{role} O6U formal charge must be zero")
    audit.require(bool(SHA256_RE.fullmatch(str(record.get("ligand_atom_mapping_sha256", "")).lower())), f"{role} ligand atom-mapping SHA-256 is unresolved")
    source_coordinates = validate_artifact(package_root, record.get("source_coordinates"), f"{role} source coordinates", audit)
    validate_artifact(package_root, record.get("curated_coordinates"), f"{role} curated coordinates", audit)
    validate_artifact(package_root, record.get("wwpdb_validation_report"), f"{role} validation report", audit)
    source_audit_path = validate_artifact(package_root, record.get("deposited_structure_audit"), f"{role} deposited-structure audit", audit)
    validate_artifact(package_root, record.get("pdb_reader_review_checklist"), f"{role} PDB Reader review checklist", audit)
    if source_coordinates is not None:
        audit.require(sha256(source_coordinates) == KCT8_CIF_SHA256, f"{role} source coordinate hash is not the pinned 8KCT mmCIF")
    if source_audit_path is not None:
        try:
            source_audit = load_json(source_audit_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read deposited-structure audit: {exc}")
        else:
            audit.require(source_audit.get("overall_status") == "deposited_structure_pass_manual_pdb_reader_blocked", "deposited-structure audit status is invalid")
            audit.require(source_audit.get("pdb_reader_approved") is False and source_audit.get("quick_bilayer_submission_allowed") is False, "deposited-structure audit must remain a non-approving source audit")
            audit.require(source_audit.get("source", {}).get("sha256") == KCT8_CIF_SHA256, "deposited-structure audit is not bound to the pinned 8KCT mmCIF")
            audit.require(source_audit.get("deposited_model", {}).get("native_component_counts") == {
                "BMA": 3, "CLR": 3, "NAG": 18, "O6U": 1, "PC1": 2,
            }, "deposited-structure audit component counts differ")
            link_record = source_audit.get("deposited_covalent_connections", {})
            audit.require(
                link_record.get("disulfide_count") == 4
                and link_record.get("protein_to_nag_count") == 12
                and link_record.get("glycan_internal_count") == 9
                and link_record.get("total_struct_conn_records") == 25,
                "deposited-structure audit covalent-link counts differ",
            )
    pdb_reader_evidence = record.get("pdb_reader_evidence")
    audit.require(isinstance(pdb_reader_evidence, list) and bool(pdb_reader_evidence), f"{role} PDB Reader evidence is missing")
    if isinstance(pdb_reader_evidence, list):
        for index, artifact in enumerate(pdb_reader_evidence):
            validate_artifact(package_root, artifact, f"{role} PDB Reader evidence {index}", audit)
    resolve_inside(package_root, str(record.get("ligand_parameter_record", "")), f"{role} ligand parameter record", audit)
    resolve_inside(package_root, str(record.get("protein_protonation_record", "")), f"{role} protein protonation record", audit)
    for key in ("pdb_reader_jobid", "curator", "independent_reviewer", "approved_at_utc"):
        value = record.get(key)
        audit.require(isinstance(value, str) and bool(value.strip()) and not contains_todo(value), f"{role}.{key} is unresolved")
    return record


def validate_realization_template(package_root: Path, audit: Audit) -> None:
    """Check the live post-stage realization template without requiring future outputs."""
    path = package_root / "templates" / "realization_record.template.json"
    if not path.is_file():
        audit.errors.append(f"realization record template is missing: {path}")
        return
    try:
        record = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"cannot read realization record template: {exc}")
        return
    audit.require(record.get("schema_version") == "2.0", "realization template schema must be 2.0")
    audit.require(record.get("system_id") == PRIMARY_ID, "realization template system ID is wrong")
    audit.require(record.get("common_build_id") == "build01", "realization template must reference build01")
    audit.require(record.get("velocity_generation_stage") == "first_NVT", "realization velocities must be generated at first NVT")
    audit.require(record.get("independent_dynamic_equilibration") is True, "each realization requires independent dynamic equilibration")
    audit.require(record.get("analysis_window_ns") == [200, 500], "realization template analysis window must be 200-500 ns")
    canary = record.get("canary_contract", {})
    audit.require(isinstance(canary, dict), "realization template lacks canary_contract")
    if isinstance(canary, dict):
        audit.require(canary.get("target_ns") == 5.0 and canary.get("target_steps_at_2fs") == 2_500_000 and canary.get("expected_endpoint_ps") == 5000.0, "realization canary endpoint/step contract is wrong")
        audit.require(canary.get("uses_original_500ns_tpr") is True and canary.get("frames_retained_as_start_of_production") is True, "realization canary must retain the original 500-ns TPR and frames")
        audit.require(canary.get("ligand_behavior_used_for_selection") is False, "realization canary cannot select by ligand behavior")
        for key in ("production_tpr_record", "checkpoint_at_5ns", "all_three_canary_validation_report"):
            artifact = canary.get(key)
            audit.require(isinstance(artifact, dict) and set(artifact) == {"path", "sha256"}, f"realization canary {key} artifact schema is missing")
    audit.require(
        record.get("failure_policy") == "inconclusive_if_any_realization_fails_qc_or_stationarity",
        "realization template must fail closed if any realization fails QC or stationarity",
    )
    forbidden_template_keys = [key for key in record if "recovery" in key.lower() or "extension" in key.lower()]
    audit.require(not forbidden_template_keys, f"realization template retains obsolete keys: {forbidden_template_keys}")


def validate_parameter_record(package_root: Path, record_path: str, audit: Audit) -> dict[str, Any] | None:
    path = resolve_inside(package_root, record_path, "nirogacestat parameter record", audit)
    if path is None:
        return None
    try:
        record = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        audit.errors.append(f"Cannot read ligand parameter record {path}: {exc}")
        return None

    audit.require(record.get("schema_version") == "1.1", "nirogacestat parameter record schema must be 1.1")
    audit.require(record.get("approval_status") == "approved", "nirogacestat parameters are not approved")
    audit.require(not contains_todo(record), "nirogacestat parameter record contains unresolved TODO/TBD values")
    audit.require(record.get("study_id") == "gamma_secretase_native_nirogacestat_rebuild_20260808", "O6U parameter record study_id is wrong")
    audit.require(str(record.get("component_id", "")).upper() == "O6U", "parameter record component must be O6U")
    audit.require(record.get("atom_mapping_complete") is True, "O6U atom mapping is incomplete")
    audit.require(record.get("stereochemistry_verified") is True, "O6U stereochemistry is not verified")
    preaudit_path = validate_artifact(package_root, record.get("preparameterization_audit"), "O6U preparameterization audit", audit)
    official_sdf_path = validate_artifact(package_root, record.get("official_source_sdf"), "official neutral O6U source SDF", audit)
    parameter_sdf_path = validate_artifact(package_root, record.get("parameterization_input_sdf"), "normalized neutral O6U parameterization SDF", audit)
    premapping_path = validate_artifact(package_root, record.get("preparameterization_atom_correspondence"), "O6U preparameterization atom correspondence", audit)
    parameterization_toolchain_path = validate_artifact(
        package_root,
        record.get("parameterization_toolchain_record"),
        "O6U parameterization toolchain record",
        audit,
    )
    parameterization_toolchain_validation_path = validate_artifact(
        package_root,
        record.get("parameterization_toolchain_independent_validation"),
        "O6U parameterization toolchain independent validation",
        audit,
    )
    water_precision_path = validate_artifact(
        package_root,
        record.get("water_probe_coordinate_precision_audit"),
        "O6U water-probe coordinate-precision audit",
        audit,
    )
    water_policy_path = validate_artifact(
        package_root,
        record.get("water_probe_disposition_policy"),
        "O6U water-probe disposition policy",
        audit,
    )
    water_chemical_prescreen_path = validate_artifact(
        package_root,
        record.get("water_probe_chemical_role_prescreen"),
        "O6U water-probe chemical-role prescreen",
        audit,
    )
    if official_sdf_path is not None:
        audit.require(sha256(official_sdf_path) == O6U_OFFICIAL_SDF_SHA256, "official neutral O6U source SDF hash differs")
    if parameter_sdf_path is not None:
        audit.require(sha256(parameter_sdf_path) == O6U_PARAMETER_INPUT_SDF_SHA256, "normalized neutral O6U parameterization SDF hash differs")
    if premapping_path is not None:
        audit.require(sha256(premapping_path) == O6U_PREMAPPING_SHA256, "O6U preparameterization atom correspondence hash differs")
    if parameterization_toolchain_path is not None:
        audit.require(
            sha256(parameterization_toolchain_path) == O6U_PARAMETERIZATION_TOOLCHAIN_RECORD_SHA256,
            "O6U parameterization toolchain record hash differs",
        )
    if parameterization_toolchain_validation_path is not None:
        audit.require(
            sha256(parameterization_toolchain_validation_path) == O6U_PARAMETERIZATION_TOOLCHAIN_VALIDATION_SHA256,
            "O6U parameterization toolchain independent-validation hash differs",
        )
        try:
            toolchain_validation = load_json(parameterization_toolchain_validation_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read O6U parameterization toolchain validation: {exc}")
        else:
            audit.require(
                toolchain_validation.get("status") == "pass_independently_reconstructed"
                and toolchain_validation.get("production_approved") is False,
                "O6U parameterization toolchain validation is not an independent PASS",
            )
            audit.require(
                toolchain_validation.get("source_record", {}).get("sha256")
                == O6U_PARAMETERIZATION_TOOLCHAIN_RECORD_SHA256,
                "O6U toolchain validation does not bind the frozen toolchain record",
            )
            audit.require(
                toolchain_validation.get("ffparam_archive", {}).get("sha256")
                == "d9508f3a1590ba9fbfb1d048e832d6726ef6131fd40a70572f5662c5bdc2cbdb"
                and toolchain_validation.get("ffparam_source_tree", {}).get("file_count") == 262,
                "O6U toolchain validation does not bind the official FFParam 1.2.0 distribution",
            )
            audit.require(
                toolchain_validation.get("python_modules", {}).get("psi4", {}).get("version_attribute") == "1.9.1",
                "O6U toolchain validation Psi4 version differs",
            )
    if water_precision_path is not None:
        audit.require(
            sha256(water_precision_path) == O6U_WATER_PRECISION_AUDIT_SHA256,
            "O6U water-probe coordinate-precision audit hash differs",
        )
        try:
            water_precision = load_json(water_precision_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read O6U water-probe coordinate-precision audit: {exc}")
        else:
            audit.require(
                water_precision.get("status") == "pass_crd_high_precision_formal_pdb_rounding_explained"
                and water_precision.get("production_approved") is False,
                "O6U water-probe coordinate-precision audit is not PASS",
            )
            audit.require(
                water_precision.get("coordinate_atom_count") == 76
                and water_precision.get("orientation_count") == 70
                and water_precision.get("orientation_identity_fields_identical") is True,
                "O6U water-probe coordinate-precision invariants differ",
            )
            audit.require(
                water_precision.get("formal_frozen_orientation_source", {}).get("orientation_da_sha256")
                == O6U_WATER_ORIENTATION_DA_SHA256,
                "O6U water-probe precision audit does not freeze the CRD-derived plan",
            )
    if water_policy_path is not None:
        audit.require(
            sha256(water_policy_path) == O6U_WATER_DISPOSITION_POLICY_SHA256,
            "O6U water-probe disposition policy hash differs",
        )
        try:
            water_policy = load_json(water_policy_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read O6U water-probe disposition policy: {exc}")
        else:
            audit.require(
                water_policy.get("schema_version") == "1.0"
                and water_policy.get("report_type") == "o6u_water_probe_disposition_policy"
                and water_policy.get("status") == "pass"
                and water_policy.get("production_approved") is False,
                "O6U water-probe disposition policy is invalid",
            )
            audit.require(
                water_policy.get("allowed_qm_status_by_disposition") == {
                    "applicable": ["pass"],
                    "weak": ["pass"],
                    "unfavourable": ["pass"],
                    "excluded": ["not_required_prespecified_exclusion"],
                },
                "O6U water-probe disposition policy changed its fail-closed QM rules",
            )
            policy_sources = water_policy.get("sources")
            audit.require(
                isinstance(policy_sources, list) and len(policy_sources) == 3,
                "O6U water-probe disposition policy lacks its three frozen sources",
            )
            if isinstance(policy_sources, list):
                for source_index, source_record in enumerate(policy_sources):
                    validate_artifact(
                        package_root,
                        source_record,
                        f"O6U water-probe policy source {source_index}",
                        audit,
                    )
    if water_chemical_prescreen_path is not None:
        audit.require(
            sha256(water_chemical_prescreen_path) == O6U_WATER_CHEMICAL_PRESCREEN_SHA256,
            "O6U water-probe chemical-role prescreen hash differs",
        )
        try:
            water_prescreen = load_json(water_chemical_prescreen_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read O6U water-probe chemical-role prescreen: {exc}")
        else:
            audit.require(
                water_prescreen.get("schema_version") == "1.0"
                and water_prescreen.get("report_type") == "o6u_water_probe_chemical_role_prescreen"
                and water_prescreen.get("status") == "pass_chemical_role_prescreen_visual_review_required"
                and water_prescreen.get("production_approved") is False,
                "O6U water-probe chemical-role prescreen is invalid",
            )
            audit.require(
                water_prescreen.get("orientation_count") == 70
                and water_prescreen.get("retained_for_visual_review_count") == 20
                and water_prescreen.get("chemically_excludable_count") == 50,
                "O6U water-probe chemical-role prescreen counts differ",
            )
            prescreen_inputs = water_prescreen.get("inputs", {})
            audit.require(
                prescreen_inputs.get("orientation_da", {}).get("sha256") == O6U_WATER_ORIENTATION_DA_SHA256
                and prescreen_inputs.get("sdf", {}).get("sha256") == O6U_PARAMETER_INPUT_SDF_SHA256
                and prescreen_inputs.get("correspondence_tsv", {}).get("sha256") == O6U_CGENFF_CORRESPONDENCE_TSV_SHA256
                and prescreen_inputs.get("policy", {}).get("sha256") == O6U_WATER_DISPOSITION_POLICY_SHA256,
                "O6U water-probe chemical-role prescreen does not bind the frozen inputs",
            )
            prescreen_rows = water_prescreen.get("orientations")
            audit.require(
                isinstance(prescreen_rows, list)
                and len(prescreen_rows) == 70
                and {item.get("orientation_id") for item in prescreen_rows if isinstance(item, dict)}
                == {f"O6U_WP_{index:03d}" for index in range(1, 71)},
                "O6U water-probe chemical-role prescreen row coverage differs",
            )
    if preaudit_path is not None:
        try:
            preaudit = load_json(preaudit_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read O6U preparameterization audit: {exc}")
        else:
            audit.require(preaudit.get("overall_status") == "local_identity_pass_parameterization_blocked", "O6U preparameterization audit status is invalid")
            audit.require(preaudit.get("md_parameterization_approved") is False, "O6U identity audit must not claim parameter approval")
            audit.require(preaudit.get("identity") == {
                "component_id": "O6U", "formal_charge": 0, "formula": "C27H41F2N5O",
                "formula_weight": 489.644, "name": "Nirogacestat",
            }, "O6U preparameterization identity differs")
            sdf_record = preaudit.get("sdf_validation", {})
            audit.require(
                sdf_record.get("atom_count") == 76
                and sdf_record.get("heavy_atom_count") == 35
                and sdf_record.get("hydrogen_count") == 41
                and sdf_record.get("bond_count") == 78
                and sdf_record.get("formal_charge") == 0
                and sdf_record.get("formula") == "C27H41F2N5O"
                and sdf_record.get("assigned_chiral_centres_zero_based") == [[10, "S"], [13, "S"]],
                "O6U preparameterization SDF invariants differ",
            )
            inputs = preaudit.get("inputs", {})
            outputs = preaudit.get("outputs", {})
            audit.require(inputs.get("ccd", {}).get("sha256") == O6U_CCD_SHA256, "O6U preparameterization audit CCD hash differs")
            audit.require(inputs.get("hydrogen_complete_sdf", {}).get("sha256") == O6U_OFFICIAL_SDF_SHA256, "O6U preparameterization audit source SDF hash differs")
            audit.require(inputs.get("native_structure", {}).get("sha256") == KCT8_CIF_SHA256, "O6U preparameterization audit 8KCT hash differs")
            parameter_input = outputs.get("parameterization_input_sdf", {})
            audit.require(
                parameter_input.get("sha256") == O6U_PARAMETER_INPUT_SDF_SHA256
                and parameter_input.get("explicit_3d") is True
                and parameter_input.get("maximum_coordinate_change_from_official_sdf_angstrom") == 0.0,
                "O6U normalized 3D parameterization input differs",
            )
            audit.require(outputs.get("atom_correspondence", {}).get("sha256") == O6U_PREMAPPING_SHA256 and outputs.get("atom_correspondence", {}).get("rows") == 76, "O6U preparameterization correspondence output differs")
    formal = finite_number(record.get("formal_charge_e"), "O6U formal charge", audit)
    topology = finite_number(record.get("topology_charge_sum_e"), "O6U topology charge sum", audit)
    if formal is not None and topology is not None:
        audit.require(abs(formal - topology) <= 0.0001, "O6U topology charge does not match approved formal charge")

    charges = record.get("atom_charge_table")
    audit.require(isinstance(charges, list) and len(charges) == 76, "O6U atom-charge table must contain exactly 76 named atoms")
    if isinstance(charges, list) and charges:
        parsed: list[float] = []
        pdb_names: list[str] = []
        gromacs_names: list[str] = []
        for index, item in enumerate(charges, start=1):
            if not isinstance(item, dict):
                audit.errors.append(f"O6U atom-charge row {index} must be an object")
                continue
            audit.require(item.get("index") == index, f"O6U atom-charge row {index} has wrong one-based index")
            pdb_name = item.get("pdb_atom_name")
            gromacs_name = item.get("gromacs_atom_name")
            audit.require(isinstance(pdb_name, str) and bool(pdb_name.strip()), f"O6U row {index} lacks pdb_atom_name")
            audit.require(isinstance(gromacs_name, str) and bool(gromacs_name.strip()), f"O6U row {index} lacks gromacs_atom_name")
            if isinstance(pdb_name, str):
                pdb_names.append(pdb_name)
            if isinstance(gromacs_name, str):
                gromacs_names.append(gromacs_name)
            parsed_value = finite_number(item.get("partial_charge_e"), f"O6U partial charge {index}", audit)
            if parsed_value is not None:
                parsed.append(parsed_value)
        audit.require(len(pdb_names) == 76 and len(set(pdb_names)) == 76, "O6U PDB atom names must be 76 unique values")
        audit.require(len(gromacs_names) == 76 and len(set(gromacs_names)) == 76, "O6U GROMACS atom names must be 76 unique values")
        if parsed:
            audit.require(any(abs(value) > 1e-8 for value in parsed), "all O6U partial charges are zero")
            if topology is not None:
                audit.require(abs(sum(parsed) - topology) <= 0.0001, "listed O6U charges do not sum to topology charge")

    max_parameter = finite_number(record.get("initial_max_parameter_penalty"), "maximum parameter penalty", audit)
    max_charge = finite_number(record.get("initial_max_charge_penalty"), "maximum charge penalty", audit)
    maximum = max(value for value in (max_parameter, max_charge) if value is not None) if any(
        value is not None for value in (max_parameter, max_charge)
    ) else None
    artifacts = record.get("validation_artifacts", [])
    if maximum is not None and maximum >= 10.0:
        audit.require(isinstance(artifacts, list) and bool(artifacts), "CGenFF penalty >=10 requires validation artifacts")
    if maximum is not None and maximum > 50.0:
        audit.require(record.get("validation_level") == "extensive", "CGenFF penalty >50 requires extensive validation")
    elif maximum is not None and maximum >= 10.0:
        audit.require(record.get("validation_level") in {"basic", "extensive"}, "CGenFF penalty 10-50 requires basic or extensive validation")
    if isinstance(artifacts, list):
        for index, artifact in enumerate(artifacts):
            validate_artifact(package_root, artifact, f"ligand validation artifact {index}", audit)
    unresolved = record.get("unresolved_high_penalty_terms")
    audit.require(unresolved == [], "unresolved high-penalty ligand terms remain")
    selection_path = validate_artifact(
        package_root,
        record.get("qm_representative_selection"),
        "frozen O6U QM representative selection",
        audit,
    )
    audit.require(
        record.get("all_frozen_qm_representatives_included") is True,
        "QM target set does not include every frozen O6U representative",
    )
    audit.require(record.get("qm_representative_count") == 5, "O6U QM representative count must be exactly five")
    if selection_path is not None:
        try:
            selection = load_json(selection_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"cannot read frozen O6U QM representative selection: {exc}")
        else:
            selected = selection.get("selected")
            selected_frames = (
                {
                    int(item.get("crest_frame_1based"))
                    for item in selected
                    if isinstance(item, dict) and isinstance(item.get("crest_frame_1based"), int)
                }
                if isinstance(selected, list)
                else set()
            )
            audit.require(selection.get("status") == "pass", "frozen O6U QM representative selection is not PASS")
            audit.require(selection.get("selected_count") == 5, "frozen O6U selection does not contain five representatives")
            method = selection.get("method")
            audit.require(
                isinstance(method, dict) and method.get("tfd_cutoff") == 0.2,
                "frozen O6U selection TFD cutoff differs from 0.20",
            )
            audit.require(
                selection.get("crest_global_minimum_frame_1based") == 1,
                "frozen O6U global-minimum frame differs from frame 1",
            )
            audit.require(
                selected_frames == {1, 342, 641, 679, 768},
                "frozen O6U representative frame identities differ",
            )
    audit.require(record.get("manual_residual_charge_spreading") is False, "manual residual-charge spreading is prohibited")
    audit.require(
        record.get("missing_duplicate_guessed_or_overwritten_parameters") is False,
        "missing, duplicate, guessed, or overwritten ligand parameters remain",
    )

    metrics = record.get("fit_metrics")
    if not isinstance(metrics, dict):
        audit.errors.append("ligand fit_metrics must be an object")
    else:
        def numeric_list(key: str) -> list[float]:
            raw = metrics.get(key)
            if not isinstance(raw, list) or not raw:
                audit.errors.append(f"fit_metrics.{key} must be a non-empty list")
                return []
            values: list[float] = []
            for index, item in enumerate(raw):
                value = finite_number(item, f"fit_metrics.{key}[{index}]", audit)
                if value is not None:
                    values.append(value)
            return values

        water_energy = numeric_list("water_interaction_energy_residuals_kcal_mol")
        if water_energy:
            audit.require(max(abs(value) for value in water_energy) <= 0.5, "a water-interaction energy residual exceeds 0.5 kcal/mol")
            audit.require(sum(abs(value) for value in water_energy) / len(water_energy) <= 0.2, "mean absolute water-interaction energy residual exceeds 0.2 kcal/mol")
        water_distance = numeric_list("water_interaction_distance_residuals_angstrom")
        if water_distance:
            audit.require(max(abs(value) for value in water_distance) <= 0.2, "a water-interaction distance residual exceeds 0.2 Angstrom")
        bonds = numeric_list("equilibrium_bond_deviations_angstrom")
        if bonds:
            audit.require(max(abs(value) for value in bonds) <= 0.03, "an equilibrium bond deviation exceeds 0.03 Angstrom")
        angles = numeric_list("equilibrium_angle_deviations_degrees")
        if angles:
            audit.require(max(abs(value) for value in angles) <= 3.0, "an equilibrium angle deviation exceeds 3 degrees")
        dipole_ratio = finite_number(metrics.get("dipole_magnitude_ratio_mm_to_qm"), "dipole magnitude ratio", audit)
        if dipole_ratio is not None:
            audit.require(1.20 <= dipole_ratio <= 1.50, "MM/QM dipole magnitude ratio must be 1.20-1.50")
        dipole_angle = finite_number(metrics.get("dipole_direction_error_degrees"), "dipole direction error", audit)
        if dipole_angle is not None:
            audit.require(
                metrics.get("dipole_direction_signed_review_passed") is True,
                "dipole direction lacks the required signed joint review with water-interaction targets",
            )
            validate_artifact(
                package_root,
                metrics.get("dipole_direction_review_artifact"),
                "dipole direction signed review",
                audit,
            )
        vibration = finite_number(
            metrics.get("vibrational_modes_below_2700_cm1_mean_absolute_relative_deviation"),
            "vibrational-mode mean absolute relative deviation",
            audit,
        )
        if vibration is not None:
            audit.require(vibration <= 0.05, "vibrational-mode mean absolute relative deviation exceeds 5%")
            audit.require(
                metrics.get("ped_low_frequency_signed_review_passed") is True,
                "the scalar vibrational average lacks the required signed low-frequency PED review",
            )
            validate_artifact(
                package_root,
                metrics.get("ped_low_frequency_review_artifact"),
                "low-frequency PED signed review",
                audit,
            )
        audit.require(
            metrics.get("water_orientation_dispositions_complete") is True,
            "applicable, excluded, weak, and unfavourable water orientations lack complete dispositions",
        )
        water_disposition_path = validate_artifact(
            package_root,
            metrics.get("water_orientation_disposition_artifact"),
            "water-orientation disposition report",
            audit,
        )
        if water_disposition_path is not None:
            try:
                water_dispositions = load_json(water_disposition_path)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                audit.errors.append(f"cannot read water-orientation disposition report: {exc}")
            else:
                audit.require(
                    water_dispositions.get("schema_version") == "1.1"
                    and water_dispositions.get("report_type") == "o6u_water_probe_prospective_disposition_table",
                    "water-orientation disposition report schema/type is invalid",
                )
                audit.require(
                    water_dispositions.get("status") == "complete_prospective_dispositions",
                    "water-orientation disposition report is not complete",
                )
                audit.require(
                    water_dispositions.get("production_approved") is False,
                    "water-orientation disposition report improperly claims production approval",
                )
                disposition_policy = water_dispositions.get("disposition_policy")
                audit.require(
                    isinstance(disposition_policy, dict)
                    and disposition_policy.get("sha256") == O6U_WATER_DISPOSITION_POLICY_SHA256,
                    "water-orientation disposition report does not bind the frozen policy",
                )
                disposition_policy_path = validate_artifact(
                    package_root,
                    disposition_policy,
                    "water-orientation disposition policy",
                    audit,
                )
                if disposition_policy_path is not None:
                    audit.require(
                        sha256(disposition_policy_path) == O6U_WATER_DISPOSITION_POLICY_SHA256,
                        "water-orientation disposition policy content differs",
                    )
                source_orientation = water_dispositions.get("source_orientation_da")
                audit.require(
                    isinstance(source_orientation, dict)
                    and source_orientation.get("sha256") == O6U_WATER_ORIENTATION_DA_SHA256,
                    "water-orientation disposition source hash differs from the frozen 70-orientation plan",
                )
                source_orientation_path = validate_artifact(
                    package_root,
                    source_orientation,
                    "frozen O6U water-orientation plan",
                    audit,
                )
                source_orientation_lines: list[str] = []
                if source_orientation_path is not None:
                    try:
                        source_orientation_lines = [
                            line.strip()
                            for line in source_orientation_path.read_text(encoding="utf-8").splitlines()
                            if line.strip()
                        ]
                    except OSError as exc:
                        audit.errors.append(f"cannot read frozen O6U water-orientation plan: {exc}")
                    else:
                        audit.require(
                            sha256(source_orientation_path) == O6U_WATER_ORIENTATION_DA_SHA256
                            and len(source_orientation_lines) == 70,
                            "frozen O6U water-orientation plan content differs",
                        )
                rows = water_dispositions.get("orientations")
                audit.require(
                    water_dispositions.get("orientation_count") == 70
                    and water_dispositions.get("pending_count") == 0
                    and isinstance(rows, list)
                    and len(rows) == 70,
                    "water-orientation disposition report must contain 70 completed rows",
                )
                if isinstance(rows, list):
                    expected_ids = {f"O6U_WP_{index:03d}" for index in range(1, 71)}
                    observed_ids: set[str] = set()
                    observed_lines: set[int] = set()
                    observed_definitions: set[str] = set()
                    observed_counts = {key: 0 for key in ("applicable", "excluded", "weak", "unfavourable")}
                    expected_probe_counts = {"A2": 2, "A31": 16, "AP": 2, "APL": 6, "D": 38, "DOP": 6}
                    observed_probe_counts = {key: 0 for key in expected_probe_counts}
                    for index, item in enumerate(rows, start=1):
                        if not isinstance(item, dict):
                            audit.errors.append(f"water-orientation row {index} must be an object")
                            continue
                        orientation_id = item.get("orientation_id")
                        line_number = item.get("source_line_number")
                        definition = item.get("source_definition")
                        disposition = item.get("prospective_disposition")
                        probe_type = item.get("probe_type")
                        audit.require(
                            isinstance(orientation_id, str) and orientation_id in expected_ids,
                            f"water-orientation row {index} has an invalid orientation_id",
                        )
                        if isinstance(orientation_id, str):
                            audit.require(orientation_id not in observed_ids, f"water orientation {orientation_id} is duplicated")
                            observed_ids.add(orientation_id)
                        audit.require(
                            isinstance(line_number, int) and 1 <= line_number <= 70,
                            f"water-orientation row {index} has an invalid source line number",
                        )
                        if isinstance(line_number, int):
                            audit.require(line_number not in observed_lines, f"water source line {line_number} is duplicated")
                            observed_lines.add(line_number)
                        audit.require(
                            isinstance(definition, str) and bool(definition.strip()),
                            f"water-orientation row {index} lacks its frozen source definition",
                        )
                        if isinstance(definition, str):
                            audit.require(definition not in observed_definitions, f"water source definition is duplicated at row {index}")
                            observed_definitions.add(definition)
                            if (
                                isinstance(line_number, int)
                                and 1 <= line_number <= len(source_orientation_lines)
                            ):
                                audit.require(
                                    definition == source_orientation_lines[line_number - 1],
                                    f"water-orientation row {index} differs from frozen source line {line_number}",
                                )
                        audit.require(
                            isinstance(disposition, str) and disposition in observed_counts,
                            f"water-orientation row {index} lacks a valid final disposition",
                        )
                        if isinstance(disposition, str) and disposition in observed_counts:
                            observed_counts[disposition] += 1
                        audit.require(
                            isinstance(probe_type, str) and probe_type in observed_probe_counts,
                            f"water-orientation row {index} has an invalid probe type",
                        )
                        if isinstance(probe_type, str) and probe_type in observed_probe_counts:
                            observed_probe_counts[probe_type] += 1
                        audit.require(
                            isinstance(item.get("disposition_rationale"), str)
                            and bool(item["disposition_rationale"].strip()),
                            f"water-orientation row {index} lacks a disposition rationale",
                        )
                        audit.require(
                            isinstance(item.get("reviewer"), str) and bool(item["reviewer"].strip()),
                            f"water-orientation row {index} lacks a reviewer",
                        )
                        audit.require(
                            isinstance(item.get("reviewed_at_utc"), str) and bool(item["reviewed_at_utc"].strip()),
                            f"water-orientation row {index} lacks a review timestamp",
                        )
                        qm_status = item.get("hf_631gd_distance_optimization_status")
                        expected_qm_by_disposition = {
                            "applicable": ["pass"],
                            "weak": ["pass"],
                            "unfavourable": ["pass"],
                            "excluded": ["not_required_prespecified_exclusion"],
                        }
                        audit.require(
                            item.get("allowed_hf_631gd_status_by_final_disposition")
                            == expected_qm_by_disposition,
                            f"water-orientation row {index} changed its frozen QM-status contract",
                        )
                        audit.require(
                            disposition in expected_qm_by_disposition
                            and qm_status in expected_qm_by_disposition[disposition],
                            f"water-orientation row {index} QM status is incompatible with its disposition",
                        )
                        selection_basis = item.get("selection_basis")
                        audit.require(
                            isinstance(selection_basis, str) and bool(selection_basis.strip()),
                            f"water-orientation row {index} lacks a selection basis",
                        )
                        if disposition == "excluded":
                            audit.require(
                                item.get("qm_input_artifact") is None
                                and item.get("qm_output_artifact") is None,
                                f"excluded water-orientation row {index} must not claim QM target artifacts",
                            )
                            validate_artifact(
                                package_root,
                                item.get("disposition_evidence_artifact"),
                                f"excluded water-orientation row {index} visual/chemical evidence",
                                audit,
                            )
                        elif disposition in {"applicable", "weak", "unfavourable"}:
                            validate_artifact(
                                package_root,
                                item.get("qm_input_artifact"),
                                f"water-orientation row {index} QM input",
                                audit,
                            )
                            validate_artifact(
                                package_root,
                                item.get("qm_output_artifact"),
                                f"water-orientation row {index} QM output",
                                audit,
                            )
                        failed_artifacts = item.get("failed_attempt_artifacts")
                        audit.require(
                            isinstance(failed_artifacts, list),
                            f"water-orientation row {index} failed-attempt artifacts must be a list",
                        )
                        if isinstance(failed_artifacts, list):
                            for failed_index, artifact in enumerate(failed_artifacts):
                                validate_artifact(
                                    package_root,
                                    artifact,
                                    f"water-orientation row {index} failed attempt {failed_index}",
                                    audit,
                                )
                    audit.require(observed_ids == expected_ids, "water-orientation IDs do not exactly cover O6U_WP_001 through O6U_WP_070")
                    audit.require(observed_lines == set(range(1, 71)), "water-orientation source lines do not exactly cover 1 through 70")
                    audit.require(observed_probe_counts == expected_probe_counts, "water-orientation probe-type counts differ from the frozen plan")
                    audit.require(
                        water_dispositions.get("final_disposition_counts") == observed_counts,
                        "water-orientation final disposition counts do not match the 70 rows",
                    )
        torsions = metrics.get("torsion_scans")
        if not isinstance(torsions, list) or not torsions:
            audit.errors.append("fit_metrics.torsion_scans must be a non-empty list")
        else:
            expected_torsion_targets = {
                "stiff_double_bond_N05_C27_C31_N07": "double_bond",
                "ring_internal_C32_N08_C31_N07": "five_membered_ring_internal",
                "coupled_ring_conjugated_N07_C31": "ring_internal_conjugated",
            }
            observed_target_ids: set[str] = set()
            for index, scan in enumerate(torsions):
                if not isinstance(scan, dict):
                    audit.errors.append(f"torsion scan {index} must be an object")
                    continue
                target_id = scan.get("target_id")
                audit.require(
                    isinstance(target_id, str) and target_id in expected_torsion_targets,
                    f"torsion scan {index} has an unknown or missing frozen target_id",
                )
                if isinstance(target_id, str):
                    audit.require(target_id not in observed_target_ids, f"torsion target {target_id} is duplicated")
                    observed_target_ids.add(target_id)
                    if target_id in expected_torsion_targets:
                        audit.require(
                            scan.get("central_bond_class") == expected_torsion_targets[target_id],
                            f"torsion target {target_id} has the wrong frozen chemical class",
                        )
                rmse = finite_number(scan.get("low_energy_rmse_kcal_mol"), f"torsion scan {index} RMSE", audit)
                minimum_error = finite_number(scan.get("minimum_position_error_degrees"), f"torsion scan {index} minimum error", audit)
                barrier_error = finite_number(scan.get("barrier_height_error_kcal_mol"), f"torsion scan {index} barrier error", audit)
                # These three values remain continuous annotations.  The O6U
                # protocol deliberately does not promote case-study values to
                # universal CGenFF cutoffs for stiff or ring-internal targets.
                if rmse is not None and minimum_error is not None and barrier_error is not None:
                    for key in (
                        "category_specific_domain_prespecified",
                        "energy_zero_alignment_verified",
                        "periodic_matching_verified",
                        "accessible_minima_reproduced",
                        "relevant_barriers_reproduced",
                        "conformer_ordering_reproduced",
                        "profile_shape_review_passed",
                        "other_targets_not_degraded",
                    ):
                        audit.require(scan.get(key) is True, f"torsion target {target_id or index} lacks {key}")
                    validate_artifact(
                        package_root,
                        scan.get("signed_review_artifact"),
                        f"torsion target {target_id or index} signed review",
                        audit,
                    )
            audit.require(
                observed_target_ids == set(expected_torsion_targets),
                "torsion validation does not exactly cover the three frozen O6U bonded targets",
            )
        regressions = metrics.get("charmm_to_gromacs_energy_component_regression")
        if not isinstance(regressions, list) or not regressions:
            audit.errors.append("CHARMM-to-GROMACS energy regression is missing")
        else:
            for index, item in enumerate(regressions):
                if not isinstance(item, dict):
                    audit.errors.append(f"energy regression {index} must be an object")
                    continue
                absolute = finite_number(item.get("absolute_difference_kcal_mol"), f"energy regression {index} absolute difference", audit)
                relative = finite_number(item.get("relative_difference_fraction"), f"energy regression {index} relative difference", audit)
                if absolute is not None and relative is not None:
                    audit.require(
                        abs(absolute) <= 0.1 or abs(relative) <= 0.001,
                        f"energy regression {index} exceeds both 0.1 kcal/mol and 0.1%",
                    )
    validate_artifact(package_root, record.get("atom_correspondence"), "O6U 76-atom correspondence table", audit)
    validate_artifact(package_root, record.get("final_topology"), "final O6U topology", audit)
    validate_artifact(package_root, record.get("final_parameters"), "final O6U parameters", audit)
    audit.require(isinstance(record.get("cgenff_version"), str) and bool(record["cgenff_version"].strip()) and not contains_todo(record["cgenff_version"]), "exact CGenFF version is unresolved")
    toolchain = record.get("toolchain_records")
    if not isinstance(toolchain, list):
        audit.errors.append("O6U toolchain_records must be a list")
    else:
        roles = [item.get("role") for item in toolchain if isinstance(item, dict)]
        audit.require(len(toolchain) == len(LIGAND_TOOLCHAIN_ROLES) and set(roles) == LIGAND_TOOLCHAIN_ROLES and len(roles) == len(set(roles)), "O6U toolchain roles are incomplete or duplicated")
        for index, item in enumerate(toolchain):
            if not isinstance(item, dict):
                audit.errors.append(f"O6U toolchain record {index} must be an object")
                continue
            role = str(item.get("role", f"row_{index}"))
            for key in ("name", "version"):
                value = item.get(key)
                audit.require(isinstance(value, str) and bool(value.strip()) and not contains_todo(value), f"O6U toolchain {role}.{key} is unresolved")
            digest = str(item.get("executable_or_container_sha256", "")).lower()
            audit.require(bool(SHA256_RE.fullmatch(digest)), f"O6U toolchain {role} executable/container digest is unresolved")
            validate_artifact(package_root, item.get("version_capture"), f"O6U toolchain {role} version capture", audit)
        cgenff_records = [item for item in toolchain if isinstance(item, dict) and item.get("role") == "cgenff_initial_assignment"]
        if len(cgenff_records) == 1:
            audit.require(cgenff_records[0].get("version") == record.get("cgenff_version"), "CGenFF version differs between summary and toolchain record")

    commands = record.get("command_records")
    if not isinstance(commands, list):
        audit.errors.append("O6U command_records must be a list")
    else:
        roles = [item.get("role") for item in commands if isinstance(item, dict)]
        audit.require(len(commands) == len(LIGAND_COMMAND_ROLES) and set(roles) == LIGAND_COMMAND_ROLES and len(roles) == len(set(roles)), "O6U command-record roles are incomplete or duplicated")
        for index, item in enumerate(commands):
            if not isinstance(item, dict):
                audit.errors.append(f"O6U command record {index} must be an object")
                continue
            validate_artifact(package_root, item.get("artifact"), f"O6U command record {item.get('role', index)}", audit)
    for key, label in (
        ("initial_penalty_inventory", "complete initial CGenFF penalty inventory"),
        ("qm_convergence_manifest", "QM convergence manifest"),
        ("parameter_change_table", "initial-to-final parameter change table"),
        ("raw_artifact_manifest", "complete ligand raw-artifact manifest"),
        ("cheap_canary_report", "cheap pre-membrane ligand canary report"),
    ):
        validate_artifact(package_root, record.get(key), label, audit)
    reviewers = record.get("reviewers")
    audit.require(isinstance(reviewers, list) and len(reviewers) >= 2 and all(isinstance(item, str) and item.strip() and not contains_todo(item) for item in reviewers), "O6U parameter record requires at least two named reviewers")
    audit.require(isinstance(record.get("approved_at_utc"), str) and bool(record["approved_at_utc"].strip()) and not contains_todo(record["approved_at_utc"]), "O6U parameter approval timestamp is unresolved")
    return record


def validate_construction_and_realizations(
    package_root: Path, manifest: dict[str, Any], stage: str, audit: Audit
) -> None:
    """Validate one immutable CHARMM-GUI construction and three velocity realizations."""
    system = manifest["systems"][0]
    construction = system.get("construction")
    if not isinstance(construction, dict):
        audit.errors.append(f"{PRIMARY_ID}.construction must be one object")
        return
    audit.require(construction.get("id") == "build01", "the sole construction ID must be build01")
    for key, label in (
        ("pdb_reader_jobid", "PDB Reader job ID"),
        ("quick_bilayer_jobid", "Quick Bilayer job ID"),
    ):
        value = construction.get(key)
        audit.require(isinstance(value, str) and bool(value.strip()) and not contains_todo(value), f"build01 lacks {label}")
    audit.require(construction.get("clone_job") is False, "build01 must be a directly audited construction, not a cloned job")
    audit.require(construction.get("ppm_applied_once") is True, "build01 must apply the PPM orientation exactly once")
    validate_orientation_record(
        package_root, construction.get("membrane_orientation_record"), construction,
        str(manifest.get("study_id", "")),
        audit,
    )
    archive = validate_artifact(package_root, construction.get("charmm_gui_archive"), "build01 CHARMM-GUI archive", audit)
    validate_artifact(package_root, construction.get("starting_coordinates"), "build01 starting coordinates", audit)
    for key, label in (
        ("topology", "topology"),
        ("index", "index"),
        ("analysis_index", "analysis index"),
        ("production_mdp", "500 ns production MDP"),
        ("minimization_mdp", "CHARMM-GUI minimization MDP"),
    ):
        validate_artifact(package_root, construction.get(key), f"build01 {label}", audit)
    common_min_value = construction.get("common_minimization_run_directory", "")
    if not isinstance(common_min_value, str) or not common_min_value.strip() or contains_todo(common_min_value):
        audit.errors.append("build01 common minimization run directory is unresolved")
    else:
        common_min = (package_root / common_min_value).resolve()
        audit.require(
            common_min != package_root and package_root in common_min.parents,
            "build01 common minimization run directory escapes the package",
        )
    directory_value = construction.get("gromacs_input_dir", "")
    if not isinstance(directory_value, str) or not directory_value.strip() or contains_todo(directory_value):
        audit.errors.append("build01 GROMACS input directory is unresolved")
    else:
        directory = (package_root / directory_value).resolve()
        audit.require(
            directory != package_root and package_root in directory.parents and directory.is_dir(),
            "build01 GROMACS input directory is missing or outside the package",
        )
        minimization_path = validate_artifact(
            package_root, construction.get("minimization_mdp"), "build01 CHARMM-GUI minimization MDP", audit
        )
        if minimization_path is not None:
            audit.require(
                minimization_path == (directory / "step6.0_minimization.mdp").resolve(),
                "build01 minimization MDP must be the exact extracted step6.0_minimization.mdp",
            )
        expected_equil = construction.get("equilibration_mdp_sha256")
        expected_names = [f"step6.{index}_equilibration.mdp" for index in range(1, 7)]
        audit.require(isinstance(expected_equil, dict) and sorted(expected_equil) == expected_names, "exact hashes for all six CHARMM-GUI equilibration MDPs are required")
        if isinstance(expected_equil, dict):
            for name in expected_names:
                expected_hash = str(expected_equil.get(name, ""))
                path = directory / name
                audit.require(bool(SHA256_RE.fullmatch(expected_hash)), f"{name} frozen SHA-256 is invalid")
                audit.require(path.is_file() and bool(SHA256_RE.fullmatch(expected_hash)) and sha256(path) == expected_hash, f"{name} differs from the exact audited CHARMM-GUI stage")
        tree_path = validate_artifact(package_root, construction.get("gromacs_input_tree_manifest"), "build01 complete GROMACS tree manifest", audit)
        if tree_path is not None and directory.is_dir():
            try:
                tree_errors = validate_tree_manifest_report(load_json(tree_path), directory, tree_path)
                audit.errors.extend(f"GROMACS tree manifest: {error}" for error in tree_errors)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                audit.errors.append(f"cannot validate GROMACS tree manifest: {exc}")
    report_path = validate_artifact(package_root, construction.get("build_validation_report"), "build01 validation report", audit)
    if report_path is not None:
        try:
            report = load_json(report_path)
            validate_build_report_binding(report, construction, manifest, audit)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            audit.errors.append(f"Cannot read build01 validation report: {exc}")
    production_artifact = construction.get("production_mdp")
    if isinstance(production_artifact, dict):
        production_path = validate_artifact(package_root, production_artifact, "build01 production MDP", audit)
        if production_path is not None:
            try:
                validate_production_mdp(production_path, manifest)
            except (ValueError, KeyError, TypeError) as exc:
                audit.errors.append(f"production MDP contract failed: {exc}")
    if stage == "production" and archive is None:
        audit.errors.append("production cannot start without the approved build01 archive")

    if stage in {"equilibration", "canary", "production"}:
        common_value = construction.get("common_minimization_run_directory", "")
        common_dir = (package_root / common_value).resolve() if isinstance(common_value, str) else package_root
        common_coordinates = common_dir / "work" / "step6.0_minimization.gro"
        common_record_path = common_dir / "common_minimization_record.json"
        if not common_coordinates.is_file() or not common_record_path.is_file():
            audit.errors.append("validated common-minimization output/record is required before independent equilibration")
        else:
            try:
                common_errors = validate_common_minimization_record(
                    package_root, manifest, construction, common_record_path
                )
                audit.errors.extend(common_errors)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                audit.errors.append(f"cannot validate common-minimization record: {exc}")
    if stage in {"builds", "equilibration", "canary", "production"}:
        validate_environment_report(package_root, construction.get("environment_validation_report"), manifest, audit)
    if stage in {"canary", "production"}:
        validate_md_stage_report(package_root, construction.get("equilibration_validation_report"), "equilibration", manifest, audit)
    if stage == "production":
        validate_md_stage_report(package_root, construction.get("canary_validation_report"), "canary", manifest, audit)

    realizations = system.get("realizations")
    audit.require(isinstance(realizations, list), f"{PRIMARY_ID}.realizations must be a list")
    if not isinstance(realizations, list):
        return
    ids = [str(item.get("id", "")) for item in realizations if isinstance(item, dict)]
    audit.require(ids == list(REALIZATION_IDS), f"realizations must be exactly {list(REALIZATION_IDS)} in order")
    audit.require(len(realizations) == 3, "exactly three velocity-seeded realizations are required")
    seeds: list[int] = []
    run_directories: list[str] = []
    forbidden_build_keys = {
        "pdb_reader_jobid", "quick_bilayer_jobid", "clone_job", "charmm_gui_archive",
        "gromacs_input_dir", "starting_coordinates", "topology", "index", "production_mdp", "minimization_mdp",
    }
    for item in realizations:
        if not isinstance(item, dict):
            audit.errors.append("a realization record is not an object")
            continue
        rid = str(item.get("id", "unknown"))
        seed = item.get("velocity_seed")
        audit.require(isinstance(seed, int) and seed > 0, f"{rid} velocity_seed must be a positive integer")
        if isinstance(seed, int):
            seeds.append(seed)
        duplicate_build_fields = sorted(forbidden_build_keys.intersection(item))
        audit.require(
            not duplicate_build_fields,
            f"{rid} duplicates construction fields {duplicate_build_fields}; all realizations must share build01",
        )
        run_value = item.get("run_directory", "")
        if not isinstance(run_value, str) or not run_value.strip() or contains_todo(run_value):
            audit.errors.append(f"{rid} run directory is unresolved")
        else:
            run = (package_root / run_value).resolve()
            audit.require(run != package_root and package_root in run.parents, f"{rid} run directory escapes the package")
            run_directories.append(str(run))
    audit.require(len(seeds) == 3 and len(seeds) == len(set(seeds)), "rep01-rep03 require three distinct velocity seeds")
    audit.require(len(run_directories) == 3 and len(run_directories) == len(set(run_directories)), "rep01-rep03 require distinct run directories")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument(
        "--stage",
        choices=("design", "chemistry", "builds", "equilibration", "canary", "production"),
        default="design",
    )
    parser.add_argument("--strict", action="store_true", help="Require the status and bound evidence for the selected stage; future-stage placeholders remain allowed")
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    package_root = manifest_path.parent.parent.resolve()
    audit = Audit()
    try:
        manifest = load_json(manifest_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: cannot read manifest: {exc}", file=sys.stderr)
        return 2

    audit.require(manifest.get("schema_version") == "3.0", "schema_version must be 3.0")
    validate_realization_template(package_root, audit)
    secret_paths = list(scan_secret_keys(manifest))
    audit.require(not secret_paths, "secret-like keys/values are prohibited in manifest: " + ", ".join(secret_paths))
    if args.strict:
        status = manifest.get("manifest_status")
        minimum = STAGE_MINIMUM_STATUS[args.stage]
        audit.require(status in RELEASE_STATUS_ORDER, f"strict mode requires a recognized staged manifest_status, found {status!r}")
        if status in RELEASE_STATUS_ORDER:
            audit.require(
                RELEASE_STATUS_ORDER.index(status) >= RELEASE_STATUS_ORDER.index(minimum),
                f"stage={args.stage} requires manifest_status at least {minimum}",
            )

    design = manifest.get("design", {})
    systems = manifest.get("systems", [])
    audit.require(design.get("type") == "single_native_holo_three_velocity_realizations", "unexpected design type")
    audit.require(design.get("biochemical_system_count") == 1, "exactly one biochemical system is required")
    audit.require(design.get("charmm_gui_construction_count") == 1, "exactly one CHARMM-GUI construction is required")
    audit.require(design.get("realization_ids") == list(REALIZATION_IDS), "design.realization_ids must be rep01-rep03")
    audit.require(
        design.get("realization_generation") == "one_common_deterministic_minimization_then_distinct_first_nvt_velocity_seeds_and_independent_dynamic_equilibration",
        "realization-generation contract is not frozen",
    )
    audit.require(design.get("sampling_unit") == "velocity_seeded_realization", "sampling unit must be velocity_seeded_realization")
    audit.require(design.get("comparative_inference_allowed") is False, "comparative inference must remain disabled")
    audit.require(design.get("allowed_primary_contrasts") == [], "no between-system contrast is allowed")
    audit.require(isinstance(systems, list), "systems must be a list")
    typed_systems = [item for item in systems if isinstance(item, dict)] if isinstance(systems, list) else []
    ids = {str(item.get("id")) for item in typed_systems}
    audit.require(ids == {PRIMARY_ID}, f"systems must contain only {PRIMARY_ID}")
    audit.require(len(typed_systems) == 1, "exactly one biochemical system is required")

    system = typed_systems[0] if len(typed_systems) == 1 else {}
    audit.require(system.get("role") == "native_inhibitor_bound_holo", "8KCT O6U system role is wrong")
    audit.require(str(system.get("pdb_id", "")).upper() == "8KCT", "native system PDB must be 8KCT")
    audit.require(str(system.get("ligand_component_id", "")).upper() == "O6U", "native system ligand must be O6U")
    audit.require(system.get("pose_provenance") == "experimentally_resolved_heavy_atom_pose", "native pose provenance is wrong")

    global_model = manifest.get("global_model", {})
    audit.require(manifest.get("study_id") == "gamma_secretase_native_nirogacestat_rebuild_20260808", "study_id is not the frozen study identifier")
    pressure = finite_number(global_model.get("pressure_bar"), "pressure_bar", audit)
    if pressure is not None:
        audit.require(abs(pressure - 1.0) < 1e-12, "pressure must be 1 bar")
    audit.require(global_model.get("membrane_upper") == global_model.get("membrane_lower"), "primary model requires symmetric leaflet composition")
    audit.require(global_model.get("membrane_upper") == "POPC", "bulk membrane must be pure POPC")
    audit.require(global_model.get("salt") == "NaCl", "solvent salt must be NaCl")
    salt_molar = finite_number(global_model.get("salt_molar"), "salt_molar", audit)
    if salt_molar is not None:
        audit.require(abs(salt_molar - 0.15) < 1e-12, "NaCl concentration must be 0.15 M")
    audit.require(global_model.get("resolved_native_lipids_retained") == {"CLR": 3, "PC1": 2}, "retain exactly the resolved native 3 CLR and 2 PC1 molecules")
    audit.require(global_model.get("nirogacestat_component_id") == "O6U", "the ligand component must be O6U")
    audit.require(
        global_model.get("nirogacestat_microstate") == "neutral PDB chemical-component microstate",
        "only the neutral PDB O6U microstate is allowed",
    )
    audit.require(global_model.get("nirogacestat_formal_charge_e") == 0, "O6U formal charge must be zero")
    asp257 = global_model.get("psen1_asp257_state")
    asp385 = global_model.get("psen1_asp385_state")
    if not contains_todo(asp257):
        audit.require(asp257 in ALLOWED_DYAD_STATES, "invalid Asp257 state")
    if not contains_todo(asp385):
        audit.require(asp385 in ALLOWED_DYAD_STATES, "invalid Asp385 state")
    audit.require(asp257 == "deprotonated", "Asp257 must match the frozen deprotonated model")
    audit.require(asp385 == "protonated_od2", "Asp385 must match the frozen OD2-protonated model")
    validate_artifact(package_root, global_model.get("dyad_rationale_record"), "dyad rationale", audit)

    simulation = manifest.get("simulation", {})
    validate_frozen_production_protocol_fields(manifest, audit)
    audit.require(simulation.get("analysis_window_ns") == [200.0, 500.0], "analysis window must be 200-500 ns")
    qc_gate = simulation.get("qc_and_stationarity_gate", {})
    audit.require(isinstance(qc_gate, dict), "simulation.qc_and_stationarity_gate must be an object")
    if isinstance(qc_gate, dict):
        audit.require(qc_gate.get("required_realization_ids") == list(REALIZATION_IDS), "QC gate must require rep01-rep03")
        audit.require(qc_gate.get("all_must_pass") is True, "all three realizations must pass QC and stationarity")
        audit.require(qc_gate.get("failure_outcome") == "analysis_inconclusive_fail_closed", "a failed realization must make analysis inconclusive")
        audit.require(qc_gate.get("trajectory_exclusion_allowed") is False, "failed realizations cannot be excluded")
        audit.require(qc_gate.get("analysis_cutoff_change_allowed") is False, "analysis cutoff cannot change after QC")
    audit.require(simulation.get("constraints") == "h-bonds", "constraints must be h-bonds")
    audit.require(simulation.get("electrostatics") == "PME", "electrostatics must be PME")
    audit.require(simulation.get("pressure_coupling") == "semiisotropic", "pressure coupling must be semiisotropic")
    audit.require(simulation.get("maxwarn") == 0, "grompp maxwarn must be 0")
    release = simulation.get("release_contract", {})
    audit.require(isinstance(release, dict), "simulation.release_contract must be an object")
    if isinstance(release, dict):
        audit.require(release.get("ordered_stages") == [
            "local_preflight", "environment_validation", "build_validation", "common_minimization",
            "three_independent_equilibrations", "three_same_tpr_5ns_canaries",
            "three_checkpoint_continuations_to_500ns",
        ], "release-stage order is not frozen")
        audit.require(release.get("canary_target_ns_per_realization") == 5.0 and release.get("canary_total_ns") == 15.0, "canary must be exactly 3 x 5 ns")
        audit.require(release.get("canary_uses_original_500ns_tpr") is True, "canary must use the original 500-ns TPR")
        audit.require(release.get("canary_frames_retained_as_production") is True, "canary frames must remain the start of production")
        audit.require(release.get("canary_selection_uses_ligand_behavior") is False, "ligand behavior cannot select canary realizations")
        audit.require(release.get("continuation_requires_all_three_canary_reports") is True, "continuation must require all three canary passes")
        audit.require(release.get("continuation_uses_same_tpr_and_checkpoint") is True, "continuation must reuse the exact TPR/checkpoint")
        audit.require(release.get("storage_headroom_fraction") == 0.30, "storage headroom must remain frozen at 30 percent")
        if args.stage in {"canary", "production"}:
            storage = finite_number(release.get("storage_budget_bytes"), "release storage_budget_bytes", audit)
            if storage is not None:
                audit.require(storage > 0 and float(storage).is_integer(), "release storage budget must be a positive whole-byte count")

    audit.require("docking_validation" not in manifest, "docking/self-redocking is outside this MD design")

    analysis_contract = manifest.get("analysis", {})
    audit.require(analysis_contract.get("role") == "descriptive_structural_stability_only", "analysis role must remain descriptive")
    audit.require(analysis_contract.get("between_system_statistics_allowed") is False, "between-system statistics must remain disabled")
    audit.require(analysis_contract.get("endpoint_energy_calculations_allowed") is False, "endpoint-energy calculations must remain disabled")
    audit.require(analysis_contract.get("pca_fel_role") == "hard_prohibited", "PCA/FEL must be hard-prohibited")
    audit.require(analysis_contract.get("pca_enabled") is False, "PCA execution must remain disabled")
    audit.require(analysis_contract.get("occupancy_derived_fel_enabled") is False, "occupancy-derived FEL must remain disabled")
    audit.require(analysis_contract.get("population_or_free_energy_surface_3d_enabled") is False, "3D population/free-energy surfaces must remain disabled")
    audit.require(analysis_contract.get("pca_fel_scripts_are_rejection_stubs") is True, "PCA/FEL executables must remain rejection stubs")
    audit.require(analysis_contract.get("pbc_distance_invariance_tolerance_nm") == 0.01, "PBC distance-invariance tolerance must be 0.01 nm")

    forbidden_serialized = json.dumps(manifest, sort_keys=True).lower()
    for forbidden in (
        "8kct_ligand_deleted",
        "matched_control_system_id",
        "mmpbsa",
        "mmgbsa",
        "self_redocking",
        "production_ns_initial",
        "recovery_extension",
        "extension_ns",
        "final_ns_if_approved",
    ):
        audit.require(forbidden not in forbidden_serialized, f"obsolete concept remains in manifest: {forbidden}")

    if args.stage in {"chemistry", "builds", "equilibration", "canary", "production"}:
        structure = validate_structure_record(package_root, str(system.get("structure_record", "")), "native system", audit)
        if structure:
            audit.require(structure.get("ligand_component_id") == "O6U", "native structure record lacks O6U")
            audit.require(structure.get("component_counts", {}).get("O6U") == 1, "native structure record must contain exactly one O6U")
        validate_parameter_record(package_root, str(system.get("ligand_parameter_record", "")), audit)

    if args.stage in {"builds", "equilibration", "canary", "production"}:
        validate_construction_and_realizations(package_root, manifest, args.stage, audit)

    if args.stage == "production":
        analysis = manifest.get("analysis", {})
        for key in (
            "analysis_plan_frozen_sha256",
            "pocket_fit_selection_record",
            "native_contact_definition_record",
            "bound_state_definition_record",
            "analysis_window_record",
        ):
            audit.require(not contains_todo(analysis.get(key)) and bool(analysis.get(key)), f"analysis.{key} is unresolved")
        audit.require(analysis.get("frames_are_independent_units") is False, "frames cannot be independent statistical units")
        digest = str(analysis.get("analysis_plan_frozen_sha256", ""))
        audit.require(bool(SHA256_RE.fullmatch(digest)), "analysis.analysis_plan_frozen_sha256 must be a lowercase SHA-256 digest")
        plan_record = analysis.get("analysis_plan")
        if not isinstance(plan_record, dict):
            audit.errors.append("analysis.analysis_plan must contain path and contract_sha256")
        else:
            plan_path = resolve_inside(package_root, plan_record.get("path", ""), "frozen analysis plan", audit)
            plan_digest = str(plan_record.get("contract_sha256", ""))
            audit.require(bool(SHA256_RE.fullmatch(plan_digest)), "analysis-plan contract SHA-256 is invalid")
            audit.require(plan_digest == digest, "analysis-plan record and analysis_plan_frozen_sha256 differ")
            if plan_path is not None:
                try:
                    plan = load_json(plan_path)
                    audit.require(analysis_plan_contract_sha256(plan) == plan_digest, "analysis-plan contract digest is stale")
                    prohibited = plan.get("prohibited_analyses")
                    audit.require(isinstance(prohibited, dict), "analysis plan must contain the hard-prohibited analysis policy")
                    if isinstance(prohibited, dict):
                        for key in (
                            "pca", "occupancy_derived_minus_kbt_ln_p", "free_energy_landscape",
                            "population_or_free_energy_surface_3d",
                        ):
                            audit.require(prohibited.get(key) is True, f"analysis plan must hard-prohibit {key}")
                        audit.require(
                            prohibited.get("policy") == "hard_prohibited_no_supplementary_exception",
                            "analysis plan allows a supplementary PCA/FEL exception",
                        )
                    for removed_key in ("pca", "grid", "support"):
                        audit.require(removed_key not in plan, f"analysis plan retains prohibited {removed_key} configuration")
                    outputs = plan.get("outputs", {})
                    audit.require(isinstance(outputs, dict), "analysis plan outputs must be an object")
                    if isinstance(outputs, dict):
                        for key in ("pca_outputs", "fel_outputs", "population_or_free_energy_surface_3d_outputs"):
                            audit.require(outputs.get(key) is False, f"analysis plan output {key} must remain disabled")
                except (OSError, ValueError, json.JSONDecodeError) as exc:
                    audit.errors.append(f"cannot validate analysis-plan contract: {exc}")
        software = manifest.get("software", {})
        audit.require(not contains_todo(software) and bool(software), "software versions/container are unresolved")

    for warning in audit.warnings:
        print(f"WARNING: {warning}")
    for error in audit.errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if audit.errors:
        print(f"Preflight failed with {len(audit.errors)} error(s).", file=sys.stderr)
        return 1
    print(f"Preflight passed for stage={args.stage}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
