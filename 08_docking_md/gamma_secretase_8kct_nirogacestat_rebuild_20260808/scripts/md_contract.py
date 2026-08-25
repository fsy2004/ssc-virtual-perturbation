#!/usr/bin/env python3
"""Shared immutable-contract helpers for the 8KCT MD release pipeline."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


SYSTEM_ID = "8kct_nirogacestat_native"
BUILD_ID = "build01"
REALIZATION_IDS = ("rep01", "rep02", "rep03")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def artifact_path(package_root: Path, artifact: Any, label: str, *, must_exist: bool = True) -> Path:
    if not isinstance(artifact, dict):
        raise ValueError(f"{label} must be an object with path and sha256")
    relative = artifact.get("path")
    expected = str(artifact.get("sha256", "")).lower()
    if not isinstance(relative, str) or not relative.strip():
        raise ValueError(f"{label}.path is unresolved")
    path = (package_root / relative).resolve()
    if path != package_root and package_root not in path.parents:
        raise ValueError(f"{label}.path escapes package root: {relative}")
    if not SHA256_RE.fullmatch(expected):
        raise ValueError(f"{label}.sha256 is not a lowercase SHA-256 digest")
    if must_exist:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"{label} is missing or empty: {path}")
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"{label} SHA-256 mismatch: expected {expected}, observed {observed}")
    return path


def validate_bound_gmx_identity(package_root: Path, manifest: dict[str, Any], gmx: str) -> dict[str, Any]:
    """Bind the invoked GROMACS binary and linked core library to the sealed environment report."""
    systems = manifest.get("systems", [])
    if not isinstance(systems, list) or len(systems) != 1:
        raise ValueError("manifest must contain exactly one system before GROMACS identity validation")
    construction = systems[0].get("construction", {})
    environment_path = artifact_path(
        package_root, construction.get("environment_validation_report"), "server environment validation report"
    )
    environment = load_json(environment_path)
    if environment.get("schema_version") != "1.0" or environment.get("report_type") != "server_environment_validation":
        raise ValueError("server environment validation report schema/type is invalid")
    if environment.get("status") != "pass" or environment.get("approval_status") != "approved":
        raise ValueError("server environment is not approved PASS")
    integrity = environment.get("integrity", {})
    if not isinstance(integrity, dict) or integrity.get("payload_sha256") != report_payload_sha256(
        environment, ("integrity", "payload_sha256")
    ):
        raise ValueError("server environment report payload checksum is invalid")
    identity_path = artifact_path(package_root, environment.get("gmx_executable"), "GROMACS executable identity")
    identity = load_json(identity_path)
    if identity.get("schema_version") != "1.0" or identity.get("record_type") != "gromacs_executable_identity":
        raise ValueError("GROMACS executable identity schema/type is invalid")
    resolved_text = shutil.which(gmx)
    if resolved_text is None and Path(gmx).is_file():
        resolved_text = str(Path(gmx))
    if resolved_text is None:
        raise FileNotFoundError(f"GROMACS executable not found: {gmx}")
    resolved = Path(resolved_text).resolve()
    if str(resolved) != identity.get("resolved_path"):
        raise ValueError(f"Invoked GROMACS path differs from the approved binary: {resolved}")
    if identity.get("bytes") != resolved.stat().st_size or identity.get("sha256") != sha256(resolved):
        raise ValueError("Invoked GROMACS binary differs from its approved hash/size")
    version_result = subprocess.run([str(resolved), "--version"], text=True, capture_output=True, check=False)
    if version_result.returncode != 0:
        raise RuntimeError("Approved GROMACS binary failed --version")
    version_match = re.search(r"^GROMACS version:\s*(\S+)", version_result.stdout + "\n" + version_result.stderr, re.MULTILINE)
    observed_version = version_match.group(1) if version_match else None
    expected_version = manifest.get("simulation", {}).get("required_version")
    if observed_version != expected_version or observed_version != environment.get("gromacs_version") or observed_version != identity.get("gromacs_version"):
        raise ValueError("Invoked GROMACS version differs from the manifest/environment/identity contract")
    linked = identity.get("linked_libraries")
    if not isinstance(linked, list) or not linked:
        raise ValueError("GROMACS executable identity lacks linked-library hashes")
    for record in linked:
        if not isinstance(record, dict):
            raise ValueError("GROMACS linked-library record is invalid")
        library = Path(str(record.get("path", ""))).resolve()
        if not library.is_file() or library.stat().st_size != record.get("bytes") or sha256(library) != record.get("sha256"):
            raise ValueError(f"Approved GROMACS linked library differs: {library}")
    return {
        "resolved_path": str(resolved),
        "sha256": sha256(resolved),
        "bytes": resolved.stat().st_size,
        "gromacs_version": observed_version,
        "identity_record": str(identity_path),
        "identity_record_sha256": sha256(identity_path),
    }


def parse_mdp(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        normalized = key.strip().lower().replace("-", "_")
        if normalized in values:
            raise ValueError(f"{path}: duplicate active MDP key {key.strip()!r}")
        values[normalized] = value.strip()
    return values


def _norm(value: str) -> str:
    return re.sub(r"[-_\s]", "", value.strip().lower())


def _bool_yes(value: str) -> bool:
    return _norm(value) in {"yes", "true"}


def _numbers(value: str, label: str) -> list[float]:
    try:
        result = [float(item) for item in value.split()]
    except ValueError as exc:
        raise ValueError(f"{label} must contain only numeric tokens") from exc
    if not result or not all(math.isfinite(item) for item in result):
        raise ValueError(f"{label} contains no finite numeric values")
    return result


def _require_numeric_equal(values: list[float], expected: list[float], label: str, tolerance: float = 1e-12) -> None:
    if len(values) != len(expected) or any(abs(left - right) > tolerance for left, right in zip(values, expected)):
        raise ValueError(f"{label} must be {expected}, found {values}")


def validate_production_mdp(path: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    """Validate the entire frozen production physics/output contract."""
    simulation = manifest.get("simulation", {})
    contract = simulation.get("production_mdp_contract", {})
    if not isinstance(contract, dict):
        raise ValueError("simulation.production_mdp_contract must be an object")
    mdp = parse_mdp(path)
    required = {
        "integrator", "dt", "nsteps", "continuation", "gen_vel", "cutoff_scheme", "nstlist",
        "rlist", "coulombtype", "rcoulomb", "pme_order", "fourierspacing", "vdwtype",
        "vdw_modifier", "rvdw_switch", "rvdw", "dispcorr", "constraints", "constraint_algorithm",
        "tcoupl", "tc_grps", "tau_t", "ref_t", "pcoupl", "pcoupltype", "tau_p", "ref_p",
        "compressibility", "comm_mode", "comm_grps", "nstcomm", "nstxout", "nstvout", "nstfout",
        "nstxout_compressed", "nstcalcenergy", "nstenergy", "nstlog", "compressed_x_precision",
        "pbc", "periodic_molecules",
    }
    missing = sorted(required.difference(mdp))
    if missing:
        raise ValueError(f"{path}: missing required production MDP fields: {missing}")

    expected_ns = float(simulation.get("production_ns"))
    expected_timestep = float(simulation.get("time_step_ps"))
    timestep = float(mdp["dt"])
    nsteps = int(mdp["nsteps"])
    if abs(timestep - expected_timestep) > 1e-12:
        raise ValueError(f"{path}: dt must match the frozen value {expected_timestep} ps")
    expected_steps = int(round(expected_ns * 1000.0 / expected_timestep))
    if nsteps != expected_steps or abs(timestep * nsteps / 1000.0 - expected_ns) > 1e-9:
        raise ValueError(
            f"{path}: nsteps must encode exactly {expected_ns:g} ns at {expected_timestep:g} ps"
        )
    if _norm(mdp["integrator"]) != _norm(str(contract.get("integrator"))):
        raise ValueError(f"{path}: integrator must be md")
    if not _bool_yes(mdp["continuation"]) or _bool_yes(mdp["gen_vel"]):
        raise ValueError(f"{path}: production requires continuation=yes and gen_vel=no")
    if _norm(mdp["pbc"]) != "xyz" or _bool_yes(mdp["periodic_molecules"]):
        raise ValueError(f"{path}: production requires pbc=xyz and periodic-molecules=no")
    if simulation.get("hydrogen_mass_repartitioning") is not True:
        raise ValueError("simulation.hydrogen_mass_repartitioning must be true for the frozen protocol")
    # HMR is encoded in the frozen topology; the production MDP need not and in
    # this release does not contain a mass-repartition-factor directive.
    if "mass_repartition_factor" in mdp and float(mdp["mass_repartition_factor"]) <= 1.0:
        raise ValueError(f"{path}: mass-repartition-factor contradicts the frozen HMR topology")

    prohibited_exact = {
        "define", "pull", "awh", "free_energy", "simulated_tempering", "annealing", "mts",
        "deform", "accelerate", "cos_acceleration", "freezegrps", "freeze_dim",
        "energygrp_excl", "rotation", "swapcoords", "adress", "qmmm", "nwall",
    }
    prohibited_prefixes = (
        "pull_", "awh_", "fep_", "coul_lambdas", "vdw_lambdas", "bonded_lambdas",
        "restraint_lambdas", "temperature_lambdas", "calc_lambda", "init_lambda", "delta_lambda",
        "sc_", "sim_temp", "annealing_", "mts_", "electric_field_", "disre", "orire", "dihre",
        "wall_",
    )
    prohibited_present = sorted(
        key for key in mdp if key in prohibited_exact or any(key.startswith(prefix) for prefix in prohibited_prefixes)
    )
    if prohibited_present:
        raise ValueError(f"{path}: production MDP contains prohibited restraint/bias/HMR-adjacent controls: {prohibited_present}")

    exact_text = {
        "cutoff_scheme": contract.get("cutoff_scheme"),
        "coulombtype": "PME",
        "vdwtype": contract.get("vdw_type"),
        "vdw_modifier": contract.get("vdw_modifier"),
        "dispcorr": contract.get("dispersion_correction"),
        "constraints": simulation.get("constraints"),
        "constraint_algorithm": contract.get("constraint_algorithm"),
        "tcoupl": contract.get("thermostat"),
        "pcoupl": contract.get("barostat"),
        "pcoupltype": simulation.get("pressure_coupling"),
        "comm_mode": contract.get("com_removal_mode"),
    }
    for key, expected in exact_text.items():
        if not isinstance(expected, str) or _norm(mdp[key]) != _norm(expected):
            raise ValueError(f"{path}: {key} must match frozen value {expected!r}, found {mdp[key]!r}")

    exact_scalar = {
        "nstlist": contract.get("neighbor_list_update_steps"),
        "rlist": contract.get("rlist_nm"),
        "rcoulomb": contract.get("rcoulomb_nm"),
        "pme_order": contract.get("pme_order"),
        "fourierspacing": contract.get("fourier_spacing_nm"),
        "rvdw_switch": contract.get("rvdw_switch_nm"),
        "rvdw": contract.get("rvdw_nm"),
        "tau_p": contract.get("barostat_tau_p_ps"),
        "nstcomm": contract.get("com_removal_interval_steps"),
    }
    for key, expected_value in exact_scalar.items():
        expected = float(expected_value)
        observed = float(mdp[key])
        if not math.isfinite(observed) or abs(observed - expected) > 1e-12:
            raise ValueError(f"{path}: {key} must be {expected_value}, found {mdp[key]}")

    thermostat_groups = contract.get("thermostat_groups")
    if not isinstance(thermostat_groups, list) or not thermostat_groups or any(
        not isinstance(item, str) or not item.strip() or "TODO" in item.upper() for item in thermostat_groups
    ):
        raise ValueError("production thermostat_groups are unresolved")
    if mdp["tc_grps"].split() != thermostat_groups:
        raise ValueError(f"{path}: tc-grps must be {thermostat_groups}, found {mdp['tc_grps'].split()}")
    tau_t = _numbers(mdp["tau_t"], "tau-t")
    ref_t = _numbers(mdp["ref_t"], "ref-t")
    _require_numeric_equal(tau_t, [float(contract.get("tau_t_ps"))] * len(thermostat_groups), "tau-t")
    _require_numeric_equal(ref_t, [float(manifest["global_model"]["temperature_k"])] * len(thermostat_groups), "ref-t")
    _require_numeric_equal(_numbers(mdp["ref_p"], "ref-p"), [1.0, 1.0], "ref-p")
    _require_numeric_equal(
        _numbers(mdp["compressibility"], "compressibility"),
        [float(value) for value in contract.get("compressibility_bar_inverse", [])],
        "compressibility",
    )
    com_groups = contract.get("com_removal_groups")
    if not isinstance(com_groups, list) or mdp["comm_grps"].split() != com_groups:
        raise ValueError(f"{path}: comm-grps must be the frozen list {com_groups}")

    cadence = contract.get("output_cadence_steps")
    if not isinstance(cadence, dict):
        raise ValueError("production output_cadence_steps must be an object")
    expected_output = {
        "nstxout": "nstxout",
        "nstvout": "nstvout",
        "nstfout": "nstfout",
        "nstxout_compressed": "nstxout_compressed",
        "nstcalcenergy": "nstcalcenergy",
        "nstenergy": "nstenergy",
        "nstlog": "nstlog",
        "compressed_x_precision": "compressed_x_precision",
    }
    for mdp_key, contract_key in expected_output.items():
        expected = int(cadence.get(contract_key))
        observed = int(mdp[mdp_key])
        if observed != expected:
            raise ValueError(f"{path}: {mdp_key} must be {expected}, found {observed}")
    if int(mdp["nstxout_compressed"]) <= 0 or int(mdp["nstenergy"]) <= 0 or int(mdp["nstlog"]) <= 0:
        raise ValueError(f"{path}: compressed trajectory, energy, and log cadence must be positive")

    return {
        "path": str(path),
        "sha256": sha256(path),
        "duration_ns": timestep * nsteps / 1000.0,
        "temperature_k": float(manifest["global_model"]["temperature_k"]),
        "pressure_bar": float(manifest["global_model"].get("pressure_bar", 1.0)),
        "pbc": "xyz",
        "hydrogen_mass_repartitioning": True,
        "output_cadence_steps": {key: int(mdp[key]) for key in expected_output},
    }


def build_contract_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    systems = manifest.get("systems", [])
    if not isinstance(systems, list) or len(systems) != 1:
        raise ValueError("manifest must contain one system")
    system = systems[0]
    construction = system.get("construction", {})
    keys = (
        "id", "pdb_reader_jobid", "quick_bilayer_jobid", "clone_job", "ppm_applied_once",
        "membrane_orientation_record", "charmm_gui_archive", "gromacs_input_dir",
        "gromacs_input_tree_manifest", "starting_coordinates", "topology", "index", "analysis_index",
        "production_mdp", "minimization_mdp", "equilibration_mdp_sha256",
    )
    return {
        "schema_version": "build_contract_v1",
        "study_id": manifest.get("study_id"),
        "system_id": system.get("id"),
        "pdb_id": system.get("pdb_id"),
        "ligand_component_id": system.get("ligand_component_id"),
        "construction": {key: copy.deepcopy(construction.get(key)) for key in keys},
    }


def build_contract_sha256(manifest: dict[str, Any]) -> str:
    return canonical_json_sha256(build_contract_payload(manifest))


def release_manifest_contract_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(manifest)
    if "manifest_status" in payload:
        payload["manifest_status"] = "STAGED_STATUS_BOUND_SEPARATELY"
    systems = payload.get("systems", [])
    if isinstance(systems, list) and len(systems) == 1 and isinstance(systems[0], dict):
        construction = systems[0].get("construction", {})
        if isinstance(construction, dict):
            for key in (
                "build_validation_report", "environment_validation_report",
                "equilibration_validation_report", "canary_validation_report",
            ):
                artifact = construction.get(key)
                if isinstance(artifact, dict):
                    artifact["sha256"] = "BOUND_AFTER_REPORT_CREATION"
    return payload


def release_manifest_contract_sha256(manifest: dict[str, Any]) -> str:
    return canonical_json_sha256(release_manifest_contract_payload(manifest))


def analysis_plan_contract_payload(plan: dict[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(plan)
    eligibility = payload.setdefault("eligibility_gate", {})
    eligibility["all_realizations_passed_qc_and_stationarity"] = "BOUND_AT_ADJUDICATION"
    eligibility["qc_and_stationarity_report_sha256"] = "BOUND_AT_ADJUDICATION"
    return payload


def analysis_plan_contract_sha256(plan: dict[str, Any]) -> str:
    return canonical_json_sha256(analysis_plan_contract_payload(plan))


def report_payload_sha256(report: dict[str, Any], field_path: tuple[str, ...]) -> str:
    payload = copy.deepcopy(report)
    cursor: Any = payload
    for key in field_path[:-1]:
        cursor = cursor[key]
    cursor[field_path[-1]] = "UNSEALED"
    return canonical_json_sha256(payload)


def validate_common_minimization_record(
    package_root: Path,
    manifest: dict[str, Any],
    construction: dict[str, Any],
    record_path: Path,
) -> list[str]:
    """Validate the immutable technical release record for the shared minimization."""
    errors: list[str] = []
    try:
        record = load_json(record_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return [f"cannot read common-minimization record: {exc}"]
    expected_identity = {
        "schema_version": "2.0",
        "report_type": "common_minimization_validation",
        "status": "pass",
        "technical_integrity_pass": True,
        "system_id": SYSTEM_ID,
        "construction_id": BUILD_ID,
        "build_contract_sha256": build_contract_sha256(manifest),
        "archive_sha256": construction.get("charmm_gui_archive", {}).get("sha256"),
        "gromacs_input_tree_manifest_sha256": construction.get("gromacs_input_tree_manifest", {}).get("sha256"),
        "minimization_mdp_sha256": construction.get("minimization_mdp", {}).get("sha256"),
        "grompp_warning_count": 0,
        "forbidden_runtime_pattern_count": 0,
    }
    for key, expected in expected_identity.items():
        if record.get(key) != expected:
            errors.append(f"common-minimization record mismatch at {key}")
    artifacts = record.get("artifacts")
    required_artifacts = ("mdp", "tpr", "energy", "log", "coordinates")
    if not isinstance(artifacts, dict) or sorted(artifacts) != sorted(required_artifacts):
        errors.append(f"common-minimization artifacts must be exactly {list(required_artifacts)}")
    else:
        resolved: dict[str, Path] = {}
        for label in required_artifacts:
            artifact = artifacts.get(label)
            try:
                path = artifact_path(package_root, artifact, f"common minimization {label}")
                resolved[label] = path
                if artifact.get("bytes") != path.stat().st_size:
                    errors.append(f"common-minimization {label} byte count is stale")
            except (OSError, ValueError) as exc:
                errors.append(str(exc))
        if "mdp" in resolved and sha256(resolved["mdp"]) != construction.get("minimization_mdp", {}).get("sha256"):
            errors.append("common-minimization MDP differs from the frozen build input")
        if "coordinates" in resolved:
            if record.get("coordinates_sha256") != sha256(resolved["coordinates"]):
                errors.append("common-minimization coordinate convenience hash is stale")
    commands = record.get("command_records")
    if not isinstance(commands, list) or len(commands) < 2:
        errors.append("common-minimization record lacks retained grompp/mdrun command provenance")
    else:
        previous_hash: str | None = None
        command_kinds: list[str] = []
        for index, artifact in enumerate(commands):
            try:
                finished_path = artifact_path(package_root, artifact, f"common minimization command record {index}")
                finished = load_json(finished_path)
                started_path = Path(str(finished.get("started_record", ""))).resolve()
                if started_path != package_root and package_root not in started_path.parents:
                    raise ValueError("common-minimization started record escapes package root")
                if not started_path.is_file() or finished.get("started_record_sha256") != sha256(started_path):
                    raise ValueError("common-minimization started record is missing or changed")
                started = load_json(started_path)
                if started.get("previous_finished_record_sha256") != previous_hash:
                    errors.append(f"common-minimization command hash chain breaks at record {index}")
                argv = started.get("argv")
                if not isinstance(argv, list) or len(argv) < 2 or not all(isinstance(item, str) for item in argv):
                    errors.append(f"common-minimization command {index} has invalid argv")
                else:
                    command_kinds.append(argv[1])
                    if "-maxwarn" in argv:
                        errors.append(f"common-minimization command {index} uses prohibited -maxwarn")
                if finished.get("returncode") != 0:
                    errors.append(f"common-minimization command {index} did not exit successfully")
                for stream in ("stdout", "stderr"):
                    stream_artifact = finished.get(stream)
                    if not isinstance(stream_artifact, dict):
                        errors.append(f"common-minimization command {index} lacks {stream} provenance")
                    else:
                        stream_path = Path(str(stream_artifact.get("path", ""))).resolve()
                        if stream_path != package_root and package_root not in stream_path.parents:
                            errors.append(f"common-minimization command {index} {stream} escapes package root")
                        elif not stream_path.is_file() or stream_artifact.get("sha256") != sha256(stream_path):
                            errors.append(f"common-minimization command {index} {stream} is missing or changed")
                previous_hash = sha256(finished_path)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                errors.append(str(exc))
        if command_kinds.count("grompp") != 1 or command_kinds.count("mdrun") < 1:
            errors.append("common minimization requires one grompp and at least one mdrun segment")
    integrity = record.get("integrity")
    if not isinstance(integrity, dict) or integrity.get("payload_sha256") != report_payload_sha256(
        record, ("integrity", "payload_sha256")
    ):
        errors.append("common-minimization record payload checksum is invalid")
    return errors
