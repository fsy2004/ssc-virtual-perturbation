#!/usr/bin/env python3
"""Restart-safe GROMACS runner for one build and three velocity realizations.

Planning is the default. Execution requires the strict preflight for the requested stage.
One deterministic minimization is run once, then copied immutably to all three
realizations; distinct velocities are generated at the first NVT stage.
All three fixed 500 ns realizations are required; this runner provides no
protocol path for changing the cutoff or omitting a realization.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from md_contract import (
    artifact_path,
    build_contract_sha256,
    parse_mdp,
    report_payload_sha256,
    sha256,
    validate_bound_gmx_identity,
    validate_common_minimization_record,
    validate_production_mdp as validate_frozen_production_mdp,
)


SYSTEM_ID = "8kct_nirogacestat_native"
REALIZATION_IDS = ("rep01", "rep02", "rep03")
EQUILIBRATION_STAGES = tuple(f"step6.{index}_equilibration" for index in range(1, 7))
GROMPP_WARNING = re.compile(r"^\s*WARNING\s+[0-9]+", re.IGNORECASE | re.MULTILINE)
FORBIDDEN_RUNTIME_PATTERNS = (
    re.compile(r"fatal error", re.IGNORECASE),
    re.compile(r"segmentation fault", re.IGNORECASE),
    re.compile(r"\b(?:nan|inf)\b", re.IGNORECASE),
    re.compile(r"lincs warning", re.IGNORECASE),
    re.compile(r"settle.*(?:error|warning)", re.IGNORECASE),
    re.compile(r"constraint.*(?:error|warning)", re.IGNORECASE),
    re.compile(r"pressure scaling more than", re.IGNORECASE),
)
PERFORMANCE_MDRUN_FLAGS = {
    "-nt": 1,
    "-ntmpi": 1,
    "-ntomp": 1,
    "-ntomp_pme": 1,
    "-pin": 1,
    "-pinoffset": 1,
    "-pinstride": 1,
    "-gpu_id": 1,
    "-gputasks": 1,
    "-nb": 1,
    "-pme": 1,
    "-pmefft": 1,
    "-bonded": 1,
    "-update": 1,
    "-npme": 1,
    "-dlb": 1,
    "-dds": 1,
    "-rdd": 1,
    "-rcon": 1,
    "-maxh": 1,
    "-v": 0,
    "-tunepme": 0,
    "-notunepme": 0,
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def inside(package_root: Path, value: str, *, must_exist: bool = True) -> Path:
    path = (package_root / value).resolve()
    if path != package_root and package_root not in path.parents:
        raise ValueError(f"Path escapes package root: {value}")
    if must_exist and not path.exists():
        raise FileNotFoundError(path)
    return path


def validate_mdrun_extra(tokens: list[str]) -> None:
    """Permit only resource/performance flags; never allow run identity or physics overrides."""
    index = 0
    while index < len(tokens):
        flag = tokens[index]
        arity = PERFORMANCE_MDRUN_FLAGS.get(flag)
        if arity is None:
            allowed = ", ".join(sorted(PERFORMANCE_MDRUN_FLAGS))
            raise ValueError(f"Unsupported --mdrun-extra token {flag!r}; allowed flags: {allowed}")
        if index + arity >= len(tokens):
            raise ValueError(f"{flag} requires {arity} value token(s)")
        for offset in range(1, arity + 1):
            value = tokens[index + offset]
            if value.startswith("-"):
                raise ValueError(f"{flag} has missing/invalid value before {value!r}")
        index += arity + 1


def patch_velocity_seed(source: Path, destination: Path, seed: int) -> dict[str, Any]:
    if seed <= 0:
        raise ValueError("velocity_seed must be a positive integer")
    settings = parse_mdp(source)
    if settings.get("integrator", "").strip().lower() not in {"md", "md-vv", "md_vv"}:
        raise ValueError(f"{source}: first seeded stage must be molecular dynamics")
    if settings.get("tcoupl", "no").strip().lower() == "no":
        raise ValueError(f"{source}: first seeded stage must be NVT with temperature coupling")
    if settings.get("pcoupl", "no").strip().lower() != "no":
        raise ValueError(f"{source}: velocity generation must occur in the first NVT stage before pressure coupling")
    if settings.get("continuation", "no").strip().lower() in {"yes", "true"}:
        raise ValueError(f"{source}: first NVT generates new velocities and must not set continuation=yes")
    lines = source.read_text(encoding="utf-8", errors="strict").splitlines()
    found_seed = False
    generation_enabled = False
    output: list[str] = []
    for line in lines:
        content = line.split(";", 1)[0]
        if "=" in content:
            key, value = (part.strip() for part in content.split("=", 1))
            normalized = key.lower().replace("-", "_")
            if normalized == "gen_seed":
                output.append(f"gen_seed = {seed} ; frozen by run_md_matrix.py from {source.name}")
                found_seed = True
                continue
            if normalized == "gen_vel":
                generation_enabled = value.lower() in {"yes", "true"}
        output.append(line)
    if not found_seed or not generation_enabled:
        raise ValueError(f"{source}: first equilibration MDP must contain gen_vel=yes and gen_seed")
    rendered = "\n".join(output) + "\n"
    if destination.exists() and destination.read_text(encoding="utf-8") != rendered:
        raise RuntimeError(f"Refusing to replace a changed frozen seed MDP: {destination}")
    if not destination.exists():
        destination.write_text(rendered, encoding="utf-8", newline="\n")
    return {
        "source": str(source),
        "source_sha256": sha256(source),
        "derived": str(destination),
        "derived_sha256": sha256(destination),
        "velocity_seed": seed,
    }


def common_minimization_plan(gmx: str) -> list[list[str]]:
    return [
        [gmx, "grompp", "-f", "step6.0_minimization.mdp", "-o", "step6.0_minimization.tpr", "-c", "step5_input.gro", "-r", "step5_input.gro", "-p", "topol.top", "-n", "index.ndx"],
        [gmx, "mdrun", "-deffnm", "step6.0_minimization"],
    ]


def realization_command_plan(gmx: str, phase: str, canary_steps: int = 2_500_000) -> list[list[str]]:
    commands: list[list[str]] = []
    if phase == "equilibration":
        previous = "step6.0_minimization"
        for stage in EQUILIBRATION_STAGES:
            mdp = "step6.1_equilibration.frozen.mdp" if stage == EQUILIBRATION_STAGES[0] else f"{stage}.mdp"
            grompp = [gmx, "grompp", "-f", mdp, "-o", f"{stage}.tpr", "-c", f"{previous}.gro"]
            if stage != EQUILIBRATION_STAGES[0]:
                grompp.extend(["-t", f"{previous}.cpt"])
            grompp.extend(["-r", "step5_input.gro", "-p", "topol.top", "-n", "index.ndx"])
            commands.extend([grompp, [gmx, "mdrun", "-deffnm", stage]])
            previous = stage
    if phase == "canary":
        commands.extend([
            [gmx, "grompp", "-f", "production.frozen.mdp", "-o", "production.tpr", "-c", "step6.6_equilibration.gro", "-t", "step6.6_equilibration.cpt", "-p", "topol.top", "-n", "index.ndx"],
            [gmx, "mdrun", "-s", "production.tpr", "-deffnm", "production", "-nsteps", str(canary_steps)],
        ])
    elif phase == "production":
        commands.append([
            gmx, "mdrun", "-s", "production.tpr", "-deffnm", "production",
            "-cpi", "production.cpt", "-append",
        ])
    return commands


def validate_staged_mdp_hashes(source: Path, construction: dict[str, Any]) -> None:
    minimization = construction.get("minimization_mdp")
    expected_minimization_hash = minimization.get("sha256") if isinstance(minimization, dict) else None
    minimization_path = source / "step6.0_minimization.mdp"
    if not isinstance(expected_minimization_hash, str) or sha256(minimization_path) != expected_minimization_hash:
        raise ValueError("CHARMM-GUI minimization MDP hash differs from the frozen build contract")
    expected = construction.get("equilibration_mdp_sha256")
    if not isinstance(expected, dict) or sorted(expected) != [f"step6.{index}_equilibration.mdp" for index in range(1, 7)]:
        raise ValueError("build01 must freeze exact SHA-256 values for all six CHARMM-GUI equilibration MDPs")
    for name, expected_hash in expected.items():
        path = source / name
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"CHARMM-GUI staged equilibration MDP hash mismatch: {path}")


def stage_common_minimization(package_root: Path, construction: dict[str, Any]) -> tuple[Path, Path]:
    source = inside(package_root, str(construction["gromacs_input_dir"]))
    run_dir = inside(package_root, str(construction["common_minimization_run_directory"]), must_exist=False)
    work = run_dir / "work"
    source_record = run_dir / "staged_source.json"
    identity = {
        "construction_id": construction["id"],
        "source": str(source),
        "archive_sha256": construction["charmm_gui_archive"]["sha256"],
        "gromacs_input_tree_manifest_sha256": construction["gromacs_input_tree_manifest"]["sha256"],
        "stage": "one_common_deterministic_minimization",
    }
    if work.exists():
        if not source_record.is_file():
            raise RuntimeError(f"Existing common-minimization work directory has no provenance: {work}")
        recorded = load_json(source_record)
        for key, value in identity.items():
            if recorded.get(key) != value:
                raise RuntimeError(f"Common-minimization provenance differs: {work}")
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, work)
        source_record.write_text(json.dumps({**identity, "staged_at_utc": utc_now()}, indent=2) + "\n", encoding="utf-8")
    validate_staged_mdp_hashes(source, construction)
    minimization = parse_mdp(work / "step6.0_minimization.mdp")
    if minimization.get("integrator", "").strip().lower() not in {"steep", "cg", "l-bfgs"}:
        raise ValueError("step6.0_minimization.mdp is not a deterministic energy-minimization stage")
    if minimization.get("gen_vel", "no").strip().lower() in {"yes", "true"}:
        raise ValueError("common minimization must not generate velocities")
    return run_dir, work


def stage_inputs(
    package_root: Path,
    manifest: dict[str, Any],
    construction: dict[str, Any],
    realization: dict[str, Any],
    run_dir: Path,
) -> Path:
    source = inside(package_root, str(construction["gromacs_input_dir"]))
    if not source.is_dir():
        raise ValueError(f"gromacs_input_dir is not a directory: {source}")
    common_run = inside(package_root, str(construction["common_minimization_run_directory"]))
    common_coordinates = common_run / "work" / "step6.0_minimization.gro"
    common_record_path = common_run / "common_minimization_record.json"
    if not common_coordinates.is_file() or common_coordinates.stat().st_size == 0:
        raise FileNotFoundError(f"Passing common minimization output is required: {common_coordinates}")
    if not common_record_path.is_file():
        raise FileNotFoundError(f"Common minimization record is required: {common_record_path}")
    common_record = load_json(common_record_path)
    common_errors = validate_common_minimization_record(
        package_root, manifest, construction, common_record_path
    )
    if common_errors or common_record.get("coordinates_sha256") != sha256(common_coordinates):
        raise RuntimeError("Common minimization record does not validate the shared coordinates: " + "; ".join(common_errors))
    work = run_dir / "work"
    source_record = run_dir / "staged_source.json"
    source_identity = {
        "construction_id": construction["id"],
        "source": str(source),
        "archive_sha256": construction["charmm_gui_archive"]["sha256"],
        "gromacs_input_tree_manifest_sha256": construction["gromacs_input_tree_manifest"]["sha256"],
        "common_minimization_sha256": sha256(common_coordinates),
    }
    if work.exists():
        if not source_record.is_file():
            raise RuntimeError(f"Existing work directory has no provenance record: {work}")
        recorded = load_json(source_record)
        for key, value in source_identity.items():
            if recorded.get(key) != value:
                raise RuntimeError(f"Existing work directory has different construction provenance: {work}")
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, work)
        source_record.write_text(json.dumps({**source_identity, "staged_at_utc": utc_now()}, indent=2) + "\n", encoding="utf-8")
    production_source = artifact_path(package_root, construction["production_mdp"], "build01 production MDP")
    production_target = work / "production.frozen.mdp"
    if production_target.exists() and sha256(production_target) != sha256(production_source):
        raise RuntimeError(f"Refusing to replace a changed frozen production MDP: {production_target}")
    if not production_target.exists():
        shutil.copy2(production_source, production_target)
    fork_target = work / "step6.0_minimization.gro"
    if fork_target.exists() and sha256(fork_target) != sha256(common_coordinates):
        raise RuntimeError(f"Realization does not use the frozen common minimization: {fork_target}")
    if not fork_target.exists():
        shutil.copy2(common_coordinates, fork_target)
    validate_staged_mdp_hashes(source, construction)
    seed_record = patch_velocity_seed(
        work / "step6.1_equilibration.mdp",
        work / "step6.1_equilibration.frozen.mdp",
        int(realization["velocity_seed"]),
    )
    for stage in EQUILIBRATION_STAGES[1:]:
        settings = parse_mdp(work / f"{stage}.mdp")
        if settings.get("gen_vel", "no").strip().lower() in {"yes", "true"}:
            raise ValueError(f"{stage}.mdp must continue the first-NVT velocities without regeneration")
        if settings.get("continuation", "no").strip().lower() not in {"yes", "true"}:
            raise ValueError(f"{stage}.mdp must set continuation=yes and consume the previous checkpoint")
    (run_dir / "velocity_seed_record.json").write_text(json.dumps(seed_record, indent=2) + "\n", encoding="utf-8")
    return work


def _next_command_sequence(log_dir: Path) -> tuple[int, str | None]:
    finished = sorted(log_dir.glob("[0-9][0-9][0-9][0-9]_*.finished.json"))
    previous_hash = sha256(finished[-1]) if finished else None
    occupied: list[int] = []
    for path in log_dir.glob("[0-9][0-9][0-9][0-9]_*"):
        try:
            occupied.append(int(path.name.split("_", 1)[0]))
        except ValueError:
            continue
    return (max(occupied, default=0) + 1, previous_hash)


def _write_new_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def run_command(argv: list[str], work: Path, log_dir: Path, extra_mdrun: list[str]) -> Path:
    command = list(argv)
    if len(command) > 1 and command[1] == "mdrun":
        command.extend(extra_mdrun)
        if "-cpi" not in command and "-deffnm" in command:
            deffnm = command[command.index("-deffnm") + 1]
            checkpoint = work / f"{deffnm}.cpt"
            if checkpoint.exists():
                command.extend(["-cpi", checkpoint.name, "-append"])
    token = command[1]
    if "-deffnm" in command:
        token += "_" + command[command.index("-deffnm") + 1]
        if "-s" in command:
            token += "_" + Path(command[command.index("-s") + 1]).stem
    elif "-o" in command:
        token += "_" + Path(command[command.index("-o") + 1]).stem
    sequence, previous_record_sha256 = _next_command_sequence(log_dir)
    safe_token = re.sub(r"[^A-Za-z0-9_.-]+", "_", token)
    prefix = f"{sequence:04d}_{safe_token}"
    started_path = log_dir / f"{prefix}.started.json"
    finished_path = log_dir / f"{prefix}.finished.json"
    stdout_path = log_dir / f"{prefix}.stdout.log"
    stderr_path = log_dir / f"{prefix}.stderr.log"
    started: dict[str, Any] = {
        "schema_version": "1.0",
        "sequence": sequence,
        "argv": command,
        "cwd": str(work),
        "started_at_utc": utc_now(),
        "previous_finished_record_sha256": previous_record_sha256,
    }
    checkpoint_output_path: Path | None = None
    if len(command) > 1 and command[1] == "mdrun":
        if "-s" in command:
            tpr_path = (work / command[command.index("-s") + 1]).resolve()
            if not tpr_path.is_file() or tpr_path.stat().st_size == 0:
                raise FileNotFoundError(f"mdrun TPR is missing or empty: {tpr_path}")
            started["tpr_input"] = {
                "path": str(tpr_path), "bytes": tpr_path.stat().st_size, "sha256": sha256(tpr_path)
            }
        if "-cpi" in command:
            checkpoint_input_path = (work / command[command.index("-cpi") + 1]).resolve()
            if not checkpoint_input_path.is_file() or checkpoint_input_path.stat().st_size == 0:
                raise FileNotFoundError(f"mdrun checkpoint input is missing or empty: {checkpoint_input_path}")
            started["checkpoint_input"] = {
                "path": str(checkpoint_input_path),
                "bytes": checkpoint_input_path.stat().st_size,
                "sha256": sha256(checkpoint_input_path),
            }
        if "-deffnm" in command:
            checkpoint_output_path = (work / f"{command[command.index('-deffnm') + 1]}.cpt").resolve()
    _write_new_json(started_path, started)
    with stdout_path.open("x", encoding="utf-8") as stdout, stderr_path.open("x", encoding="utf-8") as stderr:
        result = subprocess.run(command, cwd=work, stdout=stdout, stderr=stderr, check=False)
    finished = {
        "schema_version": "1.0",
        "sequence": sequence,
        "started_record": str(started_path),
        "started_record_sha256": sha256(started_path),
        "stdout": {"path": str(stdout_path), "sha256": sha256(stdout_path)},
        "stderr": {"path": str(stderr_path), "sha256": sha256(stderr_path)},
        "finished_at_utc": utc_now(),
        "returncode": result.returncode,
    }
    if checkpoint_output_path is not None and checkpoint_output_path.is_file() and checkpoint_output_path.stat().st_size > 0:
        finished["checkpoint_output"] = {
            "path": str(checkpoint_output_path),
            "bytes": checkpoint_output_path.stat().st_size,
            "sha256": sha256(checkpoint_output_path),
        }
    if len(command) > 1 and command[1] == "grompp" and "-o" in command:
        declared_output_path = (work / command[command.index("-o") + 1]).resolve()
        if declared_output_path.is_file() and declared_output_path.stat().st_size > 0:
            finished["declared_output"] = {
                "path": str(declared_output_path),
                "bytes": declared_output_path.stat().st_size,
                "sha256": sha256(declared_output_path),
            }
    if len(command) > 1 and command[1] == "mdrun" and "-deffnm" in command:
        deffnm = command[command.index("-deffnm") + 1]
        runtime_outputs: dict[str, dict[str, Any]] = {}
        for suffix in ("xtc", "edr", "log", "gro", "cpt"):
            output_path = (work / f"{deffnm}.{suffix}").resolve()
            if output_path.is_file() and output_path.stat().st_size > 0:
                runtime_outputs[suffix] = {
                    "path": str(output_path),
                    "bytes": output_path.stat().st_size,
                    "sha256": sha256(output_path),
                }
        finished["runtime_outputs"] = runtime_outputs
    _write_new_json(finished_path, finished)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}); see {finished_path}")
    return finished_path


def _successful_command_records(log_dir: Path) -> list[tuple[dict[str, Any], dict[str, Any], Path]]:
    records: list[tuple[dict[str, Any], dict[str, Any], Path]] = []
    previous_hash: str | None = None
    for finished_path in sorted(log_dir.glob("[0-9][0-9][0-9][0-9]_*.finished.json")):
        finished = load_json(finished_path)
        started_path = Path(str(finished.get("started_record", ""))).resolve()
        if not started_path.is_file() or finished.get("started_record_sha256") != sha256(started_path):
            raise RuntimeError(f"Invalid command provenance: {finished_path}")
        started = load_json(started_path)
        if started.get("previous_finished_record_sha256") != previous_hash:
            raise RuntimeError(f"Broken append-only command chain: {finished_path}")
        previous_hash = sha256(finished_path)
        if finished.get("returncode") == 0:
            records.append((started, finished, finished_path))
    return records


def _validate_existing_grompp_tpr(log_dir: Path, tpr: Path) -> bool:
    matches = []
    for started, finished, finished_path in _successful_command_records(log_dir):
        argv = started.get("argv")
        if not isinstance(argv, list) or len(argv) < 2 or argv[1] != "grompp" or "-o" not in argv:
            continue
        if Path(argv[argv.index("-o") + 1]).name == tpr.name:
            matches.append((finished, finished_path))
    if not matches:
        return False
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one successful grompp for {tpr.name}, found {len(matches)}")
    output = matches[0][0].get("declared_output", {})
    if not isinstance(output, dict) or output.get("sha256") != sha256(tpr):
        raise RuntimeError(f"Existing {tpr.name} is not bound to its successful grompp record")
    return True


def run_equilibration_restart_safe(
    gmx: str,
    work: Path,
    log_dir: Path,
    extra_mdrun: list[str],
) -> None:
    """Advance six stages without regenerating an existing TPR or overwriting partial data."""
    previous = "step6.0_minimization"
    for stage in EQUILIBRATION_STAGES:
        mdp = "step6.1_equilibration.frozen.mdp" if stage == EQUILIBRATION_STAGES[0] else f"{stage}.mdp"
        tpr = work / f"{stage}.tpr"
        log = work / f"{stage}.log"
        cpt = work / f"{stage}.cpt"
        gro = work / f"{stage}.gro"
        edr = work / f"{stage}.edr"
        complete = all(path.is_file() and path.stat().st_size > 0 for path in (tpr, log, cpt, gro, edr))
        if complete:
            if not _validate_existing_grompp_tpr(log_dir, tpr):
                raise RuntimeError(f"Completed {stage} lacks one immutable successful grompp record")
            log_text = log.read_text(encoding="utf-8", errors="replace")
            if "Finished mdrun" not in log_text or any(pattern.search(log_text) for pattern in FORBIDDEN_RUNTIME_PATTERNS):
                raise RuntimeError(f"Completed {stage} does not pass retained-log integrity")
            previous = stage
            continue
        if not tpr.exists():
            partial = [path for path in (log, cpt, gro, edr) if path.exists()]
            if partial:
                raise RuntimeError(f"{stage} has partial outputs but no TPR; refusing destructive restart: {partial}")
            grompp = [gmx, "grompp", "-f", mdp, "-o", tpr.name, "-c", f"{previous}.gro"]
            if stage != EQUILIBRATION_STAGES[0]:
                prior_cpt = work / f"{previous}.cpt"
                if not prior_cpt.is_file() or prior_cpt.stat().st_size == 0:
                    raise RuntimeError(f"{stage} cannot start without {prior_cpt.name}")
                grompp.extend(["-t", prior_cpt.name])
            grompp.extend(["-r", "step5_input.gro", "-p", "topol.top", "-n", "index.ndx"])
            run_command(grompp, work, log_dir, [])
        elif not _validate_existing_grompp_tpr(log_dir, tpr):
            raise RuntimeError(f"Existing {tpr.name} lacks immutable successful grompp provenance")
        mdrun = [gmx, "mdrun", "-s", tpr.name, "-deffnm", stage]
        if cpt.is_file() and cpt.stat().st_size > 0:
            mdrun.extend(["-cpi", cpt.name, "-append"])
        elif any(path.exists() for path in (log, gro, edr)):
            raise RuntimeError(f"{stage} has partial outputs without a usable checkpoint; manual quarantine is required")
        run_command(mdrun, work, log_dir, extra_mdrun)
        previous = stage


def make_common_minimization_record(
    manifest: dict[str, Any],
    construction: dict[str, Any],
    work: Path,
    command_records: list[Path],
) -> dict[str, Any]:
    """Create a checksum-sealed, fail-closed technical audit of minimization."""
    required = {
        "mdp": work / "step6.0_minimization.mdp",
        "tpr": work / "step6.0_minimization.tpr",
        "energy": work / "step6.0_minimization.edr",
        "log": work / "step6.0_minimization.log",
        "coordinates": work / "step6.0_minimization.gro",
    }
    errors: list[str] = []
    artifacts: dict[str, dict[str, Any]] = {}
    for label, path in required.items():
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty minimization {label}: {path}")
        else:
            artifacts[label] = {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
    if required["mdp"].is_file() and sha256(required["mdp"]) != construction.get("minimization_mdp", {}).get("sha256"):
        errors.append("executed minimization MDP differs from the frozen build input")

    retained_commands: list[dict[str, Any]] = []
    command_kinds: list[str] = []
    text_paths: list[Path] = [required["log"]] if required["log"].is_file() else []
    previous_finished_hash: str | None = None
    for record_path in command_records:
        if not record_path.is_file() or record_path.stat().st_size == 0:
            errors.append(f"missing minimization command record: {record_path}")
            continue
        finished = load_json(record_path)
        started_path = Path(str(finished.get("started_record", ""))).resolve()
        if not started_path.is_file() or finished.get("started_record_sha256") != sha256(started_path):
            errors.append(f"minimization command started-record provenance is invalid: {record_path}")
            continue
        started = load_json(started_path)
        if started.get("previous_finished_record_sha256") != previous_finished_hash:
            errors.append(f"minimization command hash chain is broken: {record_path}")
        argv = started.get("argv")
        if not isinstance(argv, list) or len(argv) < 2 or not all(isinstance(item, str) for item in argv):
            errors.append(f"invalid minimization command argv: {started_path}")
        else:
            command_kinds.append(argv[1])
            if "-maxwarn" in argv:
                errors.append(f"-maxwarn is prohibited: {started_path}")
        if finished.get("returncode") != 0:
            errors.append(f"minimization command failed: {record_path}")
        for stream in ("stdout", "stderr"):
            stream_artifact = finished.get(stream)
            stream_path = Path(str(stream_artifact.get("path", ""))).resolve() if isinstance(stream_artifact, dict) else Path()
            if (
                not isinstance(stream_artifact, dict) or not stream_path.is_file()
                or stream_artifact.get("sha256") != sha256(stream_path)
            ):
                errors.append(f"minimization {stream} provenance is invalid: {record_path}")
            else:
                text_paths.append(stream_path)
        retained_commands.append({"path": str(record_path), "sha256": sha256(record_path)})
        previous_finished_hash = sha256(record_path)
    if command_kinds.count("grompp") != 1 or command_kinds.count("mdrun") < 1:
        errors.append("common minimization requires exactly one retained grompp and at least one retained mdrun segment")

    combined = "\n".join(path.read_text(encoding="utf-8", errors="replace") for path in text_paths)
    warning_count = len(GROMPP_WARNING.findall(combined))
    if warning_count:
        errors.append(f"grompp reported {warning_count} numbered warning(s)")
    forbidden_count = sum(len(pattern.findall(combined)) for pattern in FORBIDDEN_RUNTIME_PATTERNS)
    if forbidden_count:
        errors.append(f"minimization logs contain {forbidden_count} forbidden runtime pattern(s)")
    if "Finished mdrun" not in combined:
        errors.append("minimization logs lack the Finished mdrun marker")

    coordinates = required["coordinates"]
    report: dict[str, Any] = {
        "schema_version": "2.0",
        "report_type": "common_minimization_validation",
        "status": "pass" if not errors else "fail",
        "technical_integrity_pass": not errors,
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "build_contract_sha256": build_contract_sha256(manifest),
        "archive_sha256": construction["charmm_gui_archive"]["sha256"],
        "gromacs_input_tree_manifest_sha256": construction["gromacs_input_tree_manifest"]["sha256"],
        "minimization_mdp_sha256": construction["minimization_mdp"]["sha256"],
        "coordinates": str(coordinates),
        "coordinates_sha256": sha256(coordinates) if coordinates.is_file() else None,
        "artifacts": artifacts,
        "command_records": retained_commands,
        "grompp_warning_count": warning_count,
        "forbidden_runtime_pattern_count": forbidden_count,
        "errors": errors,
        "completed_at_utc": utc_now(),
        "integrity": {"payload_sha256": "UNSEALED"},
    }
    report["integrity"]["payload_sha256"] = report_payload_sha256(report, ("integrity", "payload_sha256"))
    return report


def freeze_or_validate_production_tpr(
    manifest: dict[str, Any],
    manifest_path: Path,
    construction: dict[str, Any],
    realization_id: str,
    work: Path,
    log_dir: Path,
    gmx: str,
) -> dict[str, Any]:
    """Create the 500-ns TPR once; every later segment must reuse its exact hash."""
    tpr = work / "production.tpr"
    record_path = work.parent / "production_tpr_record.json"
    mdp = work / "production.frozen.mdp"
    final_gro = work / "step6.6_equilibration.gro"
    final_cpt = work / "step6.6_equilibration.cpt"
    for path in (mdp, final_gro, final_cpt, work / "topol.top", work / "index.ndx"):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Cannot freeze production TPR; missing {path}")
    identity = {
        "schema_version": "1.0",
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "realization_id": realization_id,
        "build_contract_sha256": build_contract_sha256(manifest),
        "production_mdp_sha256": sha256(mdp),
        "final_equilibration_coordinates_sha256": sha256(final_gro),
        "final_equilibration_checkpoint_sha256": sha256(final_cpt),
        "topology_sha256": sha256(work / "topol.top"),
        "index_sha256": sha256(work / "index.ndx"),
    }
    if tpr.exists() or record_path.exists():
        if not tpr.is_file() or not record_path.is_file():
            raise RuntimeError("production.tpr and production_tpr_record.json must exist together")
        recorded = load_json(record_path)
        for key, value in identity.items():
            if recorded.get(key) != value:
                raise RuntimeError(f"Existing production TPR provenance differs at {key}")
        if recorded.get("production_tpr_sha256") != sha256(tpr):
            raise RuntimeError("Existing production.tpr has changed; grompp will not be rerun")
        return recorded

    command_record = run_command(
        [gmx, "grompp", "-f", "production.frozen.mdp", "-o", "production.tpr", "-c", "step6.6_equilibration.gro", "-t", "step6.6_equilibration.cpt", "-p", "topol.top", "-n", "index.ndx"],
        work,
        log_dir,
        [],
    )
    if not tpr.is_file() or tpr.stat().st_size == 0:
        raise RuntimeError("grompp did not create production.tpr")
    record = {
        **identity,
        "production_tpr": str(tpr),
        "production_tpr_sha256": sha256(tpr),
        "grompp_finished_record": str(command_record),
        "grompp_finished_record_sha256": sha256(command_record),
        "frozen_at_utc": utc_now(),
    }
    _write_new_json(record_path, record)
    return record


def validate_production_tpr_record(
    package_root: Path, manifest: dict[str, Any], realization_id: str, work: Path
) -> dict[str, Any]:
    record_path = work.parent / "production_tpr_record.json"
    tpr = work / "production.tpr"
    if not record_path.is_file() or not tpr.is_file():
        raise FileNotFoundError("Continuation requires the frozen production TPR and its provenance record")
    record = load_json(record_path)
    expected = {
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "realization_id": realization_id,
        "build_contract_sha256": build_contract_sha256(manifest),
        "production_tpr_sha256": sha256(tpr),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise RuntimeError(f"Production TPR record mismatch at {key}")
    checkpoint = work / "production.cpt"
    if not checkpoint.is_file() or checkpoint.stat().st_size == 0:
        raise FileNotFoundError("Continuation requires the retained production canary checkpoint")
    construction = manifest["systems"][0]["construction"]
    canary_report_path = artifact_path(
        package_root, construction.get("canary_validation_report"), "strict canary validation report"
    )
    canary_report = load_json(canary_report_path)
    if (
        canary_report.get("schema_version") != "2.0"
        or canary_report.get("report_type") != "md_stage_output_validation"
        or canary_report.get("phase") != "canary"
        or canary_report.get("strict") is not True
        or canary_report.get("status") != "pass"
    ):
        raise RuntimeError("Continuation requires the bound strict-PASS canary report")
    integrity = canary_report.get("integrity", {})
    if not isinstance(integrity, dict) or integrity.get("payload_sha256") != report_payload_sha256(
        canary_report, ("integrity", "payload_sha256")
    ):
        raise RuntimeError("Canary report payload checksum is invalid")
    matching = [
        item for item in canary_report.get("runs", [])
        if isinstance(item, dict) and item.get("realization_id") == realization_id
    ]
    if len(matching) != 1:
        raise RuntimeError(f"Canary report lacks one unique run for {realization_id}")
    canary_run = matching[0]
    checkpoint_record = canary_run.get("artifacts", {}).get("checkpoint", {})
    expected_checkpoint_sha256 = checkpoint_record.get("sha256") if isinstance(checkpoint_record, dict) else None
    command_logs = work.parent / "command_logs"
    expected_current_sha256 = expected_checkpoint_sha256
    previous_finished_sha256: str | None = None
    for finished_path in sorted(command_logs.glob("[0-9][0-9][0-9][0-9]_*.finished.json")):
        finished = load_json(finished_path)
        started_path = Path(str(finished.get("started_record", ""))).resolve()
        if not started_path.is_file() or finished.get("started_record_sha256") != sha256(started_path):
            raise RuntimeError(f"Command provenance is invalid before continuation: {finished_path}")
        started = load_json(started_path)
        if started.get("previous_finished_record_sha256") != previous_finished_sha256:
            raise RuntimeError(f"Command hash chain is broken before continuation: {finished_path}")
        previous_finished_sha256 = sha256(finished_path)
        argv = started.get("argv")
        if not isinstance(argv, list) or len(argv) < 2:
            raise RuntimeError(f"Command argv is invalid before continuation: {started_path}")
        is_production_continuation = (
            argv[1] == "mdrun" and "production.tpr" in argv and "-cpi" in argv
            and "-append" in argv and "-nsteps" not in argv
        )
        if not is_production_continuation:
            continue
        if finished.get("returncode") != 0:
            raise RuntimeError(f"A prior production continuation failed: {finished_path}")
        checkpoint_input = started.get("checkpoint_input", {})
        if not isinstance(checkpoint_input, dict) or checkpoint_input.get("sha256") != expected_current_sha256:
            raise RuntimeError(f"Production checkpoint input chain is broken for {realization_id}")
        checkpoint_output = finished.get("checkpoint_output", {})
        next_sha256 = checkpoint_output.get("sha256") if isinstance(checkpoint_output, dict) else None
        if not isinstance(next_sha256, str) or len(next_sha256) != 64:
            raise RuntimeError(f"Production continuation lacks a sealed checkpoint output: {finished_path}")
        expected_current_sha256 = next_sha256
    if expected_current_sha256 != sha256(checkpoint):
        raise RuntimeError(f"Production checkpoint differs from its latest approved chain state for {realization_id}")
    if canary_run.get("final_time_ps") != 5000.0:
        raise RuntimeError(f"Canary report endpoint is not exactly 5 ns for {realization_id}")
    return record


def get_contract(manifest: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    systems = manifest.get("systems")
    if not isinstance(systems, list) or len(systems) != 1 or systems[0].get("id") != SYSTEM_ID:
        raise ValueError(f"Manifest must contain exactly one system: {SYSTEM_ID}")
    system = systems[0]
    construction = system.get("construction")
    realizations = system.get("realizations")
    if not isinstance(construction, dict) or construction.get("id") != "build01":
        raise ValueError("Manifest must contain exactly one construction named build01")
    if not isinstance(realizations, list) or [item.get("id") for item in realizations] != list(REALIZATION_IDS):
        raise ValueError(f"Manifest realizations must be exactly {list(REALIZATION_IDS)}")
    if len({item.get("velocity_seed") for item in realizations}) != 3:
        raise ValueError("All three velocity seeds must be distinct")
    return system, construction, realizations


def synthetic_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="single_system_md_runner_") as temporary:
        root = Path(temporary)
        config = root / "config"
        source = root / "builds" / "build01" / "gromacs"
        config.mkdir(parents=True)
        source.mkdir(parents=True)
        (source / "step6.0_minimization.mdp").write_text("integrator = steep\ngen_vel = no\n", encoding="utf-8")
        (source / "step6.1_equilibration.mdp").write_text(
            "integrator = md\ntcoupl = v-rescale\npcoupl = no\ngen_vel = yes\ngen_seed = -1\n",
            encoding="utf-8",
        )
        for index in range(2, 7):
            (source / f"step6.{index}_equilibration.mdp").write_text("integrator = md\ngen_vel = no\ncontinuation = yes\n", encoding="utf-8")
        for name in ("step5_input.gro", "topol.top", "index.ndx"):
            (source / name).write_text("synthetic\n", encoding="utf-8")
        archive = root / "builds" / "build01" / "archive.tgz"
        archive.write_bytes(b"synthetic archive")
        production = root / "config" / "production_500ns.mdp"
        production.write_text("""integrator = md
dt = 0.002
nsteps = 250000000
continuation = yes
gen_vel = no
pbc = xyz
periodic-molecules = no
mass-repartition-factor = 1.0
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
pme-order = 4
fourierspacing = 0.12
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw-switch = 1.0
rvdw = 1.2
DispCorr = no
constraints = h-bonds
constraint-algorithm = lincs
tcoupl = Nose-Hoover
tc-grps = Protein Membrane Solvent
tau-t = 1.0 1.0 1.0
ref-t = 310.15 310.15 310.15
pcoupl = Parrinello-Rahman
pcoupltype = semiisotropic
tau-p = 5.0
ref-p = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
comm-mode = linear
comm-grps = System
nstcomm = 100
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = 50000
nstcalcenergy = 100
nstenergy = 5000
nstlog = 5000
compressed-x-precision = 1000
""", encoding="utf-8")
        mdp_contract = {
            "integrator": "md", "thermostat": "nose-hoover",
            "thermostat_groups": ["Protein", "Membrane", "Solvent"], "tau_t_ps": 1.0,
            "barostat": "parrinello-rahman", "barostat_tau_p_ps": 5.0,
            "compressibility_bar_inverse": [4.5e-5, 4.5e-5], "cutoff_scheme": "verlet",
            "neighbor_list_update_steps": 20, "rlist_nm": 1.2, "rcoulomb_nm": 1.2,
            "vdw_type": "cut-off", "vdw_modifier": "force-switch", "rvdw_switch_nm": 1.0,
            "rvdw_nm": 1.2, "dispersion_correction": "no", "pme_order": 4,
            "fourier_spacing_nm": 0.12, "constraint_algorithm": "lincs",
            "com_removal_mode": "linear", "com_removal_groups": ["System"],
            "com_removal_interval_steps": 100,
            "output_cadence_steps": {
                "nstxout": 0, "nstvout": 0, "nstfout": 0, "nstxout_compressed": 50000,
                "nstcalcenergy": 100, "nstenergy": 5000, "nstlog": 5000,
                "compressed_x_precision": 1000,
            },
        }
        manifest = {
            "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
            "global_model": {"temperature_k": 310.15, "pressure_bar": 1.0},
            "simulation": {
                "production_ns": 500.0,
                "time_step_ps": 0.002,
                "hydrogen_mass_repartitioning": False,
                "constraints": "h-bonds",
                "pressure_coupling": "semiisotropic",
                "production_mdp_contract": mdp_contract,
                "release_contract": {
                    "canary_target_ns_per_realization": 5.0,
                    "canary_total_ns": 15.0,
                    "canary_uses_original_500ns_tpr": True,
                    "canary_frames_retained_as_production": True,
                },
            },
            "systems": [{
                "id": SYSTEM_ID,
                "construction": {
                    "id": "build01",
                    "gromacs_input_dir": "builds/build01/gromacs",
                    "minimization_mdp": {
                        "path": "builds/build01/gromacs/step6.0_minimization.mdp",
                        "sha256": sha256(source / "step6.0_minimization.mdp"),
                    },
                    "production_mdp": {"path": "config/production_500ns.mdp", "sha256": sha256(production)},
                    "gromacs_input_tree_manifest": {"path": "config/production_500ns.mdp", "sha256": sha256(production)},
                    "equilibration_mdp_sha256": {
                        f"step6.{index}_equilibration.mdp": sha256(source / f"step6.{index}_equilibration.mdp")
                        for index in range(1, 7)
                    },
                    "common_minimization_run_directory": f"runs/{SYSTEM_ID}/common_minimization",
                    "charmm_gui_archive": {"path": "builds/build01/archive.tgz", "sha256": sha256(archive)},
                },
                "realizations": [
                    {"id": rid, "velocity_seed": 26080801 + index, "run_directory": f"runs/{SYSTEM_ID}/{rid}"}
                    for index, rid in enumerate(REALIZATION_IDS)
                ],
            }],
        }
        manifest_path = config / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--manifest", str(manifest_path), "--phase", "equilibration"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Synthetic dry run failed:\n{result.stdout}\n{result.stderr}")
        planned = json.loads(result.stdout)
        if planned.get("mode") != "plan_only" or len(planned.get("runs", [])) != 3:
            raise RuntimeError("Synthetic dry run did not plan exactly three realizations")
        repeated_minimization = any(
            (command[1] == "grompp" and "step6.0_minimization.mdp" in command) or
            (command[1] == "mdrun" and "step6.0_minimization" in command)
            for run in planned["runs"] for command in run["commands"]
        )
        if planned.get("common_minimization") is not None or repeated_minimization:
            raise RuntimeError("Common minimization was not separated from realization equilibration commands")
        seeds = {run["velocity_seed"] for run in planned["runs"]}
        if len(seeds) != 3:
            raise RuntimeError("Synthetic dry run reused a velocity seed")
        frozen_seed_values: set[int] = set()
        for run in planned["runs"]:
            destination = root / f"{run['realization_id']}.first_nvt.mdp"
            patch_velocity_seed(source / "step6.1_equilibration.mdp", destination, int(run["velocity_seed"]))
            frozen_seed_values.add(int(parse_mdp(destination)["gen_seed"]))
        if frozen_seed_values != seeds:
            raise RuntimeError("First-NVT seed patching did not preserve all three unique seeds")
        commands = planned["runs"][0]["commands"]
        grompp_commands = [command for command in commands if command[1] == "grompp"]
        first_nvt = next(command for command in grompp_commands if "step6.1_equilibration.frozen.mdp" in command)
        if "-t" in first_nvt:
            raise RuntimeError("First NVT must generate velocities and must not consume a checkpoint")
        for stage_index in range(2, 7):
            current = next(command for command in grompp_commands if f"step6.{stage_index}_equilibration.mdp" in command)
            expected_checkpoint = f"step6.{stage_index - 1}_equilibration.cpt"
            if "-t" not in current or current[current.index("-t") + 1] != expected_checkpoint:
                raise RuntimeError(f"Equilibration stage {stage_index} does not consume {expected_checkpoint}")
        canary = realization_command_plan("gmx", "canary")
        production_grompp = canary[0]
        if production_grompp[production_grompp.index("-t") + 1] != "step6.6_equilibration.cpt":
            raise RuntimeError("Canary production grompp does not consume the final equilibration checkpoint")
        if canary[1][-2:] != ["-nsteps", "2500000"]:
            raise RuntimeError("Canary is not frozen to 5 ns at 2 fs")
        continuation = realization_command_plan("gmx", "production")
        if len(continuation) != 1 or "grompp" in continuation[0] or "-cpi" not in continuation[0] or "-append" not in continuation[0]:
            raise RuntimeError("Production continuation can rerun grompp or does not preserve its checkpoint")
        validate_frozen_production_mdp(production, manifest)
        tampered = config / "wrong_physics.mdp"
        tampered.write_text(production.read_text(encoding="utf-8").replace(
            "ref-t = 310.15 310.15 310.15", "ref-t = 100 100 100"
        ), encoding="utf-8")
        try:
            validate_frozen_production_mdp(tampered, manifest)
        except ValueError:
            pass
        else:
            raise RuntimeError("Wrong production temperature was accepted")
        hmr = config / "wrong_hmr.mdp"
        hmr.write_text(production.read_text(encoding="utf-8").replace(
            "mass-repartition-factor = 1.0", "mass-repartition-factor = 3.0"
        ), encoding="utf-8")
        try:
            validate_frozen_production_mdp(hmr, manifest)
        except ValueError:
            pass
        else:
            raise RuntimeError("Hydrogen mass repartitioning was accepted")
        biased = config / "wrong_bias.mdp"
        biased.write_text(production.read_text(encoding="utf-8") + "pull = yes\n", encoding="utf-8")
        try:
            validate_frozen_production_mdp(biased, manifest)
        except ValueError:
            pass
        else:
            raise RuntimeError("A biased production MDP was accepted")

        # An interrupted equilibration must reuse its one immutable TPR; a
        # changed or repeatedly generated TPR must fail before mdrun resumes.
        eq_test = root / "equilibration_restart_test"
        eq_test.mkdir()
        eq_tpr = eq_test / "step6.3_equilibration.tpr"
        eq_tpr.write_bytes(b"synthetic immutable equilibration tpr")
        eq_started = eq_test / "0001_grompp_step6.3.started.json"
        eq_started.write_text(json.dumps({
            "sequence": 1,
            "argv": ["gmx", "grompp", "-f", "step6.3_equilibration.mdp", "-o", eq_tpr.name],
            "previous_finished_record_sha256": None,
        }), encoding="utf-8")
        eq_finished = eq_test / "0001_grompp_step6.3.finished.json"
        eq_finished.write_text(json.dumps({
            "sequence": 1,
            "started_record": str(eq_started),
            "started_record_sha256": sha256(eq_started),
            "returncode": 0,
            "declared_output": {"path": str(eq_tpr), "bytes": eq_tpr.stat().st_size, "sha256": sha256(eq_tpr)},
        }), encoding="utf-8")
        if not _validate_existing_grompp_tpr(eq_test, eq_tpr):
            raise RuntimeError("Immutable equilibration TPR provenance was not recognized")
        eq_tpr.write_bytes(b"tampered equilibration tpr")
        try:
            _validate_existing_grompp_tpr(eq_test, eq_tpr)
        except RuntimeError:
            pass
        else:
            raise RuntimeError("Tampered equilibration TPR was accepted for restart")

        # Continuation must be bound to the exact checkpoint sealed at 5 ns.
        test_run = root / "runs" / SYSTEM_ID / "rep01"
        test_work = test_run / "work"
        test_work.mkdir(parents=True, exist_ok=True)
        test_tpr = test_work / "production.tpr"
        test_cpt = test_work / "production.cpt"
        test_tpr.write_bytes(b"synthetic frozen production tpr")
        test_cpt.write_bytes(b"synthetic approved 5 ns checkpoint")
        canary_path = root / "reports" / "canary.json"
        canary_path.parent.mkdir(parents=True)
        manifest["systems"][0]["construction"]["canary_validation_report"] = {
            "path": str(canary_path.relative_to(root)), "sha256": "0" * 64
        }
        build_hash = build_contract_sha256(manifest)
        (test_run / "production_tpr_record.json").write_text(json.dumps({
            "system_id": SYSTEM_ID,
            "construction_id": "build01",
            "realization_id": "rep01",
            "build_contract_sha256": build_hash,
            "production_tpr_sha256": sha256(test_tpr),
        }), encoding="utf-8")
        canary_report = {
            "schema_version": "2.0",
            "report_type": "md_stage_output_validation",
            "phase": "canary",
            "strict": True,
            "status": "pass",
            "runs": [{
                "realization_id": "rep01",
                "final_time_ps": 5000.0,
                "artifacts": {"checkpoint": {"sha256": sha256(test_cpt)}},
            }],
            "integrity": {"payload_sha256": "UNSEALED"},
        }
        canary_report["integrity"]["payload_sha256"] = report_payload_sha256(
            canary_report, ("integrity", "payload_sha256")
        )
        canary_path.write_text(json.dumps(canary_report, indent=2) + "\n", encoding="utf-8")
        manifest["systems"][0]["construction"]["canary_validation_report"]["sha256"] = sha256(canary_path)
        validate_production_tpr_record(root, manifest, "rep01", test_work)
        test_cpt.write_bytes(b"tampered checkpoint after canary approval")
        try:
            validate_production_tpr_record(root, manifest, "rep01", test_work)
        except RuntimeError:
            pass
        else:
            raise RuntimeError("A checkpoint changed after canary approval was accepted")
        test_cpt.write_bytes(b"synthetic approved 5 ns checkpoint")
        command_logs = test_run / "command_logs"
        command_logs.mkdir()
        started_path = command_logs / "0001_mdrun_production.started.json"
        started_path.write_text(json.dumps({
            "sequence": 1,
            "argv": ["gmx", "mdrun", "-s", "production.tpr", "-deffnm", "production", "-cpi", "production.cpt", "-append"],
            "previous_finished_record_sha256": None,
            "checkpoint_input": {"path": str(test_cpt), "bytes": test_cpt.stat().st_size, "sha256": sha256(test_cpt)},
        }), encoding="utf-8")
        test_cpt.write_bytes(b"synthetic approved later-segment checkpoint")
        finished_path = command_logs / "0001_mdrun_production.finished.json"
        finished_path.write_text(json.dumps({
            "sequence": 1,
            "started_record": str(started_path),
            "started_record_sha256": sha256(started_path),
            "returncode": 0,
            "checkpoint_output": {"path": str(test_cpt), "bytes": test_cpt.stat().st_size, "sha256": sha256(test_cpt)},
        }), encoding="utf-8")
        validate_production_tpr_record(root, manifest, "rep01", test_work)
        test_cpt.write_bytes(b"tampered after a later continuation segment")
        try:
            validate_production_tpr_record(root, manifest, "rep01", test_work)
        except RuntimeError:
            pass
        else:
            raise RuntimeError("A checkpoint changed after a later continuation was accepted")
    validate_mdrun_extra(["-ntomp", "8", "-gpu_id", "0", "-maxh", "23.5"])
    try:
        validate_mdrun_extra(["-deffnm", "changed"])
    except ValueError:
        pass
    else:
        raise RuntimeError("Scientific run-identity override was not rejected")
    print("SELF-TEST PASS: independent restart-safe equilibration, strict frozen production physics, exact 5 ns same-TPR canary, no-grompp checkpoint continuation, unique seeds, and performance-only mdrun extras.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--phase",
        choices=("minimization", "equilibration", "canary", "production"),
        default="minimization",
        help="production means checkpoint continuation after the mandatory all-three 5 ns canary gate",
    )
    parser.add_argument("--realization", action="append", choices=REALIZATION_IDS, help="Repeat to select realization jobs")
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--execute", action="store_true", help="Execute after the requested stage's strict upstream preflight; default is plan only")
    parser.add_argument("--mdrun-extra", action="append", default=[], help="One literal extra mdrun argument per use")
    parser.add_argument("--self-test", action="store_true", help="Run a deterministic synthetic planning test")
    args = parser.parse_args()

    validate_mdrun_extra(args.mdrun_extra)

    if args.self_test:
        if args.manifest is not None or args.execute or args.realization:
            raise SystemExit("--self-test cannot be combined with manifest or execution options")
        synthetic_self_test()
        return 0
    if args.manifest is None:
        parser.error("--manifest is required unless --self-test is used")

    manifest_path = args.manifest.resolve()
    package_root = manifest_path.parent.parent.resolve()
    manifest = load_json(manifest_path)
    _, construction, realizations = get_contract(manifest)
    production_ns = float(manifest["simulation"]["production_ns"])
    if abs(production_ns - 500.0) > 1e-12:
        raise SystemExit("production_ns must be exactly 500 ns")
    production_mdp = artifact_path(package_root, construction["production_mdp"], "build01 production MDP")
    validate_frozen_production_mdp(production_mdp, manifest)
    release = manifest["simulation"].get("release_contract", {})
    if release.get("canary_target_ns_per_realization") != 5.0 or release.get("canary_total_ns") != 15.0:
        raise SystemExit("The frozen canary must be exactly 5.0 ns per realization and 15.0 ns total")
    if release.get("canary_uses_original_500ns_tpr") is not True or release.get("canary_frames_retained_as_production") is not True:
        raise SystemExit("Canary must use the original 500-ns TPR and retain its frames as production")
    canary_steps = int(5.0 * 1000.0 / float(manifest["simulation"]["time_step_ps"]))
    if canary_steps != 2_500_000:
        raise SystemExit("The fixed 5 ns canary must equal 2,500,000 steps at 2 fs")

    # Canary/full-continuation command plans are release artifacts.  They are
    # withheld even in plan-only mode until their respective upstream gate is
    # strict-pass.  Earlier stages are checked before execution.
    preflight_stage = {
        "minimization": "builds",
        "equilibration": "equilibration",
        "canary": "canary",
        "production": "production",
    }[args.phase]
    if args.execute or args.phase in {"canary", "production"}:
        preflight = package_root / "scripts" / "validate_preflight.py"
        result = subprocess.run(
            [sys.executable, str(preflight), "--manifest", str(manifest_path), "--stage", preflight_stage, "--strict"],
            cwd=package_root,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            detail = (result.stdout + "\n" + result.stderr).strip()
            raise SystemExit(f"Strict {preflight_stage} preflight failed; no command was released:\n{detail}")
        try:
            validate_bound_gmx_identity(package_root, manifest, args.gmx)
        except (OSError, ValueError, RuntimeError) as exc:
            raise SystemExit(f"Bound GROMACS executable validation failed; no command was released: {exc}")

    if args.phase == "minimization":
        if args.realization:
            raise SystemExit("The common minimization has no realization selector")
        selected = set()
        commands = []
        common_commands = common_minimization_plan(args.gmx)
    else:
        selected = set(args.realization or REALIZATION_IDS)
        commands = realization_command_plan(args.gmx, args.phase, canary_steps)
        common_commands = []

    planned: list[dict[str, Any]] = []
    common_run_dir = inside(package_root, str(construction["common_minimization_run_directory"]), must_exist=False)
    common_plan = {
        "construction_id": "build01",
        "stage": "one_common_deterministic_minimization",
        "run_directory": str(common_run_dir),
        "commands": common_commands,
    } if common_commands else None
    if args.execute and common_commands:
        common_run_dir, common_work = stage_common_minimization(package_root, construction)
        common_record_path = common_run_dir / "common_minimization_record.json"
        if common_record_path.exists():
            raise RuntimeError("Common-minimization record already exists; refusing to rerun or overwrite it")
        common_logs = common_run_dir / "command_logs"
        common_logs.mkdir(parents=True, exist_ok=True)
        for command in common_commands:
            run_command(command, common_work, common_logs, args.mdrun_extra)
        all_common_command_records = sorted(common_logs.glob("[0-9][0-9][0-9][0-9]_*.finished.json"))
        common_report = make_common_minimization_record(
            manifest, construction, common_work, all_common_command_records
        )
        _write_new_json(common_record_path, common_report)
        if common_report["status"] != "pass":
            raise RuntimeError(f"Common minimization failed technical validation; retained report: {common_record_path}")
    for realization in realizations:
        rid = str(realization["id"])
        if rid not in selected:
            continue
        run_dir = inside(package_root, str(realization["run_directory"]), must_exist=False)
        planned.append({
            "system_id": SYSTEM_ID,
            "construction_id": "build01",
            "realization_id": rid,
            "velocity_seed": realization["velocity_seed"],
            "run_directory": str(run_dir),
            "commands": commands,
        })
        if args.execute:
            work = stage_inputs(package_root, manifest, construction, realization, run_dir)
            log_dir = run_dir / "command_logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            if args.phase == "canary":
                freeze_or_validate_production_tpr(
                    manifest, manifest_path, construction, rid, work, log_dir, args.gmx
                )
                canary_command = [
                    args.gmx, "mdrun", "-s", "production.tpr", "-deffnm", "production",
                    "-nsteps", str(canary_steps),
                ]
                run_command(canary_command, work, log_dir, args.mdrun_extra)
            elif args.phase == "production":
                validate_production_tpr_record(package_root, manifest, rid, work)
                run_command(
                    [args.gmx, "mdrun", "-s", "production.tpr", "-deffnm", "production", "-cpi", "production.cpt", "-append"],
                    work,
                    log_dir,
                    args.mdrun_extra,
                )
            else:
                run_equilibration_restart_safe(args.gmx, work, log_dir, args.mdrun_extra)
    if not planned and not common_plan:
        raise SystemExit("No realization or common-minimization stage matched the request")
    print(json.dumps({
        "mode": "execute" if args.execute else "plan_only",
        "phase": args.phase,
        "release_stage": preflight_stage,
        "canary_target_ns_per_realization": 5.0 if args.phase == "canary" else None,
        "mdrun_extra_performance_only": args.mdrun_extra,
        "common_minimization": common_plan,
        "runs": planned,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
