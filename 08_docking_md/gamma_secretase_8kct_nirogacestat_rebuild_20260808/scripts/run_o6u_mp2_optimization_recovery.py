#!/usr/bin/env python3
"""Run the frozen fail-closed O6U MP2 Cartesian recovery canary."""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

from prepare_o6u_crest_input import EXPECTED_ATOMS, load_single_sdf, validate_identity
from run_o6u_mp2_optimization_canary import sha256, write_record
from validate_o6u_crest_ensemble import read_xyz_ensemble


POLICY_STATUS = "prospective_policy_frozen_not_triggered"
EXTRACTION_STATUS = "pass_technical_restart_candidate_not_converged"
RECOVERY_OPTIONS = {
    "dynamic_level": 4,
    "geom_maxiter": 200,
    "g_convergence": "gau_tight",
}
PSI4_OPTIONS = {
    "basis": "6-31G(d)",
    "reference": "rhf",
    "scf_type": "df",
    "mp2_type": "df",
    "freeze_core": True,
    "guess": "sad",
    "e_convergence": 1.0e-8,
    "d_convergence": 1.0e-8,
    "maxiter": 200,
    **RECOVERY_OPTIONS,
}


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def verify_artifact(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    expected_hash = record.get("sha256") or record.get("sha256_at_extraction")
    expected_size = record.get("size_bytes") or record.get("bytes_at_extraction")
    if not path.is_file() or sha256(path) != expected_hash:
        raise RuntimeError(f"Artifact hash verification failed: {label}")
    if expected_size is not None and path.stat().st_size != expected_size:
        raise RuntimeError(f"Artifact size verification failed: {label}")
    return path


def pid_exists(pid: object) -> bool:
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--start-xyz", required=True, type=Path)
    parser.add_argument("--failed-record", required=True, type=Path)
    parser.add_argument("--restart-extraction-report", required=True, type=Path)
    parser.add_argument("--recovery-policy", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.threads <= 24 or not 16 <= args.memory_gib <= 128:
        raise SystemExit("Recovery resources must remain within 1-24 threads and 16-128 GiB")

    source_path = args.source_sdf.resolve()
    start_path = args.start_xyz.resolve()
    failed_path = args.failed_record.resolve()
    extraction_path = args.restart_extraction_report.resolve()
    policy_path = args.recovery_policy.resolve()
    output_dir = args.output_dir.resolve()
    for path in (source_path, start_path, failed_path, extraction_path, policy_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse recovery output directory: {output_dir}")

    failed = load_json(failed_path)
    extraction = load_json(extraction_path)
    policy = load_json(policy_path)
    if policy.get("schema_version") != "2.0" or policy.get("status") != POLICY_STATUS:
        raise SystemExit("Recovery policy is not the authoritative frozen v2 policy")
    profile = policy.get("technical_recovery_profile")
    invariants = policy.get("frozen_scientific_invariants")
    if not isinstance(profile, dict) or not isinstance(invariants, dict):
        raise SystemExit("Recovery policy is missing profile or scientific invariants")
    if (
        profile.get("dynamic_level") != 4
        or profile.get("dynamic_level_max_explicitly_set") is not False
        or profile.get("geom_maxiter") != 200
        or profile.get("model_chemistry_weakened") is not False
        or profile.get("convergence_weakened") is not False
        or invariants.get("method") != "DF-MP2/6-31G(d), frozen core, RHF reference"
        or invariants.get("geometry_convergence") != "gau_tight"
        or invariants.get("charge_e") != 0
        or invariants.get("multiplicity") != 1
    ):
        raise SystemExit("Recovery policy differs from the exact frozen profile")
    if (
        failed.get("status") != "fail"
        or failed.get("error_type") != "OptimizationConvergenceError"
        or failed.get("method") != invariants.get("method")
        or failed.get("optimizer_convergence") != "gau_tight"
        or failed.get("charge_e") != 0
        or failed.get("multiplicity") != 1
    ):
        raise SystemExit("Failed record does not pass the exact max-iteration recovery trigger")
    failed_source = failed.get("source_sdf")
    if (
        not isinstance(failed_source, dict)
        or Path(str(failed_source.get("path", ""))).resolve() != source_path
        or failed_source.get("sha256") != sha256(source_path)
    ):
        raise SystemExit("Failed run is not bound to the immutable source SDF")
    if pid_exists(failed.get("pid")):
        raise SystemExit("Failed-record PID still exists; refusing recovery")
    failed_raw = Path(str(failed.get("raw_output", ""))).resolve()
    if not failed_raw.is_file() or failed.get("raw_output_sha256") != sha256(failed_raw):
        raise SystemExit("Closed failed-run raw output is missing or hash-inconsistent")
    error_text = str(failed.get("error", "")).lower()
    if "maximum number of steps" not in error_text and "geom_maxiter" not in error_text:
        raise SystemExit("Failure is not the prespecified geometry-iteration exhaustion")

    if (
        extraction.get("status") != EXTRACTION_STATUS
        or extraction.get("role") != "failed_run_restart_candidate"
        or extraction.get("production_approved") is not False
    ):
        raise SystemExit("Restart extraction report does not pass its exact failed-run gate")
    extraction_record = extraction.get("record")
    extraction_raw = extraction.get("raw_output")
    if not isinstance(extraction_record, dict) or not isinstance(extraction_raw, dict):
        raise SystemExit("Restart extraction bindings are missing")
    if extraction_record.get("sha256") != sha256(failed_path) or extraction_raw.get("sha256_at_extraction") != sha256(failed_raw):
        raise SystemExit("Restart extraction is not bound to the closed failed run")
    declared_restart = verify_artifact(extraction.get("restart_xyz"), "extraction.restart_xyz")
    if declared_restart != start_path:
        raise SystemExit("Recovery start XYZ differs from the independently extracted restart candidate")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    if extraction.get("identity") != identity:
        raise SystemExit("Restart extraction identity differs from the immutable source")
    frames = read_xyz_ensemble(start_path)
    if len(frames) != 1:
        raise SystemExit("Recovery start must contain exactly one XYZ frame")
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    frame = frames[0]
    if frame["elements"] != expected_elements or len(expected_elements) != EXPECTED_ATOMS:
        raise SystemExit("Recovery XYZ element sequence differs from immutable O6U SDF")
    coordinates = frame["coordinates"]
    if not all(math.isfinite(float(value)) for value in coordinates.reshape(-1)):
        raise SystemExit("Recovery XYZ contains non-finite coordinates")

    if args.preflight_only:
        import psi4

        psi4.set_options(PSI4_OPTIONS)
        print(json.dumps({
            "status": "pass_recovery_runner_preflight_no_qm_executed",
            "production_approved": False,
            "qm_executed": False,
            "recovery_options": RECOVERY_OPTIONS,
        }, sort_keys=True))
        return 0

    output_dir.mkdir(parents=True, exist_ok=False)
    scratch = args.scratch_root.resolve() / f"o6u_mp2_recovery_{os.getpid()}"
    scratch.mkdir(parents=True, exist_ok=False)
    os.environ["PSI_SCRATCH"] = str(scratch)
    import psi4

    raw_output = output_dir / "o6u_mp2_631gd_optimization_recovery_canary_v3.psi4.out"
    optimized_xyz = output_dir / "o6u_mp2_631gd_optimization_recovery_canary_v3.optimized.xyz"
    record_path = output_dir / "o6u_mp2_631gd_optimization_recovery_canary_v3.json"
    geometry_lines = ["0 1"]
    geometry_lines.extend(
        f"{element:<2s} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for element, (x, y, z) in zip(expected_elements, coordinates, strict=True)
    )
    geometry_lines.extend(["units angstrom", "no_com", "no_reorient", "symmetry c1"])
    molecule = psi4.geometry("\n".join(geometry_lines))
    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gib} GiB")
    psi4.core.set_output_file(str(raw_output), False)
    psi4.set_options(PSI4_OPTIONS)
    record: dict[str, object] = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_631gd_optimization_recovery_canary",
        "role": "recovery_canary",
        "status": "running",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(), "python_version": platform.python_version(),
        "psi4_version": psi4.__version__, "identity": identity, "charge_e": 0, "multiplicity": 1,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "optimizer_convergence": "gau_tight",
        "optimizer_strategy_version": "cartesian_dynamic_level4_recovery_v2",
        "optimizer_options": RECOVERY_OPTIONS,
        "threads": args.threads, "memory_gib": args.memory_gib,
        "source_sdf": {"path": str(source_path), "sha256": sha256(source_path)},
        "start_xyz": {"path": str(start_path), "sha256": sha256(start_path)},
        "failed_v2_record": {"path": str(failed_path), "sha256": sha256(failed_path)},
        "restart_extraction_report": {"path": str(extraction_path), "sha256": sha256(extraction_path)},
        "recovery_policy": {"path": str(policy_path), "sha256": sha256(policy_path)},
        "raw_output": str(raw_output), "scratch": str(scratch), "pid": os.getpid(),
        "release_boundary": "Recovery canary only; no representative, force-field parameter, or MD is approved.",
    }
    write_record(record_path, record)
    started = time.time()
    exit_code = 1
    try:
        energy, wavefunction = psi4.optimize("mp2", molecule=molecule, return_wfn=True)
        molecule.save_xyz_file(str(optimized_xyz), True)
        record.update({
            "status": "pass_optimization_recovery_canary",
            "final_energy_hartree": float(energy),
            "wavefunction_energy_hartree": float(wavefunction.energy()),
            "optimized_xyz": {"path": str(optimized_xyz), "sha256": sha256(optimized_xyz)},
        })
        exit_code = 0
    except BaseException as exc:
        record.update({"status": "fail", "error_type": type(exc).__name__, "error": str(exc)})
    finally:
        record["elapsed_seconds"] = time.time() - started
        psi4.core.flush_outfile()
        if raw_output.is_file():
            record["raw_output_sha256"] = sha256(raw_output)
            record["raw_output_bytes"] = raw_output.stat().st_size
        write_record(record_path, record)
        print(json.dumps({"status": record["status"], "elapsed_seconds": record["elapsed_seconds"], "record_sha256": sha256(record_path)}, sort_keys=True))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
