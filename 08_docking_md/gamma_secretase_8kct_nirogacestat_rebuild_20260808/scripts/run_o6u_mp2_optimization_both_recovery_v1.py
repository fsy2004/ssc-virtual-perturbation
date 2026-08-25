#!/usr/bin/env python3
"""Run an audited O6U MP2 recovery with OptKing BOTH coordinates."""

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


RECOVERY_OPTIONS = {
    "dynamic_level": 1,
    "opt_coordinates": "both",
    "geom_maxiter": 200,
    "g_convergence": "gau_tight",
}
PSI4_OPTIONS = {
    "basis": "6-31G(d)", "reference": "rhf", "scf_type": "df", "mp2_type": "df",
    "freeze_core": True, "guess": "sad", "e_convergence": 1.0e-8,
    "d_convergence": 1.0e-8, "maxiter": 200, **RECOVERY_OPTIONS,
}


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--start-xyz", required=True, type=Path)
    parser.add_argument("--snapshot-report", required=True, type=Path)
    parser.add_argument("--plateau-report", required=True, type=Path)
    parser.add_argument("--interrupted-record", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.threads <= 24 or not 16 <= args.memory_gib <= 128:
        raise SystemExit("Resources must remain within 1-24 threads and 16-128 GiB")
    source_path, start_path = args.source_sdf.resolve(), args.start_xyz.resolve()
    snapshot_path, plateau_path = args.snapshot_report.resolve(), args.plateau_report.resolve()
    interrupted_path, output_dir = args.interrupted_record.resolve(), args.output_dir.resolve()
    for path in (source_path, start_path, snapshot_path, plateau_path, interrupted_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")

    snapshot, plateau, interrupted = map(load_json, (snapshot_path, plateau_path, interrupted_path))
    if snapshot.get("status") != "pass_technical_restart_candidate_not_converged" or snapshot.get("role") != "snapshot_only":
        raise SystemExit("Restart snapshot is not eligible")
    if plateau.get("status") != "pass_plateau_stop_authorized_no_convergence" or plateau.get("production_approved") is not False:
        raise SystemExit("Plateau termination authorization is missing")
    if plateau.get("snapshot_report", {}).get("sha256") != sha256(snapshot_path):
        raise SystemExit("Plateau report is not bound to the restart snapshot")
    if interrupted.get("status") != "stopped_plateau_sigterm" or interrupted.get("pid_confirmed_dead") is not True:
        raise SystemExit("Prior Cartesian recovery was not auditably stopped after plateau authorization")
    if interrupted.get("pid") != plateau.get("pid"):
        raise SystemExit("Interrupted PID differs from plateau authorization")
    raw_record = interrupted.get("closed_raw_output", {})
    raw_path = Path(str(raw_record.get("path", ""))).resolve()
    if not raw_path.is_file() or raw_record.get("sha256") != sha256(raw_path) or raw_record.get("size_bytes") != raw_path.stat().st_size:
        raise SystemExit("Interrupted-run raw output is not closed and hash-bound")
    restart = snapshot.get("restart_xyz", {})
    if Path(str(restart.get("path", ""))).resolve() != start_path or restart.get("sha256") != sha256(start_path):
        raise SystemExit("Start XYZ differs from the independent snapshot")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    if snapshot.get("identity") != identity:
        raise SystemExit("Snapshot identity differs from immutable source")
    frame = read_xyz_ensemble(start_path)
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(frame) != 1 or frame[0]["elements"] != expected_elements or len(expected_elements) != EXPECTED_ATOMS:
        raise SystemExit("Restart XYZ identity/order differs")
    coordinates = frame[0]["coordinates"]
    if not all(math.isfinite(float(value)) for value in coordinates.reshape(-1)):
        raise SystemExit("Restart XYZ contains non-finite coordinates")

    if args.preflight_only:
        import psi4
        psi4.set_options(PSI4_OPTIONS)
        print(json.dumps({"status": "pass_both_recovery_preflight_no_qm_executed", "qm_executed": False, "options": RECOVERY_OPTIONS}, sort_keys=True))
        return 0

    output_dir.mkdir(parents=True, exist_ok=False)
    scratch = args.scratch_root.resolve() / f"o6u_mp2_both_recovery_{os.getpid()}"
    scratch.mkdir(parents=True, exist_ok=False)
    os.environ["PSI_SCRATCH"] = str(scratch)
    import psi4

    raw_output = output_dir / "o6u_mp2_631gd_optimization_both_recovery_canary_v4.psi4.out"
    optimized_xyz = output_dir / "o6u_mp2_631gd_optimization_both_recovery_canary_v4.optimized.xyz"
    record_path = output_dir / "o6u_mp2_631gd_optimization_both_recovery_canary_v4.json"
    geometry_lines = ["0 1"] + [
        f"{element:<2s} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for element, (x, y, z) in zip(expected_elements, coordinates, strict=True)
    ] + ["units angstrom", "no_com", "no_reorient", "symmetry c1"]
    molecule = psi4.geometry("\n".join(geometry_lines))
    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gib} GiB")
    psi4.core.set_output_file(str(raw_output), False)
    psi4.set_options(PSI4_OPTIONS)
    record: dict[str, object] = {
        "schema_version": "1.0", "report_type": "o6u_mp2_631gd_optimization_both_recovery_canary",
        "role": "recovery_canary", "status": "running", "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(),
        "python_version": platform.python_version(), "psi4_version": psi4.__version__,
        "identity": identity, "charge_e": 0, "multiplicity": 1,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference", "optimizer_convergence": "gau_tight",
        "optimizer_strategy_version": "optking_both_dynamic_level1_recovery_v1",
        "optimizer_options": RECOVERY_OPTIONS, "threads": args.threads, "memory_gib": args.memory_gib,
        "source_sdf": {"path": str(source_path), "sha256": sha256(source_path)},
        "start_xyz": {"path": str(start_path), "sha256": sha256(start_path)},
        "snapshot_report": {"path": str(snapshot_path), "sha256": sha256(snapshot_path)},
        "plateau_report": {"path": str(plateau_path), "sha256": sha256(plateau_path)},
        "plateau_stop_closure": {"path": str(interrupted_path), "sha256": sha256(interrupted_path)},
        "raw_output": str(raw_output), "scratch": str(scratch), "pid": os.getpid(),
        "release_boundary": "Recovery canary only; force-field fitting, CHARMM-GUI, and MD remain blocked.",
    }
    write_record(record_path, record)
    started = time.time(); exit_code = 1
    try:
        energy, wavefunction = psi4.optimize("mp2", molecule=molecule, return_wfn=True)
        molecule.save_xyz_file(str(optimized_xyz), True)
        record.update({"status": "pass_optimization_both_recovery_canary", "final_energy_hartree": float(energy),
                       "wavefunction_energy_hartree": float(wavefunction.energy()),
                       "optimized_xyz": {"path": str(optimized_xyz), "sha256": sha256(optimized_xyz)}})
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
