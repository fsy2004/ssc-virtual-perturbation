#!/usr/bin/env python3
"""Run an audited O6U MP2 recovery canary with geomeTRIC/TRIC."""

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


EXPECTED_GEOMETRIC_VERSION = "1.1.1"
GEOMETRIC_OPTIONS = {"coordsys": "tric"}
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
    "geom_maxiter": 200,
    "g_convergence": "gau_tight",
}


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def verify_inputs(args: argparse.Namespace) -> tuple[dict[str, object], list[str], object]:
    source_path = args.source_sdf.resolve()
    start_path = args.start_xyz.resolve()
    snapshot_path = args.snapshot_report.resolve()
    stall_path = args.stall_report.resolve()
    closure_path = args.stop_closure.resolve()
    for path in (source_path, start_path, snapshot_path, stall_path, closure_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")

    snapshot, stall, closure = map(load_json, (snapshot_path, stall_path, closure_path))
    if snapshot.get("status") != "pass_technical_restart_candidate_not_converged" or snapshot.get("role") != "snapshot_only":
        raise SystemExit("Restart snapshot is not eligible")
    if stall.get("status") != "pass_displacement_stall_stop_authorized_no_convergence" or stall.get("production_approved") is not False:
        raise SystemExit("Displacement-stall authorization is missing")
    if closure.get("status") != "stopped_displacement_stall_sigterm" or closure.get("pid_confirmed_dead") is not True:
        raise SystemExit("OptKing BOTH canary was not auditably stopped")
    if stall.get("pid") != closure.get("pid"):
        raise SystemExit("Stop closure PID differs from stall authorization")
    if stall.get("snapshot_report", {}).get("sha256") != sha256(snapshot_path):
        raise SystemExit("Stall report is not bound to the restart snapshot")
    if closure.get("termination_authorization", {}).get("sha256") != sha256(stall_path):
        raise SystemExit("Stop closure is not bound to the stall authorization")
    closed_raw = closure.get("closed_raw_output", {})
    closed_raw_path = Path(str(closed_raw.get("path", ""))).resolve()
    if (
        not closed_raw_path.is_file()
        or closed_raw.get("sha256") != sha256(closed_raw_path)
        or closed_raw.get("size_bytes") != closed_raw_path.stat().st_size
    ):
        raise SystemExit("Stopped OptKing raw output is not closed and hash-bound")
    restart = snapshot.get("restart_xyz", {})
    if Path(str(restart.get("path", ""))).resolve() != start_path or restart.get("sha256") != sha256(start_path):
        raise SystemExit("Start XYZ differs from the independent snapshot")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    if snapshot.get("identity") != identity:
        raise SystemExit("Snapshot identity differs from immutable source")
    frames = read_xyz_ensemble(start_path)
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(frames) != 1 or frames[0]["elements"] != expected_elements or len(expected_elements) != EXPECTED_ATOMS:
        raise SystemExit("Restart XYZ identity/order differs")
    coordinates = frames[0]["coordinates"]
    if not all(math.isfinite(float(value)) for value in coordinates.reshape(-1)):
        raise SystemExit("Restart XYZ contains non-finite coordinates")
    return identity, expected_elements, coordinates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--start-xyz", required=True, type=Path)
    parser.add_argument("--snapshot-report", required=True, type=Path)
    parser.add_argument("--stall-report", required=True, type=Path)
    parser.add_argument("--stop-closure", required=True, type=Path)
    parser.add_argument("--geometric-sdist", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.threads <= 24 or not 16 <= args.memory_gib <= 128:
        raise SystemExit("Resources must remain within 1-24 threads and 16-128 GiB")
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    sdist_path = args.geometric_sdist.resolve()
    if not sdist_path.is_file() or sdist_path.stat().st_size == 0:
        raise SystemExit("Pinned geomeTRIC source distribution is missing")
    identity, expected_elements, coordinates = verify_inputs(args)

    import geometric
    import psi4

    if geometric.__version__ != EXPECTED_GEOMETRIC_VERSION:
        raise SystemExit(f"Unexpected geomeTRIC version: {geometric.__version__}")
    psi4.set_options(PSI4_OPTIONS)
    if args.preflight_only:
        print(json.dumps({
            "status": "pass_geometric_recovery_preflight_no_qm_executed",
            "qm_executed": False,
            "psi4_version": psi4.__version__,
            "geometric_version": geometric.__version__,
            "engine": "geometric",
            "optimizer_keywords": GEOMETRIC_OPTIONS,
            "model_chemistry": "DF-MP2/6-31G(d), frozen core, RHF reference",
            "convergence": "gau_tight",
            "geometric_sdist_sha256": sha256(sdist_path),
        }, sort_keys=True))
        return 0

    output_dir.mkdir(parents=True, exist_ok=False)
    scratch = args.scratch_root.resolve() / f"o6u_mp2_geometric_recovery_{os.getpid()}"
    scratch.mkdir(parents=True, exist_ok=False)
    os.environ["PSI_SCRATCH"] = str(scratch)
    raw_output = output_dir / "o6u_mp2_631gd_optimization_geometric_recovery_canary_v1.psi4.out"
    optimized_xyz = output_dir / "o6u_mp2_631gd_optimization_geometric_recovery_canary_v1.optimized.xyz"
    record_path = output_dir / "o6u_mp2_631gd_optimization_geometric_recovery_canary_v1.json"
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
        "schema_version": "1.0",
        "report_type": "o6u_mp2_631gd_optimization_geometric_recovery_canary",
        "role": "recovery_canary",
        "status": "running",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "psi4_version": psi4.__version__,
        "geometric_version": geometric.__version__,
        "identity": identity,
        "charge_e": 0,
        "multiplicity": 1,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "optimizer_convergence": "gau_tight",
        "optimizer_strategy_version": "geometric_tric_recovery_v1",
        "optimizer_engine": "geometric",
        "optimizer_keywords": GEOMETRIC_OPTIONS,
        "geom_maxiter": 200,
        "threads": args.threads,
        "memory_gib": args.memory_gib,
        "source_sdf": {"path": str(args.source_sdf.resolve()), "sha256": sha256(args.source_sdf.resolve())},
        "start_xyz": {"path": str(args.start_xyz.resolve()), "sha256": sha256(args.start_xyz.resolve())},
        "snapshot_report": {"path": str(args.snapshot_report.resolve()), "sha256": sha256(args.snapshot_report.resolve())},
        "stall_report": {"path": str(args.stall_report.resolve()), "sha256": sha256(args.stall_report.resolve())},
        "stop_closure": {"path": str(args.stop_closure.resolve()), "sha256": sha256(args.stop_closure.resolve())},
        "geometric_sdist": {"path": str(sdist_path), "sha256": sha256(sdist_path)},
        "references": [
            "https://psi4.github.io/psi4docs/master/optking.html",
            "https://psi4.github.io/psi4docs/master/opt.html",
            "https://doi.org/10.1063/1.4952956",
        ],
        "raw_output": str(raw_output),
        "scratch": str(scratch),
        "pid": os.getpid(),
        "release_boundary": "Recovery canary only; parameter fitting, CHARMM-GUI, and MD remain blocked.",
    }
    write_record(record_path, record)
    started = time.time()
    exit_code = 1
    try:
        energy, wavefunction = psi4.optimize(
            "mp2",
            molecule=molecule,
            return_wfn=True,
            engine="geometric",
            optimizer_keywords=GEOMETRIC_OPTIONS,
        )
        molecule.save_xyz_file(str(optimized_xyz), True)
        record.update({
            "status": "pass_optimization_geometric_recovery_canary",
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
