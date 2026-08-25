#!/usr/bin/env python3
"""Generate an auditable Psi4 sow/reap MP2 Hessian task set without running QM."""

from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

from prepare_o6u_crest_input import EXPECTED_ATOMS, load_single_sdf, validate_identity
from run_o6u_mp2_optimization_canary import sha256
from validate_o6u_crest_ensemble import read_xyz_ensemble


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
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--optimized-xyz", required=True, type=Path)
    parser.add_argument("--optimization-record", required=True, type=Path)
    parser.add_argument("--independent-validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads-per-job", type=int, default=3)
    parser.add_argument("--memory-gib-per-job", type=int, default=8)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    paths = [args.source_sdf.resolve(), args.optimized_xyz.resolve(), args.optimization_record.resolve(), args.independent_validation.resolve()]
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty input: {path}")
    if not 1 <= args.threads_per_job <= 4 or not 4 <= args.memory_gib_per_job <= 12:
        raise SystemExit("Per-job resources exceed the frozen Hessian bounds")

    record = json.loads(paths[2].read_text(encoding="utf-8"))
    validation = json.loads(paths[3].read_text(encoding="utf-8"))
    if record.get("status") != "pass_optimization_geometric_recovery_canary":
        raise SystemExit("Optimization record is not a completed geomeTRIC canary")
    if validation.get("status") != "pass_geometric_recovery_independently_reconstructed_pending_minimum_character":
        raise SystemExit("Independent validation is not the pending-minimum pass")
    if validation.get("record", {}).get("sha256") != sha256(paths[2]):
        raise SystemExit("Independent validation is not bound to the optimization record")
    if record.get("optimized_xyz", {}).get("sha256") != sha256(paths[1]):
        raise SystemExit("Optimized XYZ differs from the closed record")

    source = load_single_sdf(paths[0])
    identity = validate_identity(source)
    frames = read_xyz_ensemble(paths[1])
    elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(frames) != 1 or frames[0]["elements"] != elements or len(elements) != EXPECTED_ATOMS:
        raise SystemExit("Optimized XYZ identity/order differs")
    coordinates = frames[0]["coordinates"]
    if not all(math.isfinite(float(value)) for value in coordinates.reshape(-1)):
        raise SystemExit("Optimized coordinates are non-finite")

    output_dir.mkdir(parents=True, exist_ok=False)
    generation_output = output_dir / "sow_generation.psi4.out"
    old_cwd = Path.cwd()
    os.chdir(output_dir)
    try:
        import psi4

        psi4.set_num_threads(args.threads_per_job)
        psi4.set_memory(f"{args.memory_gib_per_job} GiB")
        psi4.core.set_output_file(str(generation_output), False)
        psi4.set_options(PSI4_OPTIONS)
        geometry_lines = ["0 1"] + [
            f"{element:<2s} {float(x): .12f} {float(y): .12f} {float(z): .12f}"
            for element, (x, y, z) in zip(elements, coordinates, strict=True)
        ] + ["units angstrom", "no_com", "no_reorient", "symmetry c1"]
        molecule = psi4.geometry("\n".join(geometry_lines))
        sow_exception: dict[str, object] | None = None
        try:
            psi4.frequency("mp2", molecule=molecule, dertype="gradient", mode="sow")
        except SystemExit as exc:
            sow_exception = {"type": "SystemExit", "code": exc.code}
        finally:
            psi4.core.flush_outfile()
    finally:
        os.chdir(old_cwd)

    input_files = sorted(output_dir.glob("*.in"))
    master_files = [path for path in input_files if "master" in path.name.lower()]
    worker_files = [path for path in input_files if path not in master_files]
    if not master_files or not worker_files:
        raise SystemExit(f"Psi4 sow did not generate master and worker inputs: {[p.name for p in input_files]}")
    if len(worker_files) < EXPECTED_ATOMS * 3:
        raise SystemExit(f"Unexpectedly small displacement task set: {len(worker_files)}")
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_hessian_sow_generation",
        "status": "pass_mp2_hessian_sow_generated_no_qm_executed",
        "qm_executed": False,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "derivative_route": "finite_difference_of_analytic_gradients",
        "threads_per_job": args.threads_per_job,
        "memory_gib_per_job": args.memory_gib_per_job,
        "optimization_record": {"path": str(paths[2]), "sha256": sha256(paths[2])},
        "independent_validation": {"path": str(paths[3]), "sha256": sha256(paths[3])},
        "optimized_xyz": {"path": str(paths[1]), "sha256": sha256(paths[1])},
        "sow_generation_output": {"path": str(generation_output), "sha256": sha256(generation_output), "bytes": generation_output.stat().st_size},
        "sow_exception": sow_exception,
        "worker_input_count": len(worker_files),
        "master_input_count": len(master_files),
        "worker_inputs": [{"name": path.name, "sha256": sha256(path), "bytes": path.stat().st_size} for path in worker_files],
        "master_inputs": [{"name": path.name, "sha256": sha256(path), "bytes": path.stat().st_size} for path in master_files],
        "references": [
            "https://psi4.github.io/psi4docs/master/freq.html",
            "https://psi4.github.io/psi4docs/master/sowreap.html",
        ],
        "release_boundary": "Input generation only; minimum character, parameter fitting, CHARMM-GUI, and MD remain blocked.",
    }
    report_path = output_dir / "O6U_MP2_HESSIAN_SOW_GENERATION.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report": str(report_path), "sha256": sha256(report_path), "worker_input_count": len(worker_files)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
