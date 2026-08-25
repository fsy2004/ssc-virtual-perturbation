#!/usr/bin/env python3
"""Run one hash-bound relaxed O6U torsion point with geomeTRIC."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

from run_o6u_mp2_torsion_canary_v1 import (
    circular_delta,
    dihedral,
    parse_mol2_graph,
    read_xyz,
    rotate_component,
    sha256,
    side_component,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--frame342", required=True, type=Path)
    parser.add_argument("--mol2", required=True, type=Path)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--signed-step-index", required=True, type=int)
    parser.add_argument("--start-xyz", required=True, type=Path)
    parser.add_argument("--start-xyz-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=44)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    if args.signed_step_index == 0:
        raise SystemExit("signed-step-index must be non-zero")

    authorization_path = args.authorization.resolve()
    authorization = json.loads(authorization_path.read_text())
    if authorization.get("status") != "pass_bonded_torsion_scope_authorized_canary_only":
        raise SystemExit("Invalid bonded/torsion authorization")

    frame342 = args.frame342.resolve()
    mol2 = args.mol2.resolve()
    start_xyz = args.start_xyz.resolve()
    if sha256(frame342) != authorization["inputs"]["frame342"]["sha256"]:
        raise SystemExit("Frame342 hash does not match authorization")
    if sha256(mol2) != authorization["inputs"]["mol2"]["sha256"]:
        raise SystemExit("MOL2 hash does not match authorization")
    if sha256(start_xyz) != args.start_xyz_sha256:
        raise SystemExit("Start XYZ hash does not match the bound previous point")

    rotor = next(
        (
            item
            for item in authorization["rotatable_heavy_atom_torsions"]
            if item["rotor_id"] == args.rotor_id
        ),
        None,
    )
    if rotor is None:
        raise SystemExit(f"Unauthorized rotor: {args.rotor_id}")

    atom_ids = [int(value) for value in rotor["one_based_ordinals"]]
    frame_elements, frame_coordinates, _ = read_xyz(frame342)
    start_elements, start_coordinates, _ = read_xyz(start_xyz)
    if start_elements != frame_elements:
        raise SystemExit("Start geometry atom identities/order differ from frame342")

    spacing = 15.0
    frame_angle = dihedral(frame_coordinates, atom_ids)
    start_angle = dihedral(start_coordinates, atom_ids)
    target_angle = frame_angle + args.signed_step_index * spacing
    previous_index = args.signed_step_index - (1 if args.signed_step_index > 0 else -1)
    expected_start_angle = frame_angle + previous_index * spacing
    start_error = abs(circular_delta(start_angle, expected_start_angle))
    if start_error > 0.5:
        raise SystemExit(
            f"Start geometry is not the previous grid point: error={start_error:.6f} deg"
        )

    graph = parse_mol2_graph(mol2)
    j_atom, k_atom = atom_ids[1], atom_ids[2]
    component = side_component(graph, j_atom, k_atom)
    direction = 1.0 if args.signed_step_index > 0 else -1.0
    choices = [
        rotate_component(start_coordinates, component, j_atom, k_atom, sign * spacing)
        for sign in (direction, -direction)
    ]
    coordinates = min(
        choices,
        key=lambda candidate: abs(
            circular_delta(dihedral(candidate, atom_ids), target_angle)
        ),
    )
    rotation_error = abs(circular_delta(dihedral(coordinates, atom_ids), target_angle))
    if rotation_error > 1e-5:
        raise SystemExit(f"Target rotation failed: error={rotation_error:.8f} deg")

    optimizer_keywords = {
        "coordsys": "tric",
        "constraints": {
            "set": [
                {
                    "type": "dihedral",
                    "indices": [value - 1 for value in atom_ids],
                    "value": target_angle,
                }
            ]
        },
        "conmethod": 1,
        "maxiter": 200,
        "convergence_set": "GAU_TIGHT",
    }
    runner_path = Path(__file__).resolve()
    preflight = {
        "status": "pass_geometric_torsion_scan_point_preflight",
        "authorization_sha256": sha256(authorization_path),
        "runner_sha256": sha256(runner_path),
        "rotor_id": rotor["rotor_id"],
        "one_based_ordinals": atom_ids,
        "signed_step_index": args.signed_step_index,
        "grid_spacing_deg": spacing,
        "frame342_angle_deg": frame_angle,
        "start_angle_deg": start_angle,
        "start_angle_error_deg": start_error,
        "target_angle_deg": target_angle,
        "frame342_sha256": sha256(frame342),
        "mol2_sha256": sha256(mol2),
        "start_xyz_sha256": sha256(start_xyz),
        "optimizer": "geomeTRIC 1.1.1 TRIC exact-set constraint conmethod=1",
        "optimizer_keywords": optimizer_keywords,
        "method": "DF-MP2/6-31G(d), frozen-core RHF, GAU_TIGHT",
        "charge": 0,
        "multiplicity": 1,
        "threads": args.threads,
        "memory_gib": args.memory_gib,
    }
    if args.preflight_only:
        print(json.dumps(preflight, sort_keys=True))
        return 0

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    pid_path = output_dir.parent / f"{output_dir.name}.pid"
    pid_path.write_text(f"{os.getpid()}\n")

    scratch = (
        args.scratch_root.resolve()
        / f"o6u_{rotor['rotor_id'].lower()}_step_{args.signed_step_index:+03d}_{os.getpid()}"
    )
    scratch.mkdir(parents=True, exist_ok=False)
    os.environ["PSI_SCRATCH"] = str(scratch)
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = str(args.threads)

    import geometric
    import psi4

    if geometric.__version__ != "1.1.1":
        raise SystemExit(f"Unexpected geomeTRIC version: {geometric.__version__}")
    if psi4.__version__ != "1.9.1":
        raise SystemExit(f"Unexpected Psi4 version: {psi4.__version__}")

    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gib} GiB")
    raw_output = output_dir / "o6u_mp2_torsion_scan_point.psi4.out"
    psi4.core.set_output_file(str(raw_output), False)
    geometry = ["0 1"] + [
        f"{element:<2s} {coordinate[0]: .10f} {coordinate[1]: .10f} {coordinate[2]: .10f}"
        for element, coordinate in zip(frame_elements, coordinates, strict=True)
    ] + ["units angstrom", "no_com", "no_reorient", "symmetry c1"]
    molecule = psi4.geometry("\n".join(geometry))
    psi4.set_options(
        {
            "basis": "6-31G(d)",
            "reference": "rhf",
            "scf_type": "df",
            "mp2_type": "df",
            "freeze_core": True,
            "guess": "sad",
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "maxiter": 200,
            "geom_maxiter": 200,
            "g_convergence": "gau_tight",
        }
    )

    report_path = output_dir / "O6U_MP2_TORSION_SCAN_POINT.json"
    report = {
        **preflight,
        "schema_version": "1.0",
        "report_type": "o6u_mp2_torsion_scan_point",
        "status": "running",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "psi4_version": psi4.__version__,
        "geometric_version": geometric.__version__,
        "pid": os.getpid(),
        "pid_file": str(pid_path),
        "raw_output": str(raw_output),
        "scratch": str(scratch),
        "production_approved": False,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    started = time.time()
    exit_code = 1
    try:
        energy = psi4.optimize(
            "mp2",
            molecule=molecule,
            engine="geometric",
            optimizer_keywords=optimizer_keywords,
        )
        optimized_xyz = output_dir / "o6u_mp2_torsion_scan_point.optimized.xyz"
        molecule.save_xyz_file(str(optimized_xyz), True)
        final_elements, final_coordinates, _ = read_xyz(optimized_xyz)
        if final_elements != frame_elements:
            raise RuntimeError("Optimized geometry atom identities/order changed")
        final_angle = dihedral(final_coordinates, atom_ids)
        target_error = abs(circular_delta(final_angle, target_angle))
        if target_error > 0.5:
            raise RuntimeError(f"Constraint error {target_error:.6f} deg")
        report.update(
            {
                "status": "pass_relaxed_mp2_torsion_scan_point",
                "final_energy_hartree": float(energy),
                "final_dihedral_deg": final_angle,
                "target_error_deg": target_error,
                "optimized_xyz": {
                    "path": str(optimized_xyz),
                    "sha256": sha256(optimized_xyz),
                },
            }
        )
        exit_code = 0
    except BaseException as exc:
        report.update(
            {
                "status": "fail",
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
        )
    finally:
        psi4.core.flush_outfile()
        report["elapsed_seconds"] = time.time() - started
        if raw_output.is_file():
            report["raw_output_sha256"] = sha256(raw_output)
            report["raw_output_bytes"] = raw_output.stat().st_size
        if exit_code == 0:
            report["scratch_cleanup"] = "removed_after_success"
            shutil.rmtree(scratch)
        else:
            report["scratch_cleanup"] = "preserved_after_failure"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "status": report["status"],
                    "report": str(report_path),
                    "sha256": sha256(report_path),
                },
                sort_keys=True,
            )
        )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
