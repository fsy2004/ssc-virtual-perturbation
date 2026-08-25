#!/usr/bin/env python3
"""Prepare a restartable, parallel QCSchema MP2 finite-difference Hessian plan."""

from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

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
    "points": 3,
    "disp_size": 0.005,
}


def last_total_gradient(path: Path, natom: int):
    text = path.read_text(encoding="utf-8", errors="replace")
    starts = [m.end() for m in re.finditer(r"-Total Gradient:\s*\n", text)]
    if not starts:
        raise SystemExit("No total-gradient block in the closed optimization output")
    rows = []
    for line in text[starts[-1] :].splitlines():
        m = re.match(
            r"\s*\d+\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
            line,
        )
        if m:
            rows.append([float(m.group(i)) for i in (1, 2, 3)])
            if len(rows) == natom:
                break
        elif rows:
            break
    if len(rows) != natom or not all(math.isfinite(x) for row in rows for x in row):
        raise SystemExit(f"Could not reconstruct the final {natom}-atom gradient")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-sdf", required=True, type=Path)
    ap.add_argument("--optimized-xyz", required=True, type=Path)
    ap.add_argument("--optimization-record", required=True, type=Path)
    ap.add_argument("--independent-validation", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--threads-per-job", type=int, default=3)
    ap.add_argument("--memory-gib-per-job", type=int, default=8)
    args = ap.parse_args()

    paths = [
        args.source_sdf.resolve(),
        args.optimized_xyz.resolve(),
        args.optimization_record.resolve(),
        args.independent_validation.resolve(),
    ]
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty input: {path}")
    out = args.output_dir.resolve()
    if out.exists():
        raise SystemExit(f"Refusing to reuse output directory: {out}")
    if not 1 <= args.threads_per_job <= 4 or not 4 <= args.memory_gib_per_job <= 12:
        raise SystemExit("Per-job resources exceed frozen bounds")

    record = json.loads(paths[2].read_text(encoding="utf-8"))
    validation = json.loads(paths[3].read_text(encoding="utf-8"))
    if record.get("status") != "pass_optimization_geometric_recovery_canary":
        raise SystemExit("Optimization record is not the closed geomeTRIC pass")
    if validation.get("status") != "pass_geometric_recovery_independently_reconstructed_pending_minimum_character":
        raise SystemExit("Independent validation is not the pending-minimum pass")
    if validation.get("record", {}).get("sha256") != sha256(paths[2]):
        raise SystemExit("Independent validation does not bind the optimization record")
    if record.get("optimized_xyz", {}).get("sha256") != sha256(paths[1]):
        raise SystemExit("Optimized XYZ differs from the closed record")
    raw = Path(record["raw_output"]).resolve()
    if not raw.is_file() or sha256(raw) != record.get("raw_output_sha256"):
        raise SystemExit("Closed optimization raw output is missing or hash-mismatched")

    source = load_single_sdf(paths[0])
    identity = validate_identity(source)
    frames = read_xyz_ensemble(paths[1])
    elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(frames) != 1 or frames[0]["elements"] != elements or len(elements) != EXPECTED_ATOMS:
        raise SystemExit("Optimized XYZ atom identity/order differs")
    coordinates = frames[0]["coordinates"]
    if not all(math.isfinite(float(x)) for x in coordinates.reshape(-1)):
        raise SystemExit("Optimized coordinates are non-finite")
    gradient = last_total_gradient(raw, EXPECTED_ATOMS)

    out.mkdir(parents=True, exist_ok=False)
    task_dir = out / "tasks"
    result_dir = out / "results"
    scratch_dir = out / "scratch"
    task_dir.mkdir()
    result_dir.mkdir()
    scratch_dir.mkdir()

    import psi4

    psi4.core.be_quiet()
    psi4.set_num_threads(args.threads_per_job)
    psi4.set_memory(f"{args.memory_gib_per_job} GiB")
    psi4.set_options(PSI4_OPTIONS)
    geometry_lines = ["0 1"] + [
        f"{el:<2s} {float(x): .12f} {float(y): .12f} {float(z): .12f}"
        for el, (x, y, z) in zip(elements, coordinates, strict=True)
    ] + ["units angstrom", "no_com", "no_reorient", "symmetry c1"]
    molecule = psi4.geometry("\n".join(geometry_lines))
    ref_gradient = psi4.core.Matrix.from_array(np.asarray(gradient, dtype=float))
    plan = psi4.hessian(
        "mp2",
        molecule=molecule,
        dertype="gradient",
        return_plan=True,
        ref_gradient=ref_gradient,
    )
    if type(plan).__name__ != "FiniteDifferenceComputer":
        raise SystemExit(f"Unexpected plan type: {type(plan)!r}")
    expected = 2 * (3 * EXPECTED_ATOMS - 6) + 1
    if len(plan.task_list) != expected:
        raise SystemExit(f"Expected {expected} projected-gradient tasks, got {len(plan.task_list)}")

    tasks = []
    for index, (label, task) in enumerate(plan.task_list.items()):
        safe = "reference" if label == "reference" else f"disp_{index:03d}"
        task_path = task_dir / f"{safe}.json"
        task_path.write_text(task.plan().json() + "\n", encoding="utf-8", newline="\n")
        tasks.append({
            "index": index,
            "label": label,
            "input": str(task_path),
            "input_sha256": sha256(task_path),
            "input_bytes": task_path.stat().st_size,
            "output": str(result_dir / f"{safe}.result.json"),
            "scratch": str(scratch_dir / safe),
        })

    manifest = {
        "schema_version": "2.0",
        "report_type": "o6u_mp2_hessian_qcschema_plan",
        "status": "pass_parallel_qcschema_plan_generated_no_qm_executed",
        "qm_executed": False,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "derivative_route": "3-point finite difference of analytic gradients",
        "displacement_bohr": 0.005,
        "translations_projected": True,
        "rotations_projected": True,
        "reference_gradient_rms_hartree_per_bohr": math.sqrt(sum(x*x for row in gradient for x in row) / (3*EXPECTED_ATOMS)),
        "threads_per_job": args.threads_per_job,
        "memory_gib_per_job": args.memory_gib_per_job,
        "task_count": len(tasks),
        "optimization_record": {"path": str(paths[2]), "sha256": sha256(paths[2])},
        "independent_validation": {"path": str(paths[3]), "sha256": sha256(paths[3])},
        "optimized_xyz": {"path": str(paths[1]), "sha256": sha256(paths[1])},
        "optimization_raw": {"path": str(raw), "sha256": sha256(raw)},
        "tasks": tasks,
        "references": [
            {"title": "Psi4 harmonic vibrational analysis documentation", "url": "https://psicode.org/psi4manual/master/freq.html", "decision": "Use finite differences of analytic gradients because Psi4 1.9.1 has no analytic MP2 Hessian."},
            {"title": "Psi4 FINDIF documentation", "url": "https://psicode.org/psi4manual/master/autodir_options_c/module__findif.html", "decision": "Use the documented 3-point stencil and 0.005 bohr displacement."},
            {"title": "Wilson, Decius and Cross, Molecular Vibrations (1955)", "url": "https://archive.org/details/molecularvibrati0000wils", "decision": "Project rigid translations and rotations for a nonlinear stationary molecule."},
        ],
        "compound_context": "Neutral-singlet, identity- and chirality-validated nirogacestat geometry; the check validates only local minimum character for parameterization.",
        "disease_context": "This technical gate does not provide SSc efficacy evidence; SSc interpretation remains restricted to disease-specific HES1-Notch evidence.",
        "release_boundary": "No parameter fitting, CHARMM-GUI, or MD release until every result and the assembled minimum-character audit pass.",
    }
    manifest_path = out / "O6U_MP2_HESSIAN_QCSCHEMA_PLAN.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": manifest["status"], "manifest": str(manifest_path), "sha256": sha256(manifest_path), "task_count": len(tasks)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
