#!/usr/bin/env python3
"""Run a non-target HF/6-31G(d) O6U single-point QM environment canary.

This verifies the frozen molecule, Psi4 runtime, integral engine, SCF stack,
threading, memory, and scratch configuration. It is not parameter target data
and must never be counted as FFParam validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import time
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=80)
    args = parser.parse_args()
    if args.threads < 1 or args.memory_gib < 2:
        raise SystemExit("threads must be >=1 and memory-gib must be >=2")

    import numpy as np
    import psi4
    from rdkit import Chem

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_output = args.output_dir / "o6u_hf_631gd_environment_canary.psi4.out"
    result_path = args.output_dir / "o6u_hf_631gd_environment_canary.json"
    mol = Chem.MolFromMolFile(str(args.sdf), removeHs=False, sanitize=True)
    if mol is None or mol.GetNumAtoms() != 76:
        raise SystemExit("SDF must contain the audited 76-atom O6U component")
    if Chem.GetFormalCharge(mol) != 0:
        raise SystemExit("O6U formal charge is not zero")
    conformer = mol.GetConformer()
    lines = ["0 1"]
    for atom in mol.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {position.x: .10f} {position.y: .10f} {position.z: .10f}")
    lines.extend(["units angstrom", "no_com", "no_reorient", "symmetry c1"])
    geometry = "\n".join(lines)

    os.environ.setdefault("PSI_SCRATCH", "/root/autodl-tmp/psi4_scratch")
    Path(os.environ["PSI_SCRATCH"]).mkdir(parents=True, exist_ok=True)
    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gib} GiB")
    psi4.core.set_output_file(str(raw_output), False)
    psi4.set_options(
        {
            "basis": "6-31G(d)",
            "reference": "rhf",
            "scf_type": "df",
            "guess": "sad",
            "e_convergence": 1.0e-8,
            "d_convergence": 1.0e-8,
            "maxiter": 200,
        }
    )
    started = time.time()
    record = {
        "schema_version": "1.0",
        "role": "runtime_canary_not_parameter_target_data",
        "production_approved": False,
        "method": "HF/6-31G(d) density-fitted single-point",
        "input_sdf": {"path": str(args.sdf.resolve()), "sha256": sha256(args.sdf)},
        "charge": 0,
        "multiplicity": 1,
        "atom_count": mol.GetNumAtoms(),
        "threads": args.threads,
        "memory_gib": args.memory_gib,
        "psi_scratch": os.environ["PSI_SCRATCH"],
        "psi4_version": psi4.__version__,
        "python_version": platform.python_version(),
        "host": platform.node(),
        "status": "running",
    }
    try:
        energy, wavefunction = psi4.energy("hf", molecule=psi4.geometry(geometry), return_wfn=True)
        dipole = np.asarray(wavefunction.variable("SCF DIPOLE"), dtype=float).reshape(-1).tolist()
        record.update(
            {
                "status": "pass_runtime_canary_only",
                "energy_hartree": float(energy),
                "scf_dipole_debye": dipole,
            }
        )
    except Exception as exc:
        record.update({"status": "fail", "error_type": type(exc).__name__, "error": str(exc)})
    finally:
        record["elapsed_seconds"] = time.time() - started
        psi4.core.flush_outfile()
        if raw_output.is_file():
            record["raw_output"] = {"path": str(raw_output.resolve()), "bytes": raw_output.stat().st_size, "sha256": sha256(raw_output)}
        result_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"status": record["status"], "elapsed_seconds": record["elapsed_seconds"], "result_sha256": sha256(result_path)}, sort_keys=True))
    return 0 if record["status"] == "pass_runtime_canary_only" else 1


if __name__ == "__main__":
    raise SystemExit(main())
