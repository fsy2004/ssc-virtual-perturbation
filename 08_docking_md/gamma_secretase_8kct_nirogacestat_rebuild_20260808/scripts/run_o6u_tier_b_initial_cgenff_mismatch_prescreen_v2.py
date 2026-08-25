#!/usr/bin/env python3
"""Compare initial-CGenFF energies for the Tier-B reference and displaced geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.5094740631


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def read_xyz(path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    lines = path.read_text().splitlines()
    count = int(lines[0].strip())
    rows = [line.split() for line in lines[2 : 2 + count]]
    if len(rows) != count:
        raise ValueError(f"XYZ atom count mismatch in {path}")
    symbols = [row[0] for row in rows]
    positions = [(float(row[1]), float(row[2]), float(row[3])) for row in rows]
    return symbols, positions


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--psf", required=True)
    parser.add_argument("--base-rtf", required=True)
    parser.add_argument("--base-prm", required=True)
    parser.add_argument("--ligand-rtf", required=True)
    parser.add_argument("--ligand-prm", required=True)
    parser.add_argument("--reference-xyz", required=True)
    parser.add_argument("--displaced-xyz", required=True)
    parser.add_argument("--reference-qm-energy-hartree", type=float, required=True)
    parser.add_argument("--displaced-qm-report", required=True)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--signed-step-index", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    inputs = {
        name: Path(value).resolve()
        for name, value in {
            "psf": args.psf,
            "base_rtf": args.base_rtf,
            "base_prm": args.base_prm,
            "ligand_rtf": args.ligand_rtf,
            "ligand_prm": args.ligand_prm,
            "reference_xyz": args.reference_xyz,
            "displaced_xyz": args.displaced_xyz,
            "displaced_qm_report": args.displaced_qm_report,
        }.items()
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    import openmm
    from openmm import Context, Platform, VerletIntegrator, unit
    from openmm.app import CharmmParameterSet, CharmmPsfFile, NoCutoff

    psf = CharmmPsfFile(str(inputs["psf"]))
    parameters = CharmmParameterSet(
        str(inputs["base_rtf"]),
        str(inputs["base_prm"]),
        str(inputs["ligand_rtf"]),
        str(inputs["ligand_prm"]),
    )
    system = psf.createSystem(
        parameters,
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
    )
    force_names: dict[int, str] = {}
    if system.getNumForces() > 31:
        raise RuntimeError("Too many forces for force-group decomposition")
    for index, force in enumerate(system.getForces()):
        force.setForceGroup(index)
        force_names[index] = force.__class__.__name__

    reference_symbols, reference_positions = read_xyz(inputs["reference_xyz"])
    displaced_symbols, displaced_positions = read_xyz(inputs["displaced_xyz"])
    psf_symbols = [atom.element.symbol for atom in psf.topology.atoms()]
    if reference_symbols != displaced_symbols or reference_symbols != psf_symbols:
        raise ValueError("PSF and XYZ atom order/element identity do not match exactly")

    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    selected_platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, selected_platform)

    def evaluate(positions: list[tuple[float, float, float]]) -> dict[str, object]:
        context.setPositions(positions * unit.angstrom)
        total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilocalories_per_mole
        )
        components: dict[str, float] = {}
        for index, name in force_names.items():
            value = context.getState(
                getEnergy=True, groups=1 << index
            ).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            components[f"{index}:{name}"] = float(value)
        if not math.isfinite(total) or not all(math.isfinite(x) for x in components.values()):
            raise ValueError("Non-finite MM energy")
        return {"total_kcal_mol": float(total), "components_kcal_mol": components}

    reference_mm = evaluate(reference_positions)
    displaced_mm = evaluate(displaced_positions)
    del context, integrator

    displaced_qm_report = json.loads(inputs["displaced_qm_report"].read_text())
    if displaced_qm_report.get("status") != "pass_relaxed_mp2_torsion_scan_point":
        raise ValueError("Displaced QM point is not a passed scan point")
    if displaced_qm_report.get("rotor_id") != args.rotor_id:
        raise ValueError("Displaced QM rotor identity mismatch")
    if int(displaced_qm_report.get("signed_step_index")) != args.signed_step_index:
        raise ValueError("Displaced QM signed-step mismatch")
    displaced_qm_energy = float(displaced_qm_report["final_energy_hartree"])
    qm_delta = (
        displaced_qm_energy - args.reference_qm_energy_hartree
    ) * HARTREE_TO_KCAL_MOL
    mm_delta = float(displaced_mm["total_kcal_mol"]) - float(
        reference_mm["total_kcal_mol"]
    )
    delta_error = mm_delta - qm_delta

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_tier_b_initial_cgenff_mismatch_prescreen",
        "status": "pass_tier_b_initial_cgenff_mismatch_prescreen",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "compound": "nirogacestat/PF-03084014/BRD-K61691541",
            "residue": "O6U",
            "rotor_id": args.rotor_id,
            "signed_step_index": args.signed_step_index,
            "purpose": "Tier-B principal-coordinate selection only; no parameter fitting or mutation",
        },
        "inputs": {name: record(path) for name, path in inputs.items()},
        "method": {
            "force_field": "initial CGenFF 5.0 plus official O6U ligand-reader parameters",
            "nonbonded_method": "NoCutoff",
            "constraints": None,
            "coordinate_minimization": False,
            "parameter_mutation": False,
            "openmm_platform": selected_platform.getName(),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "openmm": openmm.version.version,
        },
        "atom_identity": {"count": len(psf_symbols), "symbols": psf_symbols},
        "reference": {
            "qm_energy_hartree": args.reference_qm_energy_hartree,
            "mm": reference_mm,
        },
        "displaced": {
            "qm_energy_hartree": displaced_qm_energy,
            "mm": displaced_mm,
        },
        "comparison": {
            "qm_delta_kcal_mol": qm_delta,
            "initial_cgenff_mm_delta_kcal_mol": mm_delta,
            "signed_mm_minus_qm_delta_kcal_mol": delta_error,
            "absolute_mm_minus_qm_delta_kcal_mol": abs(delta_error),
        },
        "interpretation_boundary": (
            "This is a same-geometry initial-force-field prescreen for scan-scope selection. "
            "It does not validate parameters, affinity, efficacy, or production MD readiness."
        ),
    }
    report_path = output_dir / "O6U_TIER_B_INITIAL_CGENFF_MISMATCH_PRESCREEN.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(report_path)
    print(
        json.dumps(
            {
                "status": report["status"],
                "report": str(report_path),
                "sha256": sha256(report_path),
                "comparison": report["comparison"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
