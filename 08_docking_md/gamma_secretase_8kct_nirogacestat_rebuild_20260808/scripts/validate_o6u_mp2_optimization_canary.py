#!/usr/bin/env python3
"""Independently validate a completed O6U MP2 optimization target."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from prepare_o6u_crest_input import EXPECTED_CHIRAL_CENTERS, load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble, signed_chiral_volume


# Psi4 prints RHF, canonical DF-MP2, and SCS-MP2 energies in each gradient
# evaluation.  Only the unprefixed canonical MP2 line is the frozen model
# chemistry returned by psi4.optimize(); anchoring and requiring ``[Eh]``
# prevents an SCS energy from being mistaken for the canary energy.
TOTAL_ENERGY_RE = re.compile(
    r"^\s+Total Energy\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+\[Eh\]\s*$",
    re.MULTILINE,
)
EXPECTED_OPTIMIZER_STRATEGY = "cartesian_rfo_trust020_v1"
EXPECTED_OPTIMIZER_OPTIONS = {
    "opt_coordinates": "cartesian",
    "step_type": "rfo",
    "intrafrag_step_limit": 0.20,
    "intrafrag_step_limit_min": 0.01,
    "intrafrag_step_limit_max": 0.20,
    "dynamic_level": 0,
}
GAU_TIGHT_LIMITS = {
    "max_force": 1.50e-5,
    "rms_force": 1.00e-5,
    "max_displacement": 6.00e-5,
    "rms_displacement": 4.00e-5,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(message)


def parse_convergence_rows(raw: str) -> list[dict[str, float | int]]:
    """Reconstruct OptKing convergence rows without trusting status text."""
    rows: list[dict[str, float | int]] = []
    number_re = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")
    row_start = re.compile(r"^\s*\d+\s+-\d{3,}\.\d+\s+[+-]?\d+\.\d+[Ee][+-]\d+", re.IGNORECASE)
    for line in raw.splitlines():
        if not row_start.match(line):
            continue
        values = number_re.findall(line)
        if len(values) < 7:
            continue
        rows.append(
            {
                "step": int(values[0]),
                "energy_hartree": float(values[1]),
                "delta_energy_hartree": float(values[2]),
                "max_force_au": float(values[3]),
                "rms_force_au": float(values[4]),
                "max_displacement_au": float(values[5]),
                "rms_displacement_au": float(values[6]),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--expected-role",
        choices=("canary", "representative_target"),
        default="canary",
        help="Frozen role expected in the run record; defaults to the canary gate.",
    )
    args = parser.parse_args()

    record_path = args.record.resolve()
    source_path = args.source_sdf.resolve()
    report_path = args.report.resolve()
    for path in (record_path, source_path):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty required input: {path}")
    if report_path.exists():
        fail(f"Refusing to overwrite report: {report_path}")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    expected_status = {
        "canary": "pass_optimization_canary",
        "representative_target": "pass_optimization_representative_target",
    }[args.expected_role]
    if record.get("role") != args.expected_role:
        fail(f"Run-record role differs from expected role: {record.get('role')}")
    if record.get("status") != expected_status:
        fail(f"Optimization record is not a completed {args.expected_role} pass: {record.get('status')}")
    if record.get("production_approved") is not False:
        fail("Canary record improperly claims production approval")
    if record.get("charge_e") != 0 or record.get("multiplicity") != 1:
        fail("Canary record charge or multiplicity differs from frozen neutral singlet")
    if record.get("method") != "DF-MP2/6-31G(d), frozen core, RHF reference":
        fail("Canary record method differs from the frozen model chemistry")
    if record.get("psi4_version") != "1.9.1":
        fail(f"Unexpected Psi4 version: {record.get('psi4_version')}")
    if record.get("schema_version") != "1.1":
        fail(f"Unexpected repaired-canary schema: {record.get('schema_version')}")
    if record.get("optimizer_strategy_version") != EXPECTED_OPTIMIZER_STRATEGY:
        fail(f"Unexpected optimizer strategy: {record.get('optimizer_strategy_version')}")
    if record.get("optimizer_options") != EXPECTED_OPTIMIZER_OPTIONS:
        fail("Optimizer options differ from the repaired Cartesian/RFO canary")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    if record.get("source_sdf", {}).get("sha256") != sha256(source_path):
        fail("Source SDF hash does not match the run record")

    start_path = Path(record["start_xyz"]["path"])
    optimized_path = Path(record["optimized_xyz"]["path"])
    raw_path = Path(record["raw_output"])
    for path in (start_path, optimized_path, raw_path):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Declared canary artifact is missing or empty: {path}")
    if record["start_xyz"].get("sha256") != sha256(start_path):
        fail("Start XYZ hash does not match the run record")
    if record["optimized_xyz"].get("sha256") != sha256(optimized_path):
        fail("Optimized XYZ hash does not match the run record")
    if record.get("raw_output_sha256") != sha256(raw_path):
        fail("Raw Psi4 output hash does not match the run record")

    raw = raw_path.read_text(encoding="utf-8", errors="replace")
    if "Optimization is complete!" not in raw:
        fail("Raw Psi4 output lacks the independent optimization-complete marker")
    forbidden = (
        "Optimization failed",
        "Fatal Error",
        "PsiException",
        "Traceback (most recent call last)",
        "Back transformation failed",
        "Cartesian Step size too large",
        "Erasing history",
    )
    found_forbidden = [marker for marker in forbidden if marker in raw]
    if found_forbidden:
        fail(f"Raw Psi4 output contains failure markers: {found_forbidden}")
    energies = [float(value) for value in TOTAL_ENERGY_RE.findall(raw)]
    if not energies or not all(math.isfinite(value) for value in energies):
        fail("Raw Psi4 output does not contain a finite canonical DF-MP2 energy history")
    if len(energies) < 2 or len({round(value, 10) for value in energies}) < 2:
        fail("Energy history does not demonstrate progress away from the starting geometry")
    recorded_energy = float(record["final_energy_hartree"])
    if not math.isfinite(recorded_energy) or abs(energies[-1] - recorded_energy) > 1.0e-8:
        fail(
            f"Final energy mismatch: raw last {energies[-1]:.12f}, record {recorded_energy:.12f}"
        )
    convergence_rows = parse_convergence_rows(raw)
    if len(convergence_rows) < 2:
        fail("Raw Psi4 output lacks at least two independently parsed OptKing convergence rows")
    final_row = convergence_rows[-1]
    if abs(float(final_row["energy_hartree"]) - recorded_energy) > 1.0e-7:
        fail("Final OptKing convergence-row energy differs from the run record")
    observed = {
        "max_force": float(final_row["max_force_au"]),
        "rms_force": float(final_row["rms_force_au"]),
        "max_displacement": float(final_row["max_displacement_au"]),
        "rms_displacement": float(final_row["rms_displacement_au"]),
    }
    failed_limits = {
        key: {"observed": observed[key], "limit": limit}
        for key, limit in GAU_TIGHT_LIMITS.items()
        if not math.isfinite(observed[key]) or observed[key] > limit
    }
    if failed_limits:
        fail(f"Final OptKing row does not independently satisfy gau_tight: {failed_limits}")

    start_frames = read_xyz_ensemble(start_path)
    optimized_frames = read_xyz_ensemble(optimized_path)
    if len(start_frames) != 1 or len(optimized_frames) != 1:
        fail("Start and optimized XYZ files must each contain exactly one frame")
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if start_frames[0]["elements"] != expected_elements or optimized_frames[0]["elements"] != expected_elements:
        fail("Start or optimized element order differs from immutable source SDF")
    start_coordinates = np.asarray(start_frames[0]["coordinates"], dtype=float)
    optimized_coordinates = np.asarray(optimized_frames[0]["coordinates"], dtype=float)
    if not np.isfinite(optimized_coordinates).all():
        fail("Optimized XYZ contains non-finite coordinates")
    coordinate_rmsd = float(np.sqrt(np.mean((optimized_coordinates - start_coordinates) ** 2)))
    if coordinate_rmsd <= 1.0e-5:
        fail("Optimized geometry does not demonstrate finite coordinate progress")

    bonds = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in source.GetBonds())
    bond_ratios = []
    for i, j in bonds:
        initial = float(np.linalg.norm(start_coordinates[i] - start_coordinates[j]))
        final = float(np.linalg.norm(optimized_coordinates[i] - optimized_coordinates[j]))
        ratio = final / initial
        bond_ratios.append(ratio)
        if not 0.75 <= ratio <= 1.30:
            fail(f"Optimized bond {i}-{j} changed beyond topology guard: ratio {ratio:.6f}")

    chiral_review = {}
    for center in sorted(EXPECTED_CHIRAL_CENTERS):
        neighbours = sorted(neighbour.GetIdx() for neighbour in source.GetAtomWithIdx(center).GetNeighbors())
        initial_volume = signed_chiral_volume(start_coordinates, center, neighbours)
        final_volume = signed_chiral_volume(optimized_coordinates, center, neighbours)
        if abs(initial_volume) <= 1.0e-4 or abs(final_volume) <= 1.0e-4:
            fail(f"Chiral center {center} is geometrically degenerate")
        if math.copysign(1.0, initial_volume) != math.copysign(1.0, final_volume):
            fail(f"Chiral center {center} inverted during MP2 optimization")
        chiral_review[str(center)] = {
            "initial_signed_volume_angstrom3": initial_volume,
            "optimized_signed_volume_angstrom3": final_volume,
            "sign_preserved": True,
        }

    elapsed = float(record.get("elapsed_seconds", math.nan))
    if not math.isfinite(elapsed) or elapsed <= 0:
        fail("Run record has no finite positive elapsed time")
    report = {
        "schema_version": "1.1",
        "report_type": f"independent_o6u_mp2_optimization_{args.expected_role}_validation",
        "status": (
            "pass_canary_independently_reconstructed"
            if args.expected_role == "canary"
            else "pass_representative_independently_reconstructed"
        ),
        "role": args.expected_role,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "record": str(record_path),
        "record_sha256": sha256(record_path),
        "identity": identity,
        "raw_output": str(raw_path),
        "raw_output_sha256": sha256(raw_path),
        "raw_total_energy_observations": len(energies),
        "raw_distinct_total_energies_at_1e_minus_10_hartree": len({round(value, 10) for value in energies}),
        "final_energy_hartree": recorded_energy,
        "optking_convergence_rows": len(convergence_rows),
        "final_optking_convergence_row": final_row,
        "gau_tight_limits": GAU_TIGHT_LIMITS,
        "optimized_xyz": str(optimized_path),
        "optimized_xyz_sha256": sha256(optimized_path),
        "coordinate_rmsd_to_start_angstrom": coordinate_rmsd,
        "optimizer_strategy_version": EXPECTED_OPTIMIZER_STRATEGY,
        "minimum_bond_length_ratio_to_start": min(bond_ratios),
        "maximum_bond_length_ratio_to_start": max(bond_ratios),
        "chiral_review": chiral_review,
        "elapsed_seconds": elapsed,
        "release_boundary": (
            "This releases review of the remaining frozen QM representatives only; "
            "it does not approve fitted ligand parameters, a CHARMM-GUI system, or production MD."
            if args.expected_role == "canary"
            else "This validates one frozen QM representative only; parameter fitting and all MD stages remain blocked."
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report_sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
