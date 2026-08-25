#!/usr/bin/env python3
"""Independently validate the completed O6U geomeTRIC/TRIC MP2 canary."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from prepare_o6u_crest_input import EXPECTED_ATOMS, EXPECTED_CHIRAL_CENTERS, load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble, signed_chiral_volume


ENERGY_RE = re.compile(r"^\s+Total Energy\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+\[Eh\]\s*$", re.M)
NUMBER_RE = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")
LIMITS = {"max_force": 1.5e-5, "rms_force": 1.0e-5, "max_displacement": 6.0e-5, "rms_displacement": 4.0e-5}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(message)


def convergence_rows(raw: str) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in raw.splitlines():
        if "~" not in line:
            continue
        values = NUMBER_RE.findall(line)
        if len(values) != 7:
            continue
        rows.append({
            "step": int(values[0]),
            "energy_hartree": float(values[1]),
            "delta_energy_hartree": float(values[2]),
            "max_force_au": float(values[3]),
            "rms_force_au": float(values[4]),
            "max_displacement_au": float(values[5]),
            "rms_displacement_au": float(values[6]),
        })
    return rows


def final_gradient(raw: str) -> tuple[float, float, int]:
    lines = raw.splitlines()
    found: list[tuple[float, float]] = []
    for start, line in enumerate(lines):
        if "-Total Gradient:" not in line:
            continue
        values: list[float] = []
        for row in lines[start + 3 : start + 3 + EXPECTED_ATOMS]:
            fields = row.split()
            if len(fields) < 4:
                break
            try:
                values.extend(float(value) for value in fields[-3:])
            except ValueError:
                break
        if len(values) == EXPECTED_ATOMS * 3:
            found.append((max(abs(value) for value in values), math.sqrt(sum(value * value for value in values) / len(values))))
    if not found:
        fail("Raw output contains no complete gradient")
    return found[-1][0], found[-1][1], len(found)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    record_path = args.record.resolve()
    source_path = args.source_sdf.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        fail(f"Refusing to reuse output directory: {output_dir}")
    for path in (record_path, source_path):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty input: {path}")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    required = {
        "status": "pass_optimization_geometric_recovery_canary",
        "role": "recovery_canary",
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "optimizer_convergence": "gau_tight",
        "optimizer_strategy_version": "geometric_tric_recovery_v1",
        "optimizer_engine": "geometric",
        "geometric_version": "1.1.1",
        "psi4_version": "1.9.1",
        "charge_e": 0,
        "multiplicity": 1,
    }
    for key, expected in required.items():
        if record.get(key) != expected:
            fail(f"Run record field {key} differs: {record.get(key)!r}")
    if record.get("optimizer_keywords") != {"coordsys": "tric"} or record.get("production_approved") is not False:
        fail("Optimizer keywords or production boundary differs")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    if record.get("identity") != identity or record.get("source_sdf", {}).get("sha256") != sha256(source_path):
        fail("Source identity or hash differs")
    start_path = Path(record["start_xyz"]["path"])
    optimized_path = Path(record["optimized_xyz"]["path"])
    raw_path = Path(record["raw_output"])
    for path in (start_path, optimized_path, raw_path):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Declared artifact missing or empty: {path}")
    if record["start_xyz"].get("sha256") != sha256(start_path):
        fail("Start XYZ hash differs")
    if record["optimized_xyz"].get("sha256") != sha256(optimized_path):
        fail("Optimized XYZ hash differs")
    if record.get("raw_output_sha256") != sha256(raw_path) or record.get("raw_output_bytes") != raw_path.stat().st_size:
        fail("Raw output closure differs")

    raw = raw_path.read_text(encoding="utf-8", errors="replace")
    if "Optimization converged!" not in raw or "Final Energy" not in raw:
        fail("Raw output lacks geomeTRIC convergence markers")
    forbidden = [marker for marker in ("Optimization failed", "Fatal Error", "Traceback (most recent call last)") if marker in raw]
    if forbidden:
        fail(f"Raw output contains failure markers: {forbidden}")
    energies = [float(value) for value in ENERGY_RE.findall(raw)]
    if len(energies) < 2 or not all(math.isfinite(value) for value in energies):
        fail("Canonical MP2 energy history is incomplete")
    recorded_energy = float(record["final_energy_hartree"])
    if abs(energies[-1] - recorded_energy) > 1.0e-8 or abs(float(record["wavefunction_energy_hartree"]) - recorded_energy) > 1.0e-10:
        fail("Final energy cross-check failed")
    rows = convergence_rows(raw)
    if len(rows) < 2:
        fail("geomeTRIC convergence table is incomplete")
    final_row = rows[-1]
    observed = {
        "max_force": float(final_row["max_force_au"]),
        "rms_force": float(final_row["rms_force_au"]),
        "max_displacement": float(final_row["max_displacement_au"]),
        "rms_displacement": float(final_row["rms_displacement_au"]),
    }
    failed = {key: {"observed": observed[key], "limit": limit} for key, limit in LIMITS.items() if observed[key] > limit}
    if failed:
        fail(f"Final reconstructed convergence row fails gau_tight: {failed}")
    gradient_max, gradient_rms, gradient_count = final_gradient(raw)
    if gradient_max > LIMITS["max_force"] or gradient_rms > LIMITS["rms_force"]:
        fail("Final Cartesian gradient independently fails gau_tight force limits")

    start_frames, optimized_frames = read_xyz_ensemble(start_path), read_xyz_ensemble(optimized_path)
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(start_frames) != 1 or len(optimized_frames) != 1 or start_frames[0]["elements"] != expected_elements or optimized_frames[0]["elements"] != expected_elements:
        fail("XYZ identity/order differs")
    start = np.asarray(start_frames[0]["coordinates"], dtype=float)
    optimized = np.asarray(optimized_frames[0]["coordinates"], dtype=float)
    if not np.isfinite(optimized).all():
        fail("Optimized coordinates are non-finite")
    bond_ratios: list[float] = []
    for bond in source.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ratio = float(np.linalg.norm(optimized[left] - optimized[right]) / np.linalg.norm(start[left] - start[right]))
        if not 0.75 <= ratio <= 1.30:
            fail(f"Topology guard failed for bond {left}-{right}: {ratio}")
        bond_ratios.append(ratio)
    chirality: dict[str, object] = {}
    for center in sorted(EXPECTED_CHIRAL_CENTERS):
        neighbours = sorted(atom.GetIdx() for atom in source.GetAtomWithIdx(center).GetNeighbors())
        before = signed_chiral_volume(start, center, neighbours)
        after = signed_chiral_volume(optimized, center, neighbours)
        if abs(before) <= 1.0e-4 or abs(after) <= 1.0e-4 or math.copysign(1.0, before) != math.copysign(1.0, after):
            fail(f"Chirality guard failed at center {center}")
        chirality[str(center)] = {"start_signed_volume_angstrom3": before, "optimized_signed_volume_angstrom3": after, "sign_preserved": True}

    output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "1.0",
        "report_type": "independent_o6u_mp2_geometric_recovery_validation",
        "status": "pass_geometric_recovery_independently_reconstructed_pending_minimum_character",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "record": {"path": str(record_path), "sha256": sha256(record_path)},
        "identity": identity,
        "raw_output": {"path": str(raw_path), "sha256": sha256(raw_path), "bytes": raw_path.stat().st_size},
        "optimized_xyz": {"path": str(optimized_path), "sha256": sha256(optimized_path)},
        "energy_observations": len(energies),
        "final_energy_hartree": recorded_energy,
        "convergence_rows": len(rows),
        "final_convergence_row": final_row,
        "gau_tight_limits": LIMITS,
        "final_cartesian_gradient": {"max_hartree_per_bohr": gradient_max, "rms_hartree_per_bohr": gradient_rms, "complete_gradient_count": gradient_count},
        "coordinate_rmsd_to_start_angstrom": float(np.sqrt(np.mean((optimized - start) ** 2))),
        "minimum_bond_length_ratio_to_start": min(bond_ratios),
        "maximum_bond_length_ratio_to_start": max(bond_ratios),
        "chirality": chirality,
        "minimum_character_validated": False,
        "release_boundary": "Identity, topology, chirality, energy, and gau_tight convergence pass; minimum-character validation remains required before releasing downstream parameter fitting.",
    }
    out = output_dir / "O6U_MP2_GEOMETRIC_RECOVERY_INDEPENDENT_VALIDATION.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report": str(out), "sha256": sha256(out)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
