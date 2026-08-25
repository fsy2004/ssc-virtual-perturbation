#!/usr/bin/env python3
"""Extract and validate the last complete O6U geometry printed by Psi4.

The resulting XYZ is a technical restart candidate only. It does not establish
optimization convergence and cannot release parameter fitting or MD.
"""

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


HEADER = "Geometry (in Angstrom), charge = 0, multiplicity = 1:"
ATOM_ROW = re.compile(
    r"^\s*([A-Z][a-z]?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def complete_geometry_blocks(raw: str) -> list[dict[str, object]]:
    lines = raw.splitlines()
    blocks: list[dict[str, object]] = []
    for header_index, line in enumerate(lines):
        if HEADER not in line:
            continue
        elements: list[str] = []
        coordinates: list[list[float]] = []
        for candidate in lines[header_index + 1 :]:
            match = ATOM_ROW.match(candidate)
            if match:
                elements.append(match.group(1))
                coordinates.append([float(match.group(2)), float(match.group(3)), float(match.group(4))])
                if len(elements) == EXPECTED_ATOMS:
                    blocks.append(
                        {
                            "header_line_1based": header_index + 1,
                            "elements": elements,
                            "coordinates": np.asarray(coordinates, dtype=float),
                        }
                    )
                    break
            elif elements:
                break
    return blocks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", choices=("snapshot_only", "failed_run_restart_candidate"), required=True)
    args = parser.parse_args()

    record_path = args.record.resolve()
    source_path = args.source_sdf.resolve()
    output_dir = args.output_dir.resolve()
    for path in (record_path, source_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("charge_e") != 0 or record.get("multiplicity") != 1:
        raise SystemExit("Run record is not the frozen neutral singlet")
    if record.get("method") != "DF-MP2/6-31G(d), frozen core, RHF reference":
        raise SystemExit("Run record model chemistry differs")
    if record.get("optimizer_convergence") != "gau_tight":
        raise SystemExit("Run record optimization convergence differs")

    raw_path = Path(str(record.get("raw_output", ""))).resolve()
    start_path = Path(str(record.get("start_xyz", {}).get("path", ""))).resolve()
    for path in (raw_path, start_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Declared run artifact is missing or empty: {path}")
    if record.get("start_xyz", {}).get("sha256") != sha256(start_path):
        raise SystemExit("Start XYZ hash differs from the run record")
    # Read the live output exactly once.  Hashing the path and then reading it
    # again can bind a later geometry block to an earlier hash while Psi4 is
    # appending.  A byte snapshot keeps the audit internally self-consistent.
    raw_bytes = raw_path.read_bytes()
    raw_hash = sha256_bytes(raw_bytes)
    if args.role == "failed_run_restart_candidate":
        if record.get("status") != "fail" or record.get("error_type") != "OptimizationConvergenceError":
            raise SystemExit("Technical restart role requires a recorded OptimizationConvergenceError")
        if record.get("raw_output_sha256") != raw_hash:
            raise SystemExit("Failed-run raw output hash differs from the closed run record")
        error_text = str(record.get("error", "")).lower()
        if "maximum number of steps" not in error_text and "geom_maxiter" not in error_text:
            raise SystemExit("Failed run did not stop at the prespecified geometry-iteration limit")
    elif record.get("status") != "running":
        raise SystemExit("Snapshot role is restricted to a running canary")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    start_frames = read_xyz_ensemble(start_path)
    if len(start_frames) != 1 or start_frames[0]["elements"] != expected_elements:
        raise SystemExit("Start XYZ identity/order differs from the immutable source")
    start_coordinates = np.asarray(start_frames[0]["coordinates"], dtype=float)

    raw = raw_bytes.decode("utf-8", errors="replace")
    blocks = complete_geometry_blocks(raw)
    matching = [block for block in blocks if block["elements"] == expected_elements]
    if not matching:
        raise SystemExit("Raw Psi4 output contains no complete identity-matching geometry block")
    selected = matching[-1]
    coordinates = np.asarray(selected["coordinates"], dtype=float)
    if coordinates.shape != (EXPECTED_ATOMS, 3) or not np.isfinite(coordinates).all():
        raise SystemExit("Selected geometry is non-finite or has the wrong shape")

    bonds = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in source.GetBonds())
    bond_ratios: list[float] = []
    for left, right in bonds:
        initial = float(np.linalg.norm(start_coordinates[left] - start_coordinates[right]))
        current = float(np.linalg.norm(coordinates[left] - coordinates[right]))
        ratio = current / initial
        bond_ratios.append(ratio)
        if not 0.75 <= ratio <= 1.30:
            raise SystemExit(f"Restart geometry violates topology guard for bond {left}-{right}: {ratio:.6f}")

    chiral_review: dict[str, dict[str, float | bool]] = {}
    for center in sorted(EXPECTED_CHIRAL_CENTERS):
        neighbours = sorted(neighbour.GetIdx() for neighbour in source.GetAtomWithIdx(center).GetNeighbors())
        initial_volume = signed_chiral_volume(start_coordinates, center, neighbours)
        current_volume = signed_chiral_volume(coordinates, center, neighbours)
        if abs(initial_volume) <= 1.0e-4 or abs(current_volume) <= 1.0e-4:
            raise SystemExit(f"Chiral center {center} is geometrically degenerate")
        if math.copysign(1.0, initial_volume) != math.copysign(1.0, current_volume):
            raise SystemExit(f"Chiral center {center} inverted")
        chiral_review[str(center)] = {
            "start_signed_volume_angstrom3": initial_volume,
            "restart_signed_volume_angstrom3": current_volume,
            "sign_preserved": True,
        }

    output_dir.mkdir(parents=True, exist_ok=False)
    xyz_path = output_dir / "o6u_last_complete_psi4_geometry.restart_candidate.xyz"
    xyz_lines = [str(EXPECTED_ATOMS), "technical restart candidate; not a converged geometry"]
    xyz_lines.extend(
        f"{element:<2s} {x: .12f} {y: .12f} {z: .12f}"
        for element, (x, y, z) in zip(expected_elements, coordinates, strict=True)
    )
    xyz_path.write_text("\n".join(xyz_lines) + "\n", encoding="utf-8", newline="\n")
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_psi4_last_complete_geometry_extraction",
        "status": "pass_technical_restart_candidate_not_converged",
        "role": args.role,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "record": {"path": str(record_path), "sha256": sha256(record_path), "status": record.get("status")},
        "raw_output": {"path": str(raw_path), "sha256_at_extraction": raw_hash, "bytes_at_extraction": len(raw_bytes)},
        "complete_geometry_block_count": len(blocks),
        "identity_matching_geometry_block_count": len(matching),
        "selected_header_line_1based": selected["header_line_1based"],
        "restart_xyz": {"path": str(xyz_path), "sha256": sha256(xyz_path)},
        "coordinate_rmsd_to_start_angstrom": float(np.sqrt(np.mean((coordinates - start_coordinates) ** 2))),
        "minimum_bond_length_ratio_to_start": min(bond_ratios),
        "maximum_bond_length_ratio_to_start": max(bond_ratios),
        "chiral_review": chiral_review,
        "release_boundary": (
            "This file is only an identity/topology-validated technical restart candidate. It does not pass "
            "gau_tight, approve an MP2 representative, authorize parameter fitting, or authorize MD."
        ),
    }
    report_path = output_dir / "O6U_PSI4_RESTART_GEOMETRY_EXTRACTION.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "role": args.role, "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
