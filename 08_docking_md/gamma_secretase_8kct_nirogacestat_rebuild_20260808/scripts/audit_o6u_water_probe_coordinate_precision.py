#!/usr/bin/env python3
"""Explain PDB-versus-CRD water-probe differences and freeze the CRD plan."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_PDB_SHA256 = "33ae064817954275aff2b42dcfcf0597be6e15ce71d56cbec39267532c4434b7"
EXPECTED_CRD_SHA256 = "9d18e691f3e3afb29d5bb18584a6245b7e1cf3aefb2e660b2bcf08dee6857167"
EXPECTED_PDB_DA_SHA256 = "537b65a97ad9032aa2cf496ae686f3c62815a52e82f61c746d3726cdc5c1d7b6"
EXPECTED_CRD_DA_SHA256 = "5ea7e12c464750b8f35d9fa3feed875bb6b69a4091867b8fb6ddf2ad3c6272a3"
EXPECTED_TYPES = {"A2": 2, "A31": 16, "AP": 2, "APL": 6, "D": 38, "DOP": 6}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_pdb(path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(
                (
                    line[12:16].strip(),
                    (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                )
            )
    return atoms


def read_crd(path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 7 and fields[0].isdigit() and fields[1].isdigit():
            atoms.append((fields[3], tuple(float(value) for value in fields[4:7])))
    return atoms


def read_da(path: Path) -> list[list[str]]:
    return [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--crd", required=True, type=Path)
    parser.add_argument("--pdb-da", required=True, type=Path)
    parser.add_argument("--crd-da", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    paths = {key: getattr(args, key).resolve() for key in ("pdb", "crd", "pdb_da", "crd_da")}
    expected_hashes = {
        "pdb": EXPECTED_PDB_SHA256,
        "crd": EXPECTED_CRD_SHA256,
        "pdb_da": EXPECTED_PDB_DA_SHA256,
        "crd_da": EXPECTED_CRD_DA_SHA256,
    }
    for key, path in paths.items():
        if not path.is_file() or sha256(path) != expected_hashes[key]:
            raise SystemExit(f"{key} input identity differs: {path}")

    pdb_atoms = read_pdb(paths["pdb"])
    crd_atoms = read_crd(paths["crd"])
    if len(pdb_atoms) != 76 or len(crd_atoms) != 76:
        raise SystemExit(f"Expected 76 atoms, found PDB={len(pdb_atoms)} CRD={len(crd_atoms)}")
    if [item[0] for item in pdb_atoms] != [item[0] for item in crd_atoms]:
        raise SystemExit("PDB and CRD atom names/order differ")
    coordinate_differences = [
        abs(pdb_value - crd_value)
        for (_, pdb_xyz), (_, crd_xyz) in zip(pdb_atoms, crd_atoms, strict=True)
        for pdb_value, crd_value in zip(pdb_xyz, crd_xyz, strict=True)
    ]
    maximum_coordinate_difference = max(coordinate_differences)
    rms_coordinate_difference = math.sqrt(sum(value * value for value in coordinate_differences) / len(coordinate_differences))
    if maximum_coordinate_difference > 0.0005001:
        raise SystemExit("PDB/CRD coordinate difference exceeds three-decimal rounding")

    pdb_da = read_da(paths["pdb_da"])
    crd_da = read_da(paths["crd_da"])
    if len(pdb_da) != 70 or len(crd_da) != 70:
        raise SystemExit(f"Expected 70 orientations, found PDB={len(pdb_da)} CRD={len(crd_da)}")
    if dict(Counter(row[0] for row in pdb_da)) != EXPECTED_TYPES or dict(Counter(row[0] for row in crd_da)) != EXPECTED_TYPES:
        raise SystemExit("Water-probe type counts differ from the frozen plan")
    changed_rows = []
    maximum_angle_difference = 0.0
    for index, (pdb_row, crd_row) in enumerate(zip(pdb_da, crd_da, strict=True), start=1):
        numeric_start = 5 if crd_row[0] in {"A31", "AP"} else len(crd_row)
        if pdb_row[:numeric_start] != crd_row[:numeric_start] or len(pdb_row) != len(crd_row):
            raise SystemExit(f"Orientation identity differs at row {index}")
        if numeric_start == len(crd_row):
            continue
        numeric_differences = [
            abs(float(left) - float(right))
            for left, right in zip(pdb_row[numeric_start:], crd_row[numeric_start:], strict=True)
        ]
        if any(value > 0 for value in numeric_differences):
            maximum_angle_difference = max(maximum_angle_difference, *numeric_differences)
            changed_rows.append(
                {
                    "source_line_number": index,
                    "probe_type": crd_row[0],
                    "target_atom": crd_row[2],
                    "pdb_definition": "  ".join(pdb_row),
                    "crd_definition": "  ".join(crd_row),
                    "maximum_numeric_difference_degrees": max(numeric_differences),
                }
            )

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_coordinate_precision_audit",
        "status": "pass_crd_high_precision_formal_pdb_rounding_explained",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            key: {"path": str(path), "sha256": expected_hashes[key], "size_bytes": path.stat().st_size}
            for key, path in paths.items()
        },
        "coordinate_atom_count": 76,
        "atom_names_and_order_identical": True,
        "maximum_pdb_crd_coordinate_difference_angstrom": maximum_coordinate_difference,
        "rms_pdb_crd_coordinate_difference_angstrom": rms_coordinate_difference,
        "difference_explained_by_pdb_three_decimal_rounding": True,
        "orientation_count": 70,
        "orientation_identity_fields_identical": True,
        "orientation_type_counts": EXPECTED_TYPES,
        "numeric_angle_difference_row_count": len(changed_rows),
        "maximum_angle_difference_degrees": maximum_angle_difference,
        "changed_rows": changed_rows,
        "formal_frozen_orientation_source": {
            "format": "CHARMM CRD full-precision coordinates",
            "orientation_da_sha256": EXPECTED_CRD_DA_SHA256,
        },
        "historical_canary_disposition": (
            "The PDB-derived plan remains archived as a headless environment canary only. "
            "It is not an alternative QM target set."
        ),
        "release_boundary": "Orientation provenance only; no water QM target or ligand parameter is approved.",
    }
    output = args.report.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report_sha256": sha256(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
