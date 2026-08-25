#!/usr/bin/env python3
"""Build a whole-molecule exact-signature torsion correction design matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_o6u_rotor_to_charmm_dihedrals_v1 import parse_psf_text, sha256


def canonical_signature(types: list[str]) -> tuple[str, str, str, str]:
    forward = tuple(types)
    reverse = tuple(reversed(types))
    return min(forward, reverse)


def read_xyz(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2 : 2 + count]]
    if len(rows) != count:
        raise ValueError(f"XYZ atom count mismatch: {path}")
    return np.asarray([[float(row[1]), float(row[2]), float(row[3])] for row in rows])


def torsion_angle_deg(coordinates: np.ndarray, atom_ids_one_based: tuple[int, int, int, int]) -> float:
    p0, p1, p2, p3 = (coordinates[index - 1] for index in atom_ids_one_based)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def periodic_basis(angle_deg: float, reference_deg: float, periodicity: int, phase_deg: float) -> float:
    angle = math.radians(angle_deg)
    reference = math.radians(reference_deg)
    phase = math.radians(phase_deg)
    return math.cos(periodicity * angle - phase) - math.cos(periodicity * reference - phase)


def parse_parameter_line(line: str) -> dict:
    tokens = line.split()
    if len(tokens) < 7:
        raise ValueError(f"Invalid torsion parameter line: {line}")
    return {
        "pattern": tokens[:4],
        "force_constant_kcal_mol": float(tokens[4]),
        "periodicity": int(tokens[5]),
        "phase_deg": float(tokens[6]),
        "raw_line": line,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psf", required=True)
    parser.add_argument("--reference-xyz", required=True)
    parser.add_argument("--point-manifest", required=True)
    parser.add_argument("--mapping", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    psf_path = Path(args.psf).resolve()
    reference_path = Path(args.reference_xyz).resolve()
    point_manifest_path = Path(args.point_manifest).resolve()
    mapping_paths = [Path(value).resolve() for value in args.mapping]
    point_manifest = json.loads(point_manifest_path.read_text())
    if point_manifest.get("status") != "pass_hash_bound_global_design_matrix_inputs":
        raise ValueError("Point manifest did not pass")
    atoms, dihedrals = parse_psf_text(psf_path.read_text(errors="replace"))
    signatures_by_dihedral = {
        ids: canonical_signature([atoms[value]["atom_type"] for value in ids]) for ids in dihedrals
    }
    candidate_sources: dict[tuple[tuple[str, str, str, str], int, float], set[str]] = defaultdict(set)
    mapping_records = []
    for path in mapping_paths:
        report = json.loads(path.read_text())
        if report.get("status") != "pass_ccd_correspondence_bound_rotor_to_charmm_dihedral_mapping":
            raise ValueError(f"Mapping did not pass: {path}")
        mapping_records.append({"path": str(path), "sha256": sha256(path), "rotor_id": report["rotor_id"]})
        for row in report["incident_dihedrals"]:
            signature = canonical_signature(row["atom_types"])
            matches = row["matching_parameter_lines"]["ligand"] or row["matching_parameter_lines"]["base"]
            for line in matches:
                term = parse_parameter_line(line)
                candidate_sources[(signature, term["periodicity"], term["phase_deg"])].add(line)
    if not candidate_sources:
        raise ValueError("No active parameter harmonics found for mapped signatures")
    candidates = []
    for index, (key, sources) in enumerate(sorted(candidate_sources.items()), start=1):
        signature, periodicity, phase = key
        matching_dihedrals = [list(ids) for ids, found in signatures_by_dihedral.items() if found == signature]
        candidates.append({
            "column": index - 1,
            "canonical_atom_type_signature": list(signature),
            "periodicity": periodicity,
            "phase_deg": phase,
            "source_parameter_lines": sorted(sources),
            "whole_molecule_dihedral_count": len(matching_dihedrals),
            "whole_molecule_dihedrals": matching_dihedrals,
        })
    reference = read_xyz(reference_path)
    reference_angles = {ids: torsion_angle_deg(reference, ids) for ids in dihedrals}
    matrix_rows = []
    responses = []
    point_records = []
    for point in point_manifest["points"]:
        path = Path(point["xyz_path"]).resolve()
        if sha256(path) != point["xyz_sha256"]:
            raise ValueError(f"Point XYZ hash mismatch: {point['key']}")
        coordinates = read_xyz(path)
        if len(coordinates) != len(atoms):
            raise ValueError(f"Point atom count mismatch: {point['key']}")
        angles = {ids: torsion_angle_deg(coordinates, ids) for ids in dihedrals}
        row = []
        for candidate in candidates:
            signature = tuple(candidate["canonical_atom_type_signature"])
            value = sum(
                periodic_basis(angles[ids], reference_angles[ids], candidate["periodicity"], candidate["phase_deg"])
                for ids in dihedrals if signatures_by_dihedral[ids] == signature
            )
            row.append(value)
        matrix_rows.append(row)
        responses.append(float(point["qm_minus_initial_cgenff_delta_kcal_mol"]))
        point_records.append({key: point[key] for key in ("key", "rotor_id", "signed_step_index", "xyz_path", "xyz_sha256")})
    matrix = np.asarray(matrix_rows, dtype=float)
    response = np.asarray(responses, dtype=float)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0 else None
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_global_torsion_design_matrix",
        "status": "pass_global_torsion_design_matrix_built",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "psf": {"path": str(psf_path), "sha256": sha256(psf_path)},
            "reference_xyz": {"path": str(reference_path), "sha256": sha256(reference_path)},
            "point_manifest": {"path": str(point_manifest_path), "sha256": sha256(point_manifest_path)},
            "mapping_reports": mapping_records,
        },
        "point_count": len(point_records),
        "candidate_column_count": len(candidates),
        "matrix_rank": rank,
        "condition_number": condition,
        "singular_values": singular.tolist(),
        "identifiable_without_regularization": rank == len(candidates),
        "candidates": candidates,
        "points": point_records,
        "design_matrix": matrix.tolist(),
        "response_qm_minus_initial_cgenff_kcal_mol": response.tolist(),
        "parameter_mutation": False,
        "production_md_approved": False,
        "interpretation_boundary": (
            "Design-matrix and identifiability preparation only. Columns are exact-signature "
            "corrections using existing active periodicities/phases; no force constants are fitted or promoted."
        ),
    }
    output_dir.mkdir(parents=True)
    path = output_dir / "O6U_GLOBAL_TORSION_DESIGN_MATRIX.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    print(json.dumps({
        "status": report["status"],
        "point_count": report["point_count"],
        "candidate_column_count": report["candidate_column_count"],
        "matrix_rank": rank,
        "condition_number": condition,
        "identifiable_without_regularization": report["identifiable_without_regularization"],
        "sha256": sha256(path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
