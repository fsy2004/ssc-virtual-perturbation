#!/usr/bin/env python3
"""Analyze whether scanned O6U torsion type signatures collide elsewhere."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_o6u_rotor_to_charmm_dihedrals_v1 import parse_psf_text, sha256


def canonical_signature(types: list[str]) -> tuple[str, str, str, str]:
    forward = tuple(types)
    reverse = tuple(reversed(types))
    return min(forward, reverse)


def has_nonlocal_collision(occurrences: list[dict], target_bonds: set[frozenset[int]]) -> bool:
    return any(frozenset(row["central_bond"]) not in target_bonds for row in occurrences)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psf", required=True)
    parser.add_argument("--mapping", action="append", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    psf_path = Path(args.psf).resolve()
    mapping_paths = [Path(value).resolve() for value in args.mapping]
    atoms, dihedrals = parse_psf_text(psf_path.read_text(errors="replace"))
    occurrences: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for ids in dihedrals:
        types = [atoms[value]["atom_type"] for value in ids]
        signature = canonical_signature(types)
        occurrences[signature].append({
            "atom_ids": list(ids),
            "atom_names": [atoms[value]["atom_name"] for value in ids],
            "atom_types": types,
            "central_bond": [ids[1], ids[2]],
            "central_bond_atom_names": [atoms[ids[1]]["atom_name"], atoms[ids[2]]["atom_name"]],
        })
    mappings = []
    target_bonds_by_signature: dict[tuple[str, str, str, str], set[frozenset[int]]] = defaultdict(set)
    rotors_by_signature: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for path in mapping_paths:
        report = json.loads(path.read_text())
        if report.get("status") != "pass_ccd_correspondence_bound_rotor_to_charmm_dihedral_mapping":
            raise ValueError(f"Mapping report did not pass: {path}")
        mappings.append({"path": str(path), "sha256": sha256(path), "rotor_id": report["rotor_id"]})
        for row in report["incident_dihedrals"]:
            signature = canonical_signature(row["atom_types"])
            target_bonds_by_signature[signature].add(frozenset(row["atom_ids"][1:3]))
            rotors_by_signature[signature].add(report["rotor_id"])
    rows = []
    for signature in sorted(target_bonds_by_signature):
        found = occurrences[signature]
        rows.append({
            "canonical_atom_type_signature": list(signature),
            "scanned_rotors": sorted(rotors_by_signature[signature]),
            "target_central_bonds": [sorted(value) for value in sorted(target_bonds_by_signature[signature], key=lambda item: sorted(item))],
            "global_occurrence_count": len(found),
            "global_occurrences": found,
            "shared_across_scanned_rotors": len(rotors_by_signature[signature]) > 1,
            "nonlocal_central_bond_collision": has_nonlocal_collision(found, target_bonds_by_signature[signature]),
        })
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_torsion_type_conflict_analysis",
        "status": "pass_torsion_type_conflict_analysis",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "psf": {"path": str(psf_path), "sha256": sha256(psf_path), "size_bytes": psf_path.stat().st_size},
        "mapping_reports": mappings,
        "signature_count": len(rows),
        "signatures": rows,
        "nonlocal_collision_count": sum(row["nonlocal_central_bond_collision"] for row in rows),
        "shared_signature_count": sum(row["shared_across_scanned_rotors"] for row in rows),
        "parameter_mutation": False,
        "production_md_approved": False,
        "interpretation_boundary": (
            "Conflict analysis only. A nonlocal collision means a direct atom-type parameter "
            "override would affect additional torsions and requires joint fitting or retyping review."
        ),
    }
    output_dir.mkdir(parents=True)
    path = output_dir / "O6U_TORSION_TYPE_CONFLICT_ANALYSIS.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    print(json.dumps({
        "status": report["status"],
        "signature_count": report["signature_count"],
        "nonlocal_collision_count": report["nonlocal_collision_count"],
        "shared_signature_count": report["shared_signature_count"],
        "sha256": sha256(path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
