#!/usr/bin/env python3
"""Map an O6U rotor central bond to CHARMM atom types and parameter lines."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_psf_text(text: str) -> tuple[dict[int, dict], list[tuple[int, int, int, int]]]:
    lines = text.splitlines()
    atom_marker = next(index for index, line in enumerate(lines) if "!NATOM" in line)
    atom_count = int(lines[atom_marker].split()[0])
    atoms: dict[int, dict] = {}
    for line in lines[atom_marker + 1 : atom_marker + 1 + atom_count]:
        fields = line.split()
        atom_id = int(fields[0])
        atoms[atom_id] = {
            "atom_id": atom_id,
            "segment": fields[1],
            "residue_id": fields[2],
            "residue_name": fields[3],
            "atom_name": fields[4],
            "atom_type": fields[5],
            "charge": float(fields[6]),
            "mass": float(fields[7]),
        }
    dihedral_marker = next(index for index, line in enumerate(lines) if "!NPHI" in line)
    dihedral_count = int(lines[dihedral_marker].split()[0])
    values: list[int] = []
    for line in lines[dihedral_marker + 1 :]:
        if "!" in line and values:
            break
        for token in line.split():
            try:
                values.append(int(token))
            except ValueError:
                break
        if len(values) >= 4 * dihedral_count:
            break
    if len(values) < 4 * dihedral_count:
        raise ValueError("PSF dihedral section is truncated")
    dihedrals = [tuple(values[index : index + 4]) for index in range(0, 4 * dihedral_count, 4)]
    return atoms, dihedrals


def map_rotor(
    atoms: dict[int, dict],
    dihedrals: list[tuple[int, int, int, int]],
    central_a: str,
    central_b: str,
) -> list[dict]:
    ids_a = [atom_id for atom_id, atom in atoms.items() if atom["atom_name"] == central_a]
    ids_b = [atom_id for atom_id, atom in atoms.items() if atom["atom_name"] == central_b]
    if len(ids_a) != 1 or len(ids_b) != 1:
        raise ValueError(f"Central atom names are not unique: {central_a}={ids_a}, {central_b}={ids_b}")
    pair = (ids_a[0], ids_b[0])
    mapped = []
    for ids in dihedrals:
        if (ids[1], ids[2]) not in {pair, pair[::-1]}:
            continue
        mapped.append({
            "atom_ids": list(ids),
            "atom_names": [atoms[value]["atom_name"] for value in ids],
            "atom_types": [atoms[value]["atom_type"] for value in ids],
        })
    if not mapped:
        raise ValueError(f"No PSF dihedrals found around central bond {central_a}-{central_b}")
    return mapped


def parameter_lines(path: Path) -> list[str]:
    lines = path.read_text(errors="replace").splitlines()
    in_dihedrals = False
    result = []
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("DIHEDRALS") or upper == "DIHEDRAL":
            in_dihedrals = True
            continue
        if in_dihedrals and upper.startswith(("IMPROPER", "NONBONDED", "CMAP", "HBOND", "END")):
            break
        if in_dihedrals and stripped and not stripped.startswith(("!", "*")) and len(stripped.split()) >= 7:
            result.append(stripped)
    return result


def match_parameter_lines(atom_types: list[str], lines: list[str]) -> list[str]:
    forward = tuple(atom_types)
    reverse = tuple(reversed(atom_types))
    matches = []
    for line in lines:
        pattern = tuple(line.split()[:4])
        for target in (forward, reverse):
            if all(expected == "X" or expected == observed for expected, observed in zip(pattern, target)):
                matches.append(line)
                break
    return matches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psf", required=True)
    parser.add_argument("--base-prm", required=True)
    parser.add_argument("--ligand-prm", required=True)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--central-a", required=True)
    parser.add_argument("--central-b", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    inputs = {name: Path(value).resolve() for name, value in {
        "psf": args.psf,
        "base_prm": args.base_prm,
        "ligand_prm": args.ligand_prm,
    }.items()}
    atoms, dihedrals = parse_psf_text(inputs["psf"].read_text(errors="replace"))
    mapped = map_rotor(atoms, dihedrals, args.central_a, args.central_b)
    base_lines = parameter_lines(inputs["base_prm"])
    ligand_lines = parameter_lines(inputs["ligand_prm"])
    for row in mapped:
        row["matching_parameter_lines"] = {
            "base": match_parameter_lines(row["atom_types"], base_lines),
            "ligand": match_parameter_lines(row["atom_types"], ligand_lines),
        }
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_rotor_to_charmm_dihedral_mapping",
        "status": "pass_rotor_to_charmm_dihedral_mapping",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rotor_id": args.rotor_id,
        "central_bond_atom_names": [args.central_a, args.central_b],
        "inputs": {name: {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size} for name, path in inputs.items()},
        "incident_dihedral_count": len(mapped),
        "incident_dihedrals": mapped,
        "parameter_mutation": False,
        "production_md_approved": False,
        "interpretation_boundary": (
            "Topology and existing-parameter mapping only. This report does not select, fit, "
            "or promote new bonded/torsion parameters."
        ),
    }
    output_dir.mkdir(parents=True)
    report_path = output_dir / "O6U_ROTOR_TO_CHARMM_DIHEDRAL_MAPPING.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(report_path)
    print(json.dumps({
        "status": report["status"],
        "rotor_id": args.rotor_id,
        "incident_dihedral_count": len(mapped),
        "sha256": sha256(report_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
