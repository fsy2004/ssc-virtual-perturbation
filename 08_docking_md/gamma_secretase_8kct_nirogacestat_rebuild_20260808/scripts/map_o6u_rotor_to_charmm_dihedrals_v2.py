#!/usr/bin/env python3
"""Map a CCD-named O6U rotor to CHARMM dihedrals through the bound correspondence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from map_o6u_rotor_to_charmm_dihedrals_v1 import (
    map_rotor,
    match_parameter_lines,
    parameter_lines,
    parse_psf_text,
    sha256,
)


def parse_correspondence_text(text: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        ccd = row["ccd_atom_id"]
        cgenff = row["cgenff_atom_name"]
        if ccd in mapping:
            raise ValueError(f"Duplicate CCD atom in correspondence: {ccd}")
        mapping[ccd] = cgenff
    if not mapping:
        raise ValueError("Empty atom correspondence")
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psf", required=True)
    parser.add_argument("--base-prm", required=True)
    parser.add_argument("--ligand-prm", required=True)
    parser.add_argument("--correspondence", required=True)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--ccd-central-a", required=True)
    parser.add_argument("--ccd-central-b", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    inputs = {name: Path(value).resolve() for name, value in {
        "psf": args.psf,
        "base_prm": args.base_prm,
        "ligand_prm": args.ligand_prm,
        "correspondence": args.correspondence,
    }.items()}
    correspondence = parse_correspondence_text(inputs["correspondence"].read_text())
    missing = [name for name in (args.ccd_central_a, args.ccd_central_b) if name not in correspondence]
    if missing:
        raise ValueError(f"CCD central atoms missing from correspondence: {missing}")
    cgenff_a = correspondence[args.ccd_central_a]
    cgenff_b = correspondence[args.ccd_central_b]
    atoms, dihedrals = parse_psf_text(inputs["psf"].read_text(errors="replace"))
    mapped = map_rotor(atoms, dihedrals, cgenff_a, cgenff_b)
    base_lines = parameter_lines(inputs["base_prm"])
    ligand_lines = parameter_lines(inputs["ligand_prm"])
    for row in mapped:
        row["matching_parameter_lines"] = {
            "base": match_parameter_lines(row["atom_types"], base_lines),
            "ligand": match_parameter_lines(row["atom_types"], ligand_lines),
        }
    report = {
        "schema_version": "2.0",
        "report_type": "o6u_rotor_to_charmm_dihedral_mapping",
        "status": "pass_ccd_correspondence_bound_rotor_to_charmm_dihedral_mapping",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rotor_id": args.rotor_id,
        "central_bond": {
            "ccd_atom_names": [args.ccd_central_a, args.ccd_central_b],
            "cgenff_atom_names": [cgenff_a, cgenff_b],
        },
        "inputs": {name: {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size} for name, path in inputs.items()},
        "incident_dihedral_count": len(mapped),
        "incident_dihedrals": mapped,
        "parameter_mutation": False,
        "production_md_approved": False,
        "interpretation_boundary": (
            "CCD-to-CGenFF correspondence-bound topology mapping only. This report does not "
            "select, fit, or promote new bonded/torsion parameters."
        ),
    }
    output_dir.mkdir(parents=True)
    report_path = output_dir / "O6U_ROTOR_TO_CHARMM_DIHEDRAL_MAPPING_V2.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(report_path)
    print(json.dumps({
        "status": report["status"],
        "rotor_id": args.rotor_id,
        "cgenff_central_bond": [cgenff_a, cgenff_b],
        "incident_dihedral_count": len(mapped),
        "sha256": sha256(report_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
