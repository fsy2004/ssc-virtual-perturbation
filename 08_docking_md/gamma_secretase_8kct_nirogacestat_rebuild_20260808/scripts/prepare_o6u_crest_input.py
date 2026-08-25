#!/usr/bin/env python3
"""Prepare an identity-locked neutral O6U XYZ input for CREST sampling."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


EXPECTED_FORMULA = "C27H41F2N5O"
EXPECTED_ATOMS = 76
EXPECTED_HEAVY_ATOMS = 35
EXPECTED_BONDS = 78
EXPECTED_CHIRAL_CENTERS = {10: "S", 13: "S"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_single_sdf(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    molecules = [mol for mol in supplier if mol is not None]
    if len(molecules) != 1:
        raise ValueError(f"Expected exactly one valid SDF molecule, found {len(molecules)}")
    mol = molecules[0]
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError("O6U must contain exactly one explicitly three-dimensional conformer")
    return mol


def validate_identity(mol: Chem.Mol) -> dict[str, object]:
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    centers = dict(
        Chem.FindMolChiralCenters(
            mol,
            includeUnassigned=True,
            includeCIP=True,
            useLegacyImplementation=False,
        )
    )
    observed = {
        "formula": formula,
        "atom_count": mol.GetNumAtoms(),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "bond_count": mol.GetNumBonds(),
        "formal_charge_e": formal_charge,
        "chiral_centers_zero_based": {str(index): label for index, label in sorted(centers.items())},
    }
    failures = []
    if formula != EXPECTED_FORMULA:
        failures.append(f"formula {formula} != {EXPECTED_FORMULA}")
    if mol.GetNumAtoms() != EXPECTED_ATOMS:
        failures.append(f"atom count {mol.GetNumAtoms()} != {EXPECTED_ATOMS}")
    if mol.GetNumHeavyAtoms() != EXPECTED_HEAVY_ATOMS:
        failures.append(f"heavy-atom count {mol.GetNumHeavyAtoms()} != {EXPECTED_HEAVY_ATOMS}")
    if mol.GetNumBonds() != EXPECTED_BONDS:
        failures.append(f"bond count {mol.GetNumBonds()} != {EXPECTED_BONDS}")
    if formal_charge != 0:
        failures.append(f"formal charge {formal_charge} != 0")
    if centers != EXPECTED_CHIRAL_CENTERS:
        failures.append(f"chiral centers {centers} != {EXPECTED_CHIRAL_CENTERS}")
    if failures:
        raise ValueError("; ".join(failures))
    return observed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdf", required=True, type=Path)
    parser.add_argument("--output-xyz", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    source = args.sdf.resolve()
    output = args.output_xyz.resolve()
    report_path = args.report.resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise SystemExit(f"Missing or empty source SDF: {source}")
    if output.exists() or report_path.exists():
        raise SystemExit("Refusing to overwrite an existing XYZ or report")

    mol = load_single_sdf(source)
    identity = validate_identity(mol)
    conformer = mol.GetConformer()
    rows = []
    for atom in mol.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        coordinates = (float(position.x), float(position.y), float(position.z))
        if not all(math.isfinite(value) for value in coordinates):
            raise SystemExit(f"Non-finite coordinate at atom {atom.GetIdx()}")
        rows.append((atom.GetSymbol(), *coordinates))

    output.parent.mkdir(parents=True, exist_ok=True)
    source_hash = sha256(source)
    text = [
        str(len(rows)),
        f"O6U neutral singlet; immutable_source_sha256={source_hash}",
    ]
    text.extend(f"{symbol:<2s} {x: .10f} {y: .10f} {z: .10f}" for symbol, x, y, z in rows)
    output.write_text("\n".join(text) + "\n", encoding="utf-8", newline="\n")

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_crest_input_preparation",
        "status": "pass",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sdf": str(source),
        "source_sdf_sha256": source_hash,
        "output_xyz": str(output),
        "output_xyz_sha256": sha256(output),
        "identity": identity,
        "element_sequence": [row[0] for row in rows],
        "crest_model": {
            "charge_e": 0,
            "unpaired_electrons": 0,
            "role": "conformer_start_generation_only",
            "parameter_fit_target": False,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "xyz_sha256": report["output_xyz_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
