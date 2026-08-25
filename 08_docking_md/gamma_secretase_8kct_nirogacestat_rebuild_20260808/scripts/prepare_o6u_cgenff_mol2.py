#!/usr/bin/env python3
"""Create and audit the neutral O6U MOL2 submitted for initial CGenFF typing."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path

from openbabel import openbabel as ob
from rdkit import Chem, rdBase


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mapping(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) == 76, "O6U correspondence table must contain exactly 76 atoms")
    require([int(row["sdf_zero_based_index"]) for row in rows] == list(range(76)), "SDF mapping is not contiguous")
    names = [row["ccd_atom_id"].strip() for row in rows]
    require(len(set(names)) == 76 and all(names), "CCD atom names are empty or non-unique")
    return rows


def read_obmol(path: Path, fmt: str) -> ob.OBMol:
    conversion = ob.OBConversion()
    require(conversion.SetInFormat(fmt), f"Open Babel lacks {fmt} input support")
    molecule = ob.OBMol()
    require(conversion.ReadFile(molecule, str(path)), f"Open Babel could not read {path}")
    return molecule


def graph_signature(molecule: ob.OBMol) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int, int, bool], ...]]:
    atoms = tuple((atom.GetAtomicNum(), atom.GetFormalCharge()) for atom in ob.OBMolAtomIter(molecule))
    bonds = []
    for bond in ob.OBMolBondIter(molecule):
        a, b = sorted((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1))
        aromatic = bool(bond.IsAromatic())
        # MOL2 and SDF may choose different but equivalent Kekule forms.  For
        # aromatic bonds, compare aromatic membership rather than the arbitrary
        # alternating single/double assignment; retain exact order elsewhere.
        bonds.append((a, b, 0 if aromatic else int(bond.GetBondOrder()), aromatic))
    return atoms, tuple(sorted(bonds))


def ob_canonical_isomeric_smiles(molecule: ob.OBMol) -> str:
    conversion = ob.OBConversion()
    require(conversion.SetOutFormat("can"), "Open Babel lacks canonical SMILES output")
    value = conversion.WriteString(molecule).strip().split()[0]
    require(value, "Open Babel produced an empty canonical SMILES")
    return value


def mol2_atom_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    in_atoms = False
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if raw.startswith("@<TRIPOS>"):
            if in_atoms:
                break
            continue
        if in_atoms and raw.strip():
            rows.append(raw.split())
    return rows


def rewrite_mol2_atom_metadata(path: Path, atom_names: list[str]) -> None:
    """Apply immutable CCD names after Open Babel assigns generic MOL2 names."""
    lines = path.read_text(encoding="utf-8").splitlines()
    in_atoms = False
    atom_index = 0
    rewritten: list[str] = []
    for raw in lines:
        if raw.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            rewritten.append(raw)
            continue
        if raw.startswith("@<TRIPOS>"):
            if in_atoms:
                in_atoms = False
            rewritten.append(raw)
            continue
        if in_atoms and raw.strip():
            fields = raw.split()
            require(len(fields) >= 9, "Malformed MOL2 atom record")
            require(atom_index < len(atom_names), "MOL2 contains more atoms than the CCD map")
            fields[1] = atom_names[atom_index]
            fields[7] = "O6U"
            rewritten.append(" ".join(fields))
            atom_index += 1
        else:
            rewritten.append(raw)
    require(atom_index == len(atom_names), "MOL2 atom metadata rewrite did not cover all atoms")
    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8", newline="\n")


def prepare(input_sdf: Path, mapping_tsv: Path, output_mol2: Path, audit_json: Path) -> dict[str, object]:
    rows = load_mapping(mapping_tsv)
    source = read_obmol(input_sdf, "sdf")
    require(source.NumAtoms() == 76 and source.NumBonds() == 78, "Unexpected O6U source atom/bond count")
    require(source.GetTotalCharge() == 0, "O6U source is not neutral")
    require(sum(1 for atom in ob.OBMolAtomIter(source) if atom.GetAtomicNum() == 1) == 41, "O6U must contain 41 explicit hydrogens")
    require([atom.GetAtomicNum() for atom in ob.OBMolAtomIter(source)] == [ob.GetAtomicNum(row["element"]) for row in rows], "CCD/SDF element order differs")

    source.SetTitle("O6U")
    residue = source.NewResidue()
    residue.SetName("O6U")
    residue.SetNum(1)
    for atom, row in zip(ob.OBMolAtomIter(source), rows, strict=True):
        residue.AddAtom(atom)
        residue.SetAtomID(atom, row["ccd_atom_id"])

    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".mol2", delete=False, dir=output_mol2.parent) as temporary:
        temporary_path = Path(temporary.name)
    try:
        conversion = ob.OBConversion()
        require(conversion.SetOutFormat("mol2"), "Open Babel lacks MOL2 output support")
        require(conversion.WriteFile(source, str(temporary_path)), "Open Babel could not write MOL2")
        conversion.CloseOutFile()
        rewrite_mol2_atom_metadata(temporary_path, [row["ccd_atom_id"] for row in rows])
        os.replace(temporary_path, output_mol2)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()

    generated = read_obmol(output_mol2, "mol2")
    require(generated.NumAtoms() == 76 and generated.NumBonds() == 78, "MOL2 atom/bond counts changed")
    require(generated.GetTotalCharge() == 0, "MOL2 formal charge changed")
    require(graph_signature(source) == graph_signature(generated), "MOL2 graph/formal-charge signature changed")
    max_coordinate_change = max(
        ((a.GetX() - b.GetX()) ** 2 + (a.GetY() - b.GetY()) ** 2 + (a.GetZ() - b.GetZ()) ** 2) ** 0.5
        for a, b in zip(ob.OBMolAtomIter(source), ob.OBMolAtomIter(generated), strict=True)
    )
    require(max_coordinate_change <= 1.0e-6, f"MOL2 coordinates changed by {max_coordinate_change:.9g} A")

    atom_rows = mol2_atom_rows(output_mol2)
    require(len(atom_rows) == 76, "Could not parse 76 MOL2 atom records")
    mol2_names = [row[1] for row in atom_rows]
    expected_names = [row["ccd_atom_id"] for row in rows]
    require(mol2_names == expected_names, "MOL2 atom names do not exactly match CCD names/order")
    partial_charges = [float(row[8]) for row in atom_rows]

    sdf_rdkit = Chem.MolFromMolFile(str(input_sdf), sanitize=True, removeHs=False)
    require(sdf_rdkit is not None, "RDKit sanitization of the immutable SDF failed")
    sdf_smiles = Chem.MolToSmiles(Chem.RemoveHs(sdf_rdkit), canonical=True, isomericSmiles=True)
    ob_source_smiles = ob_canonical_isomeric_smiles(source)
    ob_mol2_smiles = ob_canonical_isomeric_smiles(generated)
    require(ob_source_smiles == ob_mol2_smiles, "Open Babel canonical isomeric SMILES changed during MOL2 conversion")

    payload: dict[str, object] = {
        "schema_version": 1,
        "status": "pass_for_initial_cgenff_submission_only",
        "input_sdf": {"path": str(input_sdf.resolve()), "sha256": sha256(input_sdf)},
        "mapping_tsv": {"path": str(mapping_tsv.resolve()), "sha256": sha256(mapping_tsv)},
        "output_mol2": {"path": str(output_mol2.resolve()), "sha256": sha256(output_mol2)},
        "software": {"openbabel": ob.OBReleaseVersion(), "rdkit": rdBase.rdkitVersion},
        "identity": {
            "residue_name": "O6U",
            "atom_count": 76,
            "heavy_atom_count": 35,
            "explicit_hydrogen_count": 41,
            "bond_count": 78,
            "formal_charge": 0,
            "ccd_atom_names_and_order_exact": True,
            "graph_and_formal_charges_exact": True,
            "rdkit_source_canonical_isomeric_smiles": sdf_smiles,
            "openbabel_source_and_mol2_canonical_isomeric_smiles": ob_source_smiles,
            "max_coordinate_change_angstrom": max_coordinate_change,
        },
        "mol2_partial_charge_field": {
            "sum": sum(partial_charges),
            "interpretation": "transport-format values only; CGenFF must assign and report the force-field charges",
        },
        "release_boundary": "This MOL2 audit does not approve CGenFF penalties, QM targets, fitted parameters, or MD use.",
    }
    audit_json.parent.mkdir(parents=True, exist_ok=True)
    audit_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-sdf", type=Path, required=True)
    parser.add_argument("--mapping-tsv", type=Path, required=True)
    parser.add_argument("--output-mol2", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.input_sdf, args.mapping_tsv):
        require(path.is_file(), f"Missing input: {path}")
    payload = prepare(args.input_sdf, args.mapping_tsv, args.output_mol2, args.audit_json)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
