#!/usr/bin/env python3
"""Fail-closed identity audit for the neutral O6U parameterization input.

This script verifies chemical identity and atom correspondence only.  It does
not assign atom types, partial charges, or force-field parameters, and it
cannot approve a ligand model for MD.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tempfile
from collections import Counter
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


EXPECTED_FORMULA = "C27H41F2N5O"
EXPECTED_ELEMENT_COUNTS = {"C": 27, "H": 41, "F": 2, "N": 5, "O": 1}
EXPECTED_FORMAL_CHARGE = 0
EXPECTED_ATOMS = 76
EXPECTED_HEAVY_ATOMS = 35
EXPECTED_BONDS = 78
EXPECTED_CIP = [(10, "S"), (13, "S")]
EXPECTED_COMPONENT = "O6U"
EXPECTED_NATIVE_SITE = {"chain": "B", "residue_number": 502}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def cif_scalar(block: gemmi.cif.Block, tag: str) -> str:
    value = block.find_value(tag)
    require(value is not None and str(value) not in {"", ".", "?"}, f"Missing CCD field {tag}")
    return str(value).strip().strip('"')


def read_ccd(path: Path) -> dict[str, object]:
    block = gemmi.cif.read_file(str(path)).sole_block()
    require(cif_scalar(block, "_chem_comp.id") == EXPECTED_COMPONENT, "CCD component is not O6U")
    formula_spaced = cif_scalar(block, "_chem_comp.formula")
    formula = formula_spaced.replace(" ", "")
    require(formula == EXPECTED_FORMULA, f"Unexpected CCD formula: {formula_spaced}")
    formal_charge = int(cif_scalar(block, "_chem_comp.pdbx_formal_charge"))
    require(formal_charge == EXPECTED_FORMAL_CHARGE, f"Unexpected CCD formal charge: {formal_charge}")
    formula_weight = float(cif_scalar(block, "_chem_comp.formula_weight"))
    require(abs(formula_weight - 489.644) <= 0.001, f"Unexpected CCD formula weight: {formula_weight}")

    atom_tags = [
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.alt_atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_aromatic_flag",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
        "_chem_comp_atom.pdbx_stereo_config",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
        "_chem_comp_atom.pdbx_ordinal",
    ]
    atoms: list[dict[str, object]] = []
    for row in block.find(atom_tags):
        atoms.append(
            {
                "atom_id": str(row[0]),
                "alt_atom_id": str(row[1]),
                "element": str(row[2]).upper(),
                "formal_charge": int(str(row[3])),
                "aromatic": str(row[4]),
                "leaving": str(row[5]),
                "stereo": str(row[6]),
                "ideal_x": float(str(row[7])),
                "ideal_y": float(str(row[8])),
                "ideal_z": float(str(row[9])),
                "ordinal": int(str(row[10])),
            }
        )
    require(len(atoms) == EXPECTED_ATOMS, f"CCD has {len(atoms)} atoms, expected {EXPECTED_ATOMS}")
    require(len({str(atom["atom_id"]) for atom in atoms}) == EXPECTED_ATOMS, "CCD atom IDs are not unique")
    require([int(atom["ordinal"]) for atom in atoms] == list(range(1, EXPECTED_ATOMS + 1)), "CCD atom ordinals are not contiguous")
    observed_elements = Counter(str(atom["element"]) for atom in atoms)
    require(dict(observed_elements) == EXPECTED_ELEMENT_COUNTS, f"CCD element counts differ: {dict(observed_elements)}")
    require(sum(int(atom["formal_charge"]) for atom in atoms) == EXPECTED_FORMAL_CHARGE, "CCD atom charges do not sum to zero")
    require(all(str(atom["leaving"]) == "N" for atom in atoms), "CCD contains a leaving atom")
    require(all(math.isfinite(float(atom[key])) for atom in atoms for key in ("ideal_x", "ideal_y", "ideal_z")), "CCD has non-finite ideal coordinates")

    bond_tags = [
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
        "_chem_comp_bond.pdbx_stereo_config",
    ]
    bonds = [
        {
            "atom_id_1": str(row[0]),
            "atom_id_2": str(row[1]),
            "order": str(row[2]),
            "aromatic": str(row[3]),
            "stereo": str(row[4]),
        }
        for row in block.find(bond_tags)
    ]
    require(len(bonds) == EXPECTED_BONDS, f"CCD has {len(bonds)} bonds, expected {EXPECTED_BONDS}")
    atom_ids = {str(atom["atom_id"]) for atom in atoms}
    require(all(str(bond["atom_id_1"]) in atom_ids and str(bond["atom_id_2"]) in atom_ids for bond in bonds), "CCD bond references an unknown atom")
    require(all(str(bond["order"]) in {"SING", "DOUB", "TRIP"} for bond in bonds), "CCD has an unsupported or undefined bond order")
    return {
        "name": cif_scalar(block, "_chem_comp.name"),
        "formula": formula,
        "formal_charge": formal_charge,
        "formula_weight": formula_weight,
        "atoms": atoms,
        "bonds": bonds,
    }


def read_sdf(path: Path, ccd: dict[str, object]) -> tuple[Chem.Mol, dict[str, object]]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    require(len(supplier) == 1 and supplier[0] is not None, "SDF must contain exactly one readable molecule")
    molecule = supplier[0]
    require(molecule.GetNumAtoms() == EXPECTED_ATOMS, f"SDF has {molecule.GetNumAtoms()} atoms")
    require(molecule.GetNumHeavyAtoms() == EXPECTED_HEAVY_ATOMS, f"SDF has {molecule.GetNumHeavyAtoms()} heavy atoms")
    require(molecule.GetNumBonds() == EXPECTED_BONDS, f"SDF has {molecule.GetNumBonds()} bonds")
    require(Chem.GetFormalCharge(molecule) == EXPECTED_FORMAL_CHARGE, "SDF formal charge is not zero")
    require(len(Chem.GetMolFrags(molecule)) == 1, "SDF contains disconnected components")
    require(all(atom.GetNumRadicalElectrons() == 0 for atom in molecule.GetAtoms()), "SDF contains radical electrons")
    require(rdMolDescriptors.CalcMolFormula(molecule) == EXPECTED_FORMULA, "SDF formula differs from the O6U CCD")
    # Average molecular weights vary slightly with the atomic-weight table used
    # by the CCD and RDKit.  Formula and element counts are the exact identity
    # gates; 0.02 Da safely covers the observed table-version difference.
    require(abs(Descriptors.MolWt(molecule) - float(ccd["formula_weight"])) <= 0.02, "SDF molecular weight differs from the O6U CCD")
    centers = sorted(Chem.FindMolChiralCenters(molecule, includeUnassigned=True, includeCIP=True))
    require(centers == EXPECTED_CIP, f"SDF stereochemistry differs: {centers}")

    ccd_atoms = list(ccd["atoms"])
    sdf_elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    require(sdf_elements == [str(atom["element"]) for atom in ccd_atoms], "SDF atom order/elements differ from CCD ordinals")
    atom_id_to_index = {str(atom["atom_id"]): index for index, atom in enumerate(ccd_atoms)}
    ccd_bonds: dict[tuple[int, int], tuple[str, str]] = {}
    for bond in list(ccd["bonds"]):
        pair = tuple(sorted((atom_id_to_index[str(bond["atom_id_1"])], atom_id_to_index[str(bond["atom_id_2"])])))
        ccd_bonds[pair] = (str(bond["order"]), str(bond["aromatic"]))
    sdf_bonds: dict[tuple[int, int], Chem.Bond] = {}
    for bond in molecule.GetBonds():
        pair = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        sdf_bonds[pair] = bond
    require(set(sdf_bonds) == set(ccd_bonds), "SDF and CCD connectivity differ")
    order_map = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB", Chem.BondType.TRIPLE: "TRIP"}
    for pair, (ccd_order, ccd_aromatic) in ccd_bonds.items():
        sdf_bond = sdf_bonds[pair]
        if ccd_aromatic == "Y":
            require(sdf_bond.GetIsAromatic(), f"SDF loses CCD aromaticity at bond {pair}")
        else:
            require(not sdf_bond.GetIsAromatic(), f"SDF introduces aromaticity at bond {pair}")
            require(order_map.get(sdf_bond.GetBondType()) == ccd_order, f"SDF bond order differs at bond {pair}")
    return molecule, {
        "atom_count": molecule.GetNumAtoms(),
        "heavy_atom_count": molecule.GetNumHeavyAtoms(),
        "hydrogen_count": molecule.GetNumAtoms() - molecule.GetNumHeavyAtoms(),
        "bond_count": molecule.GetNumBonds(),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": Descriptors.MolWt(molecule),
        "assigned_chiral_centres_zero_based": [[index, label] for index, label in centers],
        "canonical_isomeric_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True),
    }


def read_native_structure(path: Path, ccd: dict[str, object]) -> tuple[dict[str, tuple[float, float, float]], dict[str, object]]:
    structure = gemmi.read_structure(str(path))
    require(len(structure) == 1, "8KCT coordinate file must contain exactly one model")
    residues = [(chain.name, residue) for chain in structure[0] for residue in chain if residue.name == EXPECTED_COMPONENT]
    require(len(residues) == 1, f"Expected one O6U residue, observed {len(residues)}")
    chain_name, residue = residues[0]
    require(chain_name == EXPECTED_NATIVE_SITE["chain"] and residue.seqid.num == EXPECTED_NATIVE_SITE["residue_number"], "Native O6U is not 8KCT chain B residue 502")
    require(len(residue) == EXPECTED_HEAVY_ATOMS, f"Native O6U has {len(residue)} atoms")
    native: dict[str, tuple[float, float, float]] = {}
    for atom in residue:
        atom_id = atom.name.strip()
        require(atom_id not in native, f"Duplicate native O6U atom name: {atom_id}")
        require(atom.altloc in {"\x00", " ", "A"}, f"Unresolved alternate location for O6U atom {atom_id}: {atom.altloc!r}")
        require(atom.occ >= 0.99, f"Native O6U atom {atom_id} occupancy is below 0.99")
        xyz = (float(atom.pos.x), float(atom.pos.y), float(atom.pos.z))
        require(all(math.isfinite(value) for value in xyz), f"Native O6U atom {atom_id} has non-finite coordinates")
        native[atom_id] = xyz
    ccd_heavy = {str(atom["atom_id"]): str(atom["element"]) for atom in list(ccd["atoms"]) if str(atom["element"]) != "H"}
    require(set(native) == set(ccd_heavy), "Native O6U atom names do not exactly match the CCD heavy-atom set")
    for atom in residue:
        require(atom.element.name.upper() == ccd_heavy[atom.name.strip()], f"Native/CCD element mismatch for {atom.name.strip()}")
    return native, {
        "pdb_id": "8KCT",
        "model_count": len(structure),
        "component_id": EXPECTED_COMPONENT,
        "chain": chain_name,
        "residue_number": residue.seqid.num,
        "heavy_atom_count": len(residue),
        "minimum_occupancy": min(float(atom.occ) for atom in residue),
        "alternate_locations": sorted({atom.altloc for atom in residue if atom.altloc not in {"\x00", " "}}),
    }


def write_outputs(ccd_path: Path, sdf_path: Path, structure_path: Path, output_dir: Path) -> Path:
    ccd = read_ccd(ccd_path)
    molecule, sdf_record = read_sdf(sdf_path, ccd)
    native, native_record = read_native_structure(structure_path, ccd)
    output_dir.mkdir(parents=True, exist_ok=True)

    # The RCSB SDF contains genuine 3D coordinates but leaves the molfile
    # dimension marker blank, which makes RDKit warn that a nominally 2D
    # molecule has nonzero Z coordinates.  Preserve that file as the immutable
    # source and write a derived, explicitly 3D SDF without changing atom order,
    # coordinates, connectivity, formal charge, or stereochemistry.
    approved_input = output_dir / "O6U_neutral_hydrogen_complete_3D.sdf"
    molecule.GetConformer().Set3D(True)
    molecule.SetProp("_Name", "O6U_neutral_hydrogen_complete_3D")
    molecule.SetProp("SOURCE_COMPONENT", "RCSB_O6U")
    molecule.SetProp("SOURCE_SDF_SHA256", sha256(sdf_path))
    writer3d = Chem.SDWriter(str(approved_input))
    writer3d.write(molecule)
    writer3d.close()
    normalized, normalized_record = read_sdf(approved_input, ccd)
    require(normalized.GetConformer().Is3D(), "Derived parameterization SDF is not explicitly marked 3D")
    source_conf = molecule.GetConformer()
    normalized_conf = normalized.GetConformer()
    maximum_coordinate_change = max(
        abs(source_conf.GetAtomPosition(index)[axis] - normalized_conf.GetAtomPosition(index)[axis])
        for index in range(EXPECTED_ATOMS)
        for axis in range(3)
    )
    require(maximum_coordinate_change <= 1e-6, f"3D SDF normalization changed coordinates by {maximum_coordinate_change} Angstrom")

    table_path = output_dir / "O6U_preparameterization_atom_correspondence.tsv"
    fields = [
        "ccd_ordinal", "ccd_atom_id", "ccd_alt_atom_id", "element", "hydrogen_or_heavy",
        "ccd_formal_charge", "ccd_aromatic", "ccd_stereo", "sdf_zero_based_index",
        "sdf_element", "native_8kct_present", "native_x_angstrom", "native_y_angstrom",
        "native_z_angstrom", "cgenff_atom_name", "cgenff_atom_type", "gromacs_atom_name",
    ]
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, ccd_atom in enumerate(list(ccd["atoms"])):
            atom_id = str(ccd_atom["atom_id"])
            xyz = native.get(atom_id)
            writer.writerow(
                {
                    "ccd_ordinal": ccd_atom["ordinal"],
                    "ccd_atom_id": atom_id,
                    "ccd_alt_atom_id": ccd_atom["alt_atom_id"],
                    "element": ccd_atom["element"],
                    "hydrogen_or_heavy": "hydrogen" if ccd_atom["element"] == "H" else "heavy",
                    "ccd_formal_charge": ccd_atom["formal_charge"],
                    "ccd_aromatic": ccd_atom["aromatic"],
                    "ccd_stereo": ccd_atom["stereo"],
                    "sdf_zero_based_index": index,
                    "sdf_element": molecule.GetAtomWithIdx(index).GetSymbol(),
                    "native_8kct_present": "yes" if xyz else "no_hydrogen_not_resolved",
                    "native_x_angstrom": "" if xyz is None else f"{xyz[0]:.3f}",
                    "native_y_angstrom": "" if xyz is None else f"{xyz[1]:.3f}",
                    "native_z_angstrom": "" if xyz is None else f"{xyz[2]:.3f}",
                    "cgenff_atom_name": "",
                    "cgenff_atom_type": "",
                    "gromacs_atom_name": "",
                }
            )

    audit = {
        "schema_version": "1.0",
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "audit_scope": "chemical identity, explicit hydrogens, stereochemistry, connectivity, and native heavy-atom correspondence only",
        "overall_status": "local_identity_pass_parameterization_blocked",
        "md_parameterization_approved": False,
        "parameterization_blockers": [
            "frozen CGenFF version and complete initial penalty inventory are absent",
            "targeted FFParam/QM validation and final residual-penalty review are absent",
            "final CHARMM atom types and partial charges are absent",
            "CHARMM-to-GROMACS energy-component regression is absent",
            "independent parameter and topology approvals are absent",
        ],
        "prohibited_inferences": [
            "This audit does not validate any force-field parameter or partial charge.",
            "The official SDF is not a production topology.",
            "Blank CGenFF and GROMACS columns may not be filled without retained tool output and review.",
        ],
        "identity": {
            "component_id": EXPECTED_COMPONENT,
            "name": ccd["name"],
            "formula": ccd["formula"],
            "formal_charge": ccd["formal_charge"],
            "formula_weight": ccd["formula_weight"],
        },
        "sdf_validation": sdf_record,
        "native_structure_validation": native_record,
        "inputs": {
            "ccd": {"path": str(ccd_path.resolve()), "bytes": ccd_path.stat().st_size, "sha256": sha256(ccd_path)},
            "hydrogen_complete_sdf": {"path": str(sdf_path.resolve()), "bytes": sdf_path.stat().st_size, "sha256": sha256(sdf_path)},
            "native_structure": {"path": str(structure_path.resolve()), "bytes": structure_path.stat().st_size, "sha256": sha256(structure_path)},
        },
        "outputs": {
            "parameterization_input_sdf": {
                "path": str(approved_input.resolve()),
                "bytes": approved_input.stat().st_size,
                "sha256": sha256(approved_input),
                "explicit_3d": normalized.GetConformer().Is3D(),
                "maximum_coordinate_change_from_official_sdf_angstrom": maximum_coordinate_change,
                "validated_identity": normalized_record,
            },
            "atom_correspondence": {"path": str(table_path.resolve()), "bytes": table_path.stat().st_size, "sha256": sha256(table_path), "rows": EXPECTED_ATOMS},
        },
    }
    audit_path = output_dir / "O6U_PREPARAMETERIZATION_AUDIT.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    return audit_path


def self_test(ccd_path: Path, sdf_path: Path, structure_path: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="o6u_audit_") as tmp:
        root = Path(tmp)
        audit_path = write_outputs(ccd_path, sdf_path, structure_path, root / "valid")
        payload = json.loads(audit_path.read_text(encoding="utf-8"))
        require(payload["overall_status"] == "local_identity_pass_parameterization_blocked", "Valid audit did not retain the parameterization block")
        require(payload["md_parameterization_approved"] is False, "Identity audit incorrectly approved parameters")
        bad_mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)[0]
        # N04 (CCD ordinal 30, zero-based index 29) accepts a formal positive
        # charge without creating an invalid carbon valence, so the negative
        # test reaches the explicit charge gate rather than merely failing SDF
        # sanitization.
        bad_mol.GetAtomWithIdx(29).SetFormalCharge(1)
        bad_path = root / "tampered_charge.sdf"
        writer = Chem.SDWriter(str(bad_path))
        writer.write(bad_mol)
        writer.close()
        rejected = False
        try:
            write_outputs(ccd_path, bad_path, structure_path, root / "invalid")
        except ValueError:
            rejected = True
        require(rejected, "Tampered formal charge was not rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ccd", type=Path, required=True)
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    for path in (args.ccd, args.sdf, args.structure):
        require(path.is_file(), f"Missing input: {path}")
    if args.self_test:
        self_test(args.ccd, args.sdf, args.structure)
        print(json.dumps({"status": "pass", "test": "valid input plus formal-charge tamper rejection"}))
        return 0
    require(args.output_dir is not None, "--output-dir is required unless --self-test is used")
    audit_path = write_outputs(args.ccd, args.sdf, args.structure, args.output_dir)
    print(json.dumps({"status": "pass", "audit": str(audit_path), "sha256": sha256(audit_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
