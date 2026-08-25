#!/usr/bin/env python3
"""Prepare auditable 8KCT native-pose redocking inputs without inventing chemistry."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import gemmi
from rdkit import Chem


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pdb_fields(line: str) -> tuple[str, str, str, str]:
    return line[12:16].strip(), line[17:20].strip(), line[21:22], line[22:26].strip()


def ccd_atom_ids(ccd_path: Path) -> list[str]:
    block = gemmi.cif.read_file(str(ccd_path)).sole_block()
    table = block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"])
    atom_ids: list[str] = []
    for row in table:
        atom_id, element = str(row[0]), str(row[1]).upper()
        if element != "H":
            atom_ids.append(atom_id)
    if len(atom_ids) != 35 or len(set(atom_ids)) != 35:
        raise ValueError(f"O6U CCD must contain 35 unique heavy atoms, observed {len(atom_ids)}")
    return atom_ids


def extract_pdb(source: Path, protein_path: Path, native_ligand_pdb: Path) -> dict[str, object]:
    protein_lines: list[str] = []
    ligand_lines: list[str] = []
    chains: list[str] = []
    last_chain: str | None = None
    ligand_ids: set[tuple[str, str]] = set()
    for raw in source.read_text(encoding="ascii", errors="strict").splitlines():
        record = raw[:6].strip()
        if record == "ATOM":
            _, _, chain, _ = pdb_fields(raw)
            if last_chain is not None and chain != last_chain:
                protein_lines.append("TER")
            protein_lines.append(raw)
            if chain not in chains:
                chains.append(chain)
            last_chain = chain
        elif record == "HETATM":
            atom_name, resname, chain, resid = pdb_fields(raw)
            if resname == "O6U":
                ligand_lines.append(raw)
                ligand_ids.add((chain, resid))
    if protein_lines:
        protein_lines.extend(["TER", "END"])
    if len(ligand_ids) != 1:
        raise ValueError(f"Expected exactly one O6U residue, observed {sorted(ligand_ids)}")
    if len(ligand_lines) != 35:
        raise ValueError(f"Expected 35 deposited O6U heavy atoms, observed {len(ligand_lines)}")
    protein_path.write_text("\n".join(protein_lines) + "\n", encoding="ascii", newline="\n")
    native_ligand_pdb.write_text("\n".join(ligand_lines + ["TER", "END"]) + "\n", encoding="ascii", newline="\n")
    return {
        "protein_atom_count": len([line for line in protein_lines if line.startswith("ATOM")]),
        "protein_chains": chains,
        "o6u_heavy_atom_count": len(ligand_lines),
        "o6u_residue": {"chain": next(iter(ligand_ids))[0], "resid": next(iter(ligand_ids))[1]},
    }


def native_reference_sdf(ideal_sdf: Path, ccd_path: Path, native_pdb: Path, output_sdf: Path) -> dict[str, object]:
    supplier = Chem.SDMolSupplier(str(ideal_sdf), removeHs=False, sanitize=True)
    source = supplier[0] if len(supplier) else None
    if source is None:
        raise ValueError("Could not read O6U ideal SDF")
    if source.GetNumAtoms() != 76 or source.GetNumHeavyAtoms() != 35 or Chem.GetFormalCharge(source) != 0:
        raise ValueError("O6U ideal SDF identity differs from 76 atoms / 35 heavy atoms / formal charge zero")
    centers = Chem.FindMolChiralCenters(source, includeUnassigned=True, includeCIP=True)
    if sorted(label for _, label in centers) != ["S", "S"]:
        raise ValueError(f"O6U stereochemistry differs from two assigned S centres: {centers}")

    atom_ids = ccd_atom_ids(ccd_path)
    native_coordinates: dict[str, tuple[float, float, float, str]] = {}
    for line in native_pdb.read_text(encoding="ascii").splitlines():
        if not line.startswith("HETATM"):
            continue
        atom_name = line[12:16].strip()
        native_coordinates[atom_name] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]), line[76:78].strip().upper()
        )
    if set(native_coordinates) != set(atom_ids):
        raise ValueError("Deposited O6U atom names do not exactly match the CCD heavy-atom set")

    reference = Chem.RemoveHs(source)
    if reference.GetNumAtoms() != len(atom_ids):
        raise ValueError("Hydrogen removal changed the expected O6U heavy-atom count")
    conformer = reference.GetConformer()
    for index, atom_id in enumerate(atom_ids):
        x, y, z, element = native_coordinates[atom_id]
        observed = reference.GetAtomWithIdx(index).GetSymbol().upper()
        if observed != element:
            raise ValueError(f"CCD/SDF atom-order element mismatch at {atom_id}: {observed} != {element}")
        conformer.SetAtomPosition(index, (x, y, z))
        reference.GetAtomWithIdx(index).SetProp("atom_id", atom_id)
    reference.SetProp("_Name", "8KCT_native_O6U_heavy_atom_reference")
    reference.SetProp("provenance", "RCSB 8KCT deposited heavy-atom coordinates plus RCSB O6U CCD bond order/stereochemistry")
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(reference)
    writer.close()
    return {
        "formula": "C27H41F2N5O",
        "formal_charge": 0,
        "hydrogen_complete_atom_count": 76,
        "heavy_atom_count": 35,
        "chiral_centres": centers,
        "isomeric_smiles": Chem.MolToSmiles(Chem.RemoveHs(source), isomericSmiles=True),
        "ccd_heavy_atom_order": atom_ids,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--ccd", type=Path, required=True)
    parser.add_argument("--ideal-sdf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.pdb, args.ccd, args.ideal_sdf):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    protein = args.output_dir / "8KCT_protein_only.pdb"
    native_pdb = args.output_dir / "8KCT_native_O6U.pdb"
    native_sdf = args.output_dir / "8KCT_native_O6U_heavy.sdf"
    structure = extract_pdb(args.pdb, protein, native_pdb)
    ligand = native_reference_sdf(args.ideal_sdf, args.ccd, native_pdb, native_sdf)
    inputs = {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in (args.pdb, args.ccd, args.ideal_sdf)}
    outputs = {path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in (protein, native_pdb, native_sdf)}
    manifest = {
        "schema_version": "1.0",
        "purpose": "native-pose self-redocking protocol QA only; not affinity inference",
        "structure": structure,
        "ligand_identity": ligand,
        "inputs": inputs,
        "outputs": outputs,
    }
    manifest_path = args.output_dir / "native_redocking_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "manifest": str(manifest_path), "manifest_sha256": sha256(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
