#!/usr/bin/env python3
"""Create an auditable coordinate-only 8KCT PDB Reader input.

The only atom-record metadata edit is reassignment of the mature PSEN1 CTF
(author chain B, residues 377-467) to generated chain Z. A TER record is added
between PSEN1 residues 291 and 377. Coordinates, serials, residue numbers,
occupancies, B factors, atom names, elements, covalent LINK/SSBOND records, and
CONECT records remain unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom_record(line: str) -> dict[str, object]:
    return {
        "record": line[:6].strip(),
        "serial": int(line[6:11]),
        "atom_name": line[12:16],
        "altloc": line[16:17],
        "resname": line[17:20].strip(),
        "chain": line[21:22],
        "resid": int(line[22:26]),
        "icode": line[26:27],
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        "occupancy": line[54:60],
        "bfactor": line[60:66],
        "element": line[76:78].strip(),
        "raw": line,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True, type=Path)
    parser.add_argument("--output-pdb", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    args = parser.parse_args()
    lines = args.input_pdb.read_text(encoding="ascii").splitlines()
    source_atoms = [atom_record(line) for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if len(source_atoms) != 10839:
        raise SystemExit(f"FAIL: expected 10839 deposited atoms, observed {len(source_atoms)}")

    retained_prefixes = ("SSBOND", "LINK  ", "ATOM  ", "HETATM", "TER   ", "CONECT", "END   ", "END")
    output = [line.ljust(80) for line in (
        "REMARK 900 GENERATED FOR CHARMM-GUI PDB READER; SOURCE 8KCT ASSEMBLY 1",
        "REMARK 900 PSEN1 NTF AUTH B 76-291; CTF GENERATED CHAIN Z 377-467",
        "REMARK 900 NO RESIDUES MODELLED; DEPOSITED ATOM COORDINATES UNCHANGED",
    )]
    inserted_ter = False
    changed_ctf_atoms = 0
    previous_atom = None
    for line in lines:
        if not line.startswith(retained_prefixes):
            continue
        if line.startswith("ATOM  "):
            record = atom_record(line)
            if record["chain"] == "B" and record["resid"] == 377 and previous_atom and not inserted_ter:
                if previous_atom["chain"] != "B" or previous_atom["resid"] != 291:
                    raise SystemExit("FAIL: PSEN1 377 does not immediately follow PSEN1 291")
                output.append("TER   99999      THR B 291".ljust(80))
                inserted_ter = True
            if record["chain"] == "B" and 377 <= int(record["resid"]) <= 467:
                line = line[:21] + "Z" + line[22:]
                changed_ctf_atoms += 1
            previous_atom = record
        elif line.startswith("TER   ") and len(line) >= 26 and line[21:22] == "B" and line[22:26].strip() == "467":
            line = line[:21] + "Z" + line[22:]
        output.append(line)
    if not inserted_ter:
        raise SystemExit("FAIL: PSEN1 NTF/CTF TER was not inserted")
    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    args.output_pdb.write_text("\n".join(output) + "\n", encoding="ascii", newline="\n")

    generated_lines = args.output_pdb.read_text(encoding="ascii").splitlines()
    generated_atoms = [atom_record(line) for line in generated_lines if line.startswith(("ATOM  ", "HETATM"))]
    if len(generated_atoms) != len(source_atoms):
        raise SystemExit("FAIL: generated atom count changed")
    changed_fields = Counter()
    max_coordinate_change = 0.0
    for source, generated in zip(source_atoms, generated_atoms):
        if source["serial"] != generated["serial"]:
            raise SystemExit("FAIL: atom serial/order changed")
        max_coordinate_change = max(max_coordinate_change, math.dist(source["xyz"], generated["xyz"]))
        for field in ("record", "atom_name", "altloc", "resname", "resid", "icode", "occupancy", "bfactor", "element"):
            if source[field] != generated[field]:
                changed_fields[field] += 1
        if source["chain"] != generated["chain"]:
            changed_fields["chain"] += 1
            if not (source["record"] == "ATOM" and source["chain"] == "B" and generated["chain"] == "Z" and 377 <= int(source["resid"]) <= 467):
                raise SystemExit("FAIL: an unapproved chain-ID edit was introduced")

    component_counts = Counter(atom["resname"] for atom in generated_atoms if atom["record"] == "HETATM")
    component_residue_counts = Counter(
        resname
        for resname, chain, resid, icode in {
            (str(atom["resname"]), str(atom["chain"]), int(atom["resid"]), str(atom["icode"]))
            for atom in generated_atoms
            if atom["record"] == "HETATM"
        }
    )
    chain_residues = {}
    for chain in ("A", "B", "C", "D", "Z"):
        residues = sorted({int(atom["resid"]) for atom in generated_atoms if atom["record"] == "ATOM" and atom["chain"] == chain})
        chain_residues[chain] = {"count": len(residues), "first": min(residues) if residues else None, "last": max(residues) if residues else None}
    c291 = next(atom["xyz"] for atom in generated_atoms if atom["record"] == "ATOM" and atom["chain"] == "B" and atom["resid"] == 291 and atom["atom_name"].strip() == "C")
    n377 = next(atom["xyz"] for atom in generated_atoms if atom["record"] == "ATOM" and atom["chain"] == "Z" and atom["resid"] == 377 and atom["atom_name"].strip() == "N")
    source_atom_chains = {str(atom["chain"]) for atom in source_atoms}
    inserted_ntf_ter = [line for line in generated_lines if line.startswith("TER") and len(line) >= 26 and line[17:20].strip() == "THR" and line[21:22] == "B" and line[22:26].strip() == "291"]
    ctf_terminal_ter = [line for line in generated_lines if line.startswith("TER") and len(line) >= 26 and line[17:20].strip() == "ILE" and line[21:22] == "Z" and line[22:26].strip() == "467"]
    checks = {
        "atom_count_preserved": len(generated_atoms) == 10839,
        "coordinates_exact": max_coordinate_change == 0.0,
        "only_chain_metadata_changed": set(changed_fields) <= {"chain"},
        "psen1_ctf_chain_edit_count_expected": changed_ctf_atoms == changed_fields["chain"] and changed_ctf_atoms > 0,
        "psen1_ntf_range_exact": chain_residues["B"] == {"count": 216, "first": 76, "last": 291},
        "generated_chain_z_unused_in_source": "Z" not in source_atom_chains,
        "psen1_ctf_range_exact": chain_residues["Z"] == {"count": 91, "first": 377, "last": 467},
        "one_inserted_psen1_ntf_ter": len(inserted_ntf_ter) == 1,
        "one_rewritten_psen1_ctf_terminal_ter": len(ctf_terminal_ter) == 1,
        "four_disulfides_retained": sum(line.startswith("SSBOND") for line in generated_lines) == 4,
        "twenty_one_link_records_retained": sum(line.startswith("LINK  ") for line in generated_lines) == 21,
        "one_o6u_retained": component_residue_counts["O6U"] == 1 and component_counts["O6U"] == 35,
        "three_clr_retained": component_residue_counts["CLR"] == 3 and component_counts["CLR"] == 84,
        "two_pc1_retained": component_residue_counts["PC1"] == 2 and component_counts["PC1"] == 78,
        "eighteen_nag_retained": component_residue_counts["NAG"] == 18 and component_counts["NAG"] == 18 * 14,
        "three_bma_retained": component_residue_counts["BMA"] == 3 and component_counts["BMA"] == 3 * 11,
    }
    status = "pass_for_pdb_reader_submission_after_ligand_parameter_approval" if all(checks.values()) else "fail"
    audit = {
        "schema_version": "1.0",
        "status": status,
        "submission_approved": False,
        "release_boundary": "The coordinate derivative is ready, but PDB Reader submission remains blocked until the independently approved O6U parameter record exists.",
        "source": {"path": str(args.input_pdb), "bytes": args.input_pdb.stat().st_size, "sha256": sha256(args.input_pdb)},
        "output": {"path": str(args.output_pdb), "bytes": args.output_pdb.stat().st_size, "sha256": sha256(args.output_pdb)},
        "atom_count": len(generated_atoms),
        "changed_ctf_atom_records": changed_ctf_atoms,
        "changed_fields": dict(changed_fields),
        "maximum_coordinate_change_angstrom": max_coordinate_change,
        "psen1_291_c_to_377_n_distance_angstrom": math.dist(c291, n377),
        "protein_chain_residue_ranges": chain_residues,
        "component_atom_counts": dict(sorted(component_counts.items())),
        "component_residue_counts": dict(sorted(component_residue_counts.items())),
        "checks": checks,
    }
    args.audit_json.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "output_sha256": sha256(args.output_pdb), "audit_sha256": sha256(args.audit_json), "changed_ctf_atoms": changed_ctf_atoms}, sort_keys=True))
    return 0 if status != "fail" else 1


if __name__ == "__main__":
    raise SystemExit(main())
