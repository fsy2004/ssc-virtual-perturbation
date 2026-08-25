#!/usr/bin/env python3
"""Fail-closed audit of the official CHARMM-GUI/CGenFF O6U assignment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from itertools import permutations
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_mol2(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    atoms = []
    bonds = []
    section = None
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            section = line.removeprefix("@<TRIPOS>")
            continue
        if not line.strip():
            continue
        if section == "ATOM":
            fields = line.split()
            atoms.append(
                {
                    "ordinal": int(fields[0]),
                    "name": fields[1],
                    "element": re.match(r"[A-Za-z]+", fields[5]).group(0).split(".")[0],
                    "xyz": tuple(float(x) for x in fields[2:5]),
                }
            )
        elif section == "BOND":
            fields = line.split()
            bonds.append((int(fields[1]), int(fields[2]), fields[3]))
    return atoms, bonds


def parse_pdb(path: Path):
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            element = line[76:78].strip() or re.match(r"[A-Za-z]+", name).group(0)[0]
            atoms.append(
                {
                    "ordinal": int(line[6:11]),
                    "name": name,
                    "element": element,
                    "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                }
            )
    return atoms


def parse_rtf(path: Path):
    atoms = []
    bonds = []
    impropers = []
    summary = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("RESI "):
            match = re.search(
                r"RESI\s+(\S+)\s+([-+]?\d+(?:\.\d+)?)\s+!\s+param penalty=\s*([\d.]+)\s*;\s*charge penalty=\s*([\d.]+)",
                line,
            )
            if not match:
                raise ValueError("Cannot parse RESI summary")
            summary = {
                "residue": match.group(1),
                "declared_charge_e": float(match.group(2)),
                "max_parameter_penalty": float(match.group(3)),
                "max_charge_penalty": float(match.group(4)),
            }
        elif line.startswith("ATOM "):
            match = re.match(r"ATOM\s+(\S+)\s+(\S+)\s+([-+]?\d+(?:\.\d+)?)\s*!\s*([\d.]+)", line)
            if not match:
                raise ValueError(f"Cannot parse ATOM line: {line}")
            atoms.append(
                {
                    "name": match.group(1),
                    "type": match.group(2),
                    "charge_e": float(match.group(3)),
                    "charge_penalty": float(match.group(4)),
                }
            )
        elif line.startswith("BOND "):
            fields = line.split()
            for left, right in zip(fields[1::2], fields[2::2]):
                bonds.append((left, right))
        elif line.startswith("IMPR "):
            impropers.append(tuple(line.split()[1:5]))
    if summary is None:
        raise ValueError("RTF lacks RESI summary")
    return summary, atoms, bonds, impropers


def parse_prm(path: Path):
    section = None
    terms = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if stripped in {"BONDS", "ANGLES", "DIHEDRALS", "IMPROPERS"}:
            section = stripped.lower()
            continue
        if section and stripped and not stripped.startswith(("!", "*", "END")):
            match = re.search(r"penalty=\s*([\d.]+)", line)
            if match:
                ntypes = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}[section]
                fields = stripped.split()
                terms.append(
                    {
                        "section": section,
                        "atom_types": fields[:ntypes],
                        "penalty": float(match.group(1)),
                        "line_number": line_number,
                        "raw_line": line,
                    }
                )
    return terms


def enumerate_paths(names, adjacency, length):
    paths = set()

    def walk(path):
        if len(path) == length:
            rev = tuple(reversed(path))
            paths.add(min(tuple(path), rev))
            return
        for nxt in adjacency[path[-1]]:
            if nxt not in path:
                walk(path + [nxt])

    for name in names:
        walk([name])
    return sorted(paths)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol2", type=Path, required=True)
    parser.add_argument("--correspondence", type=Path, required=True)
    parser.add_argument("--rtf", type=Path, required=True)
    parser.add_argument("--prm", type=Path, required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    mol2_atoms, mol2_bonds = parse_mol2(args.mol2)
    pdb_atoms = parse_pdb(args.pdb)
    summary, rtf_atoms, rtf_bonds, rtf_impropers = parse_rtf(args.rtf)
    prm_terms = parse_prm(args.prm)
    if not (len(mol2_atoms) == len(pdb_atoms) == len(rtf_atoms) == 76):
        raise SystemExit("FAIL: expected 76 atoms in MOL2, PDB, and RTF")

    with args.correspondence.open(encoding="utf-8", newline="") as handle:
        correspondence = list(csv.DictReader(handle, delimiter="\t"))
    if len(correspondence) != 76:
        raise SystemExit("FAIL: expected 76 CCD correspondence rows")

    src_edges = {tuple(sorted((a, b))) for a, b, _ in mol2_bonds}
    src_adjacency = defaultdict(set)
    for left, right in src_edges:
        src_adjacency[left].add(right)
        src_adjacency[right].add(left)
    rtf_adjacency = defaultdict(set)
    for left, right in rtf_bonds:
        rtf_adjacency[left].add(right)
        rtf_adjacency[right].add(left)

    # Heavy atoms are preserved in submitted order by Ligand Reader. Hydrogen
    # names/order can be regenerated. Resolve those within each heavy-atom
    # parent by the minimum coordinate-displacement assignment, rather than
    # making an unsafe global ordinal assumption.
    output_by_name = {rtf["name"]: (rtf, pdb) for rtf, pdb in zip(rtf_atoms, pdb_atoms)}
    name_to_ordinal = {}
    for index in range(35):
        src, pdb, rtf = mol2_atoms[index], pdb_atoms[index], rtf_atoms[index]
        if src["element"].upper() != pdb["element"].upper() or src["element"].upper() == "H":
            raise SystemExit(f"FAIL: heavy-atom order/element mismatch at ordinal {index + 1}")
        name_to_ordinal[rtf["name"]] = src["ordinal"]
    for output_parent, source_parent in list(name_to_ordinal.items()):
        source_h = sorted(x for x in src_adjacency[source_parent] if mol2_atoms[x - 1]["element"].upper() == "H")
        output_h = sorted(x for x in rtf_adjacency[output_parent] if output_by_name[x][1]["element"].upper() == "H")
        if len(source_h) != len(output_h):
            raise SystemExit(f"FAIL: hydrogen count differs at parent {output_parent}")
        if not source_h:
            continue
        best = min(
            permutations(output_h),
            key=lambda ordering: sum(
                math.dist(mol2_atoms[source_ordinal - 1]["xyz"], output_by_name[output_name][1]["xyz"])
                for source_ordinal, output_name in zip(source_h, ordering)
            ),
        )
        for source_ordinal, output_name in zip(source_h, best):
            name_to_ordinal[output_name] = source_ordinal
    if len(name_to_ordinal) != 76 or len(set(name_to_ordinal.values())) != 76:
        raise SystemExit("FAIL: atom mapping is incomplete or non-bijective")

    ordinal_to_name = {ordinal: name for name, ordinal in name_to_ordinal.items()}
    mapped = []
    for index, (src, ccd) in enumerate(zip(mol2_atoms, correspondence), 1):
        if src["ordinal"] != index or ccd["ccd_atom_id"] != src["name"]:
            raise SystemExit(f"FAIL: CCD/MOL2 name or order mismatch at ordinal {index}")
        output_name = ordinal_to_name[index]
        rtf, pdb = output_by_name[output_name]
        if src["element"].upper() != pdb["element"].upper():
            raise SystemExit(f"FAIL: mapped element mismatch at ordinal {index}")
        mapped.append(
            {
                "ordinal": index,
                "ccd_atom_id": ccd["ccd_atom_id"],
                "ccd_alt_atom_id": ccd["ccd_alt_atom_id"],
                "element": ccd["element"],
                "cgenff_atom_name": rtf["name"],
                "cgenff_atom_type": rtf["type"],
                "charge_e": rtf["charge_e"],
                "charge_penalty": rtf["charge_penalty"],
                "ligand_reader_coordinate_shift_angstrom": math.dist(src["xyz"], pdb["xyz"]),
            }
        )

    rtf_edges = {tuple(sorted((name_to_ordinal[a], name_to_ordinal[b]))) for a, b in rtf_bonds}
    if src_edges != rtf_edges:
        raise SystemExit("FAIL: CGenFF RTF connectivity differs from submitted MOL2")

    atom_by_name = {row["cgenff_atom_name"]: row for row in mapped}
    adjacency = defaultdict(set)
    for left, right in rtf_bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    concrete = {2: enumerate_paths(atom_by_name, adjacency, 2), 3: enumerate_paths(atom_by_name, adjacency, 3), 4: enumerate_paths(atom_by_name, adjacency, 4)}
    for term in prm_terms:
        ntypes = len(term["atom_types"])
        matches = []
        for path in concrete[ntypes]:
            types = [atom_by_name[name]["cgenff_atom_type"] for name in path]
            if types == term["atom_types"] or list(reversed(types)) == term["atom_types"]:
                matches.append(
                    {
                        "cgenff_atom_names": list(path),
                        "ccd_atom_ids": [atom_by_name[name]["ccd_atom_id"] for name in path],
                    }
                )
        term["matched_o6u_terms"] = matches

    charge_sum = round(math.fsum(row["charge_e"] for row in mapped), 12)
    high_charge = [row for row in mapped if row["charge_penalty"] > 10.0]
    high_bonded = [term for term in prm_terms if term["penalty"] > 10.0]
    checks = {
        "atom_count_76": len(mapped) == 76,
        "heavy_atom_count_35": sum(row["element"] != "H" for row in mapped) == 35,
        "hydrogen_count_41": sum(row["element"] == "H" for row in mapped) == 41,
        "bond_count_78": len(src_edges) == 78,
        "connectivity_exact": src_edges == rtf_edges,
        "topology_charge_within_0_0001_e": abs(charge_sum) <= 0.0001,
        "partial_charge_vector_not_all_zero": any(abs(row["charge_e"]) > 1e-12 for row in mapped),
        "all_parameter_terms_mapped_to_o6u": all(term["matched_o6u_terms"] for term in prm_terms),
    }
    status = "pass_initial_assignment_only_requires_targeted_qm_validation" if all(checks.values()) else "fail"
    audit = {
        "schema_version": "1.0",
        "status": status,
        "production_approved": False,
        "release_boundary": "Initial CGenFF assignment only; high-penalty charges and bonded terms remain unvalidated.",
        "cgenff_program_version": "4.0",
        "cgenff_topology_parameter_release": "5.0",
        "summary": summary,
        "computed_topology_charge_e": charge_sum,
        "checks": checks,
        "coordinate_mapping": {
            "method": "immutable heavy-atom order plus exact graph consistency; regenerated hydrogens assigned within each mapped heavy-atom parent by minimum coordinate displacement",
            "maximum_ligand_reader_coordinate_shift_angstrom": max(row["ligand_reader_coordinate_shift_angstrom"] for row in mapped),
            "rms_ligand_reader_coordinate_shift_angstrom": math.sqrt(sum(row["ligand_reader_coordinate_shift_angstrom"] ** 2 for row in mapped) / len(mapped)),
            "maximum_heavy_atom_shift_angstrom": max(row["ligand_reader_coordinate_shift_angstrom"] for row in mapped if row["element"] != "H"),
            "maximum_hydrogen_shift_angstrom": max(row["ligand_reader_coordinate_shift_angstrom"] for row in mapped if row["element"] == "H"),
        },
        "high_charge_penalty_atoms_above_10": high_charge,
        "high_bonded_penalty_terms_above_10": high_bonded,
        "all_initial_bonded_penalty_terms": prm_terms,
        "improper_terms": rtf_impropers,
        "inputs": {str(path): {"sha256": sha256(path), "bytes": path.stat().st_size} for path in [args.mol2, args.correspondence, args.rtf, args.prm, args.pdb]},
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fieldnames = list(mapped[0])
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(mapped)
    print(json.dumps({"status": status, "charge_sum_e": charge_sum, "high_charge_atoms": len(high_charge), "high_bonded_parameter_lines": len(high_bonded), "output_json_sha256": sha256(args.output_json), "output_tsv_sha256": sha256(args.output_tsv)}, sort_keys=True))
    return 0 if status != "fail" else 1


if __name__ == "__main__":
    raise SystemExit(main())
