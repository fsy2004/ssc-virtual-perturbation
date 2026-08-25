#!/usr/bin/env python3
"""Generate builds/analysis.ndx with the frozen analysis index groups.

Groups required by make_analysis_trajectories.py:
  System, Protein_O6U, PSEN1_Core, Protein_Heavy, O6U_Heavy
Generated deterministically from step5_input.pdb (trajectory atom order):
  - System        : all atoms
  - Protein_O6U   : protein (PRO*) + O6U atoms
  - PSEN1_Core    : the 13 Guo-2025 contact-residue C-alpha atoms (fit group)
  - Protein_Heavy : protein heavy atoms (C,N,O,S)
  - O6U_Heavy     : O6U non-hydrogen atoms
Writes the ndx file and its SHA-256 sidecar.
"""
from __future__ import annotations

import hashlib
import json
import sys
import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
STEP5 = ROOT / "analysis_config_work" / "step5_input.pdb"
OUT = ROOT / "builds" / "analysis.ndx"

CONTACT_RESIDUES = {77, 261, 268, 271, 272, 282, 287, 379, 380, 381, 425, 431, 432}


def parse_step5(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            fld = line.split()
            atoms.append({
                "index": len(atoms),
                "segid": fld[10] if len(fld) > 10 else "?",
                "resname": fld[3],
                "resid": int(fld[4]),
                "name": fld[2],
            })
    return atoms


def write_group(handle, name: str, indices: list[int]) -> None:
    handle.write(f"[ {name} ]\n")
    for start in range(0, len(indices), 15):
        # Internal atom indices are zero-based; GROMACS .ndx files are one-based.
        handle.write(" ".join(str(i + 1) for i in indices[start:start + 15]) + "\n")
    handle.write("\n")


def name_is_hydrogen(name: str) -> bool:
    return name.strip().upper().lstrip("0123456789").startswith("H")


def _psen1_segment(segid: str) -> bool:
    return segid in {"PROD", "PROE"} or segid.endswith("_PROD") or segid.endswith("_PROE")


def groups_from_pdb(path: Path) -> dict[str, list[int]]:
    atoms = parse_step5(path)
    system = [a["index"] for a in atoms]
    protein = [a["index"] for a in atoms if a["segid"].startswith("PRO")]
    o6u = [a["index"] for a in atoms if a["resname"] == "O6U"]
    return {
        "System": system,
        "Protein_O6U": protein + o6u,
        "PSEN1_Core": [
            a["index"] for a in atoms
            if _psen1_segment(a["segid"]) and a["name"] == "CA" and a["resid"] in CONTACT_RESIDUES
        ],
        "Protein_Heavy": [
            a["index"] for a in atoms
            if a["segid"].startswith("PRO") and not name_is_hydrogen(a["name"])
        ],
        "O6U_Heavy": [
            a["index"] for a in atoms
            if a["resname"] == "O6U" and not name_is_hydrogen(a["name"])
        ],
    }


def complete_fragment_indices(seed_group: Any) -> list[int]:
    return sorted({int(index) for fragment in seed_group.fragments for index in fragment.indices})


def groups_from_topology(path: Path, atom_order_source: Path = STEP5) -> dict[str, list[int]]:
    import MDAnalysis as mda

    universe = mda.Universe(str(path))
    atom_order_groups = groups_from_pdb(atom_order_source)
    if len(universe.atoms) != len(atom_order_groups["System"]):
        raise ValueError(
            f"TPR/PDB atom-count mismatch: {len(universe.atoms)} != {len(atom_order_groups['System'])}"
        )
    protein = universe.select_atoms("protein")
    o6u = universe.select_atoms("resname O6U")
    seed = protein + o6u
    return {
        "System": [int(index) for index in universe.atoms.indices],
        # gmx trjconv -pbc cluster requires every selected molecule to be complete.
        # Fragment closure retains covalently attached glycans while excluding solvent,
        # free lipids and ions that do not belong to protein/O6U fragments.
        "Protein_O6U": complete_fragment_indices(seed),
        # TPR residue numbers are globally renumbered by GROMACS. These groups use
        # the hash-bound CHARMM-GUI PDB atom order, which is identical to the TPR,
        # so PSEN1 residue identities remain the original structural numbering.
        "PSEN1_Core": atom_order_groups["PSEN1_Core"],
        "Protein_Heavy": atom_order_groups["Protein_Heavy"],
        "O6U_Heavy": atom_order_groups["O6U_Heavy"],
    }


def main(trajectory_topology: Path = STEP5, output: Path = OUT) -> int:
    trajectory_topology = trajectory_topology.resolve()
    output = output.resolve()
    if trajectory_topology.suffix.lower() == ".tpr":
        groups = groups_from_topology(trajectory_topology)
        method = "production_tpr_fragment_closure"
    else:
        groups = groups_from_pdb(trajectory_topology)
        method = "pdb_segment_fallback_not_for_formal_cluster_release"
    system = groups["System"]
    protein_o6u = groups["Protein_O6U"]
    psen1_core = groups["PSEN1_Core"]
    protein_heavy = groups["Protein_Heavy"]
    o6u_heavy = groups["O6U_Heavy"]
    print("atoms:", len(system))
    print(f"System {len(system)} | Protein_O6U {len(protein_o6u)} | "
          f"PSEN1_Core {len(psen1_core)} | Protein_Heavy {len(protein_heavy)} | "
          f"O6U_Heavy {len(o6u_heavy)}")
    if len(psen1_core) != len(CONTACT_RESIDUES):
        raise ValueError(f"PSEN1_Core count {len(psen1_core)} != {len(CONTACT_RESIDUES)}")
    if len(o6u_heavy) != 35:
        raise ValueError(f"O6U_Heavy count {len(o6u_heavy)} != 35")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii") as handle:
        write_group(handle, "System", system)
        write_group(handle, "Protein_O6U", protein_o6u)
        write_group(handle, "PSEN1_Core", psen1_core)
        write_group(handle, "Protein_Heavy", protein_heavy)
        write_group(handle, "O6U_Heavy", o6u_heavy)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(".ndx.sha256").write_text(f"{digest}  {output.name}\n", encoding="ascii")
    source_hash = hashlib.sha256(trajectory_topology.read_bytes()).hexdigest()
    provenance = {
        "schema_version": "1.0",
        "report_type": "analysis_index_build",
        "status": "formal" if method == "production_tpr_fragment_closure" else "local_fallback",
        "method": method,
        "trajectory_topology": {"path": str(trajectory_topology), "sha256": source_hash},
        "atom_order_source": {
            "path": str(STEP5.resolve()),
            "sha256": hashlib.sha256(STEP5.read_bytes()).hexdigest(),
        },
        "group_counts": {name: len(indices) for name, indices in groups.items()},
        "gromacs_index_sha256": digest,
        "gromacs_indexing": "one_based",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    provenance_path = output.with_suffix(".ndx.provenance.json")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {output} sha256={digest}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory-topology", type=Path, default=STEP5)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    sys.exit(main(args.trajectory_topology, args.output))
