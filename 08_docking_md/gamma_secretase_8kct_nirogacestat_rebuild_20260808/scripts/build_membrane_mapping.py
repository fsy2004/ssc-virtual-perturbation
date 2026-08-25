#!/usr/bin/env python3
"""Build the frozen membrane QC mapping record (config/membrane_qc_mapping.json).

All metric settings and gate values are predeclared from the constructed
geometry (step5_input.pdb) BEFORE any production trajectory is inspected:
  - phosphate slab from POPC/DSPC phosphorus z distribution
  - hydrophobic-core half thickness from POPC acyl-chain terminal carbon z
  - protein XY exclusion from protein heavy-atom XY spread around the membrane axis
  - standard literature values for bandwidth/cluster cutoff (recorded here)
The two external membrane metrics (protein-aware APL, POPC order parameters)
remain status=not_available (hard pre-production NO-GO) until a source-hashed
validated APL@Voro/FATSLiM and gorder route exists; this record does not fake
that validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
STEP5 = ROOT / "analysis_config_work" / "step5_input.pdb"
OUT = ROOT / "config" / "membrane_qc_mapping.json"


def parse_step5(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append({
                "index": len(atoms),
                "segid": line[72:76].strip(),
                "resname": line[17:21].strip(),
                "resid": int(line[22:27]),
                "name": line[12:16].strip(),
                "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
            })
    return atoms


def topology_identity_records(path: Path) -> list[dict[str, Any]]:
    import MDAnalysis as mda

    universe = mda.Universe(str(path))
    records: list[dict[str, Any]] = []
    for atom in universe.atoms:
        record: dict[str, Any] = {
            "index": int(atom.index),
            "name": str(atom.name),
            "resname": str(atom.resname),
            "resid": int(atom.resid),
        }
        for field in ("segid", "chainID"):
            try:
                value = getattr(atom, field)
            except Exception:
                continue
            if value not in (None, ""):
                record[field] = str(value)
        records.append(record)
    return records


def main(trajectory_topology: Path | None = None) -> int:
    atoms = parse_step5(STEP5)
    topology_source = (trajectory_topology or STEP5).resolve()
    trajectory_identity = topology_identity_records(topology_source)
    if len(trajectory_identity) != len(atoms):
        raise ValueError(
            f"Trajectory topology atom count differs from step5 index source: {len(trajectory_identity)} != {len(atoms)}"
        )
    print("atoms:", len(atoms))

    phosphates = [a for a in atoms if a["resname"] in ("POPC", "DSPC") and a["name"] == "P"]
    zs = sorted(a["z"] for a in phosphates)
    mid = (zs[0] + zs[-1]) / 2.0
    lower = [a for a in phosphates if a["z"] < mid]
    upper = [a for a in phosphates if a["z"] >= mid]
    lower_mean = sum(a["z"] for a in lower) / len(lower)
    upper_mean = sum(a["z"] for a in upper) / len(upper)
    print(f"phosphates: {len(phosphates)} (lower {len(lower)}, upper {len(upper)}) "
          f"means {lower_mean:.3f}/{upper_mean:.3f} A")

    def atom_record(a: dict[str, Any]) -> dict[str, Any]:
        return dict(trajectory_identity[a["index"]])

    # hydrophobic core: POPC acyl-chain carbon z distribution relative to membrane centre
    membrane_center_z = (lower_mean + upper_mean) / 2.0
    # acyl carbons = all POPC C* except the ester/backbone carbons (C1, C11, C12, C13, C31)
    tail_names = ("C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C210", "C211",
                  "C212", "C213", "C214", "C215", "C216", "C217", "C218",
                  "C32", "C33", "C34", "C35", "C36", "C37", "C38", "C39", "C310", "C311",
                  "C312", "C313", "C314", "C315", "C316", "C317", "C318")
    tail_z = [a["z"] for a in atoms if a["resname"] == "POPC" and a["name"] in tail_names]
    tail_z.sort()
    p025 = tail_z[int(0.025 * len(tail_z))]
    p975 = tail_z[int(0.975 * len(tail_z)) - 1]
    core_half = max(abs(p025 - membrane_center_z), abs(p975 - membrane_center_z))
    print(f"POPC acyl carbons: n={len(tail_z)} p025={p025:.3f} p975={p975:.3f} "
          f"core half = {core_half:.3f} A")

    # protein XY spread around the membrane axis (centre = protein-heavy-atom XY median)
    protein_heavy = [a for a in atoms if a["segid"].startswith("PRO") and a["name"][0] in "CNOS"]
    xs = sorted(a["x"] for a in protein_heavy)
    ys = sorted(a["y"] for a in protein_heavy)
    cx = xs[len(xs) // 2]
    cy = ys[len(ys) // 2]
    radii = sorted(((a["x"] - cx) ** 2 + (a["y"] - cy) ** 2) ** 0.5 for a in protein_heavy)
    r95 = radii[int(0.95 * len(radii))]
    r99 = radii[int(0.99 * len(radii))]
    print(f"protein heavy atoms: {len(protein_heavy)}, XY median ({cx:.2f},{cy:.2f}), "
          f"r95={r95:.3f} A, r99={r99:.3f} A")

    # protein tilt anchors: TM-core upper/lower residues (from Guo 2025 contacts)
    upper_anchor_residues = [379, 380, 381, 425, 431, 432]
    lower_anchor_residues = [261, 268, 271, 272, 282, 287]
    tilt_upper = [atom_record(a) for a in atoms if a["segid"] in ("PROD", "PROE")
                  and a["name"] == "CA" and a["resid"] in upper_anchor_residues]
    tilt_lower = [atom_record(a) for a in atoms if a["segid"] in ("PROD", "PROE")
                  and a["name"] == "CA" and a["resid"] in lower_anchor_residues]
    print(f"tilt anchors: upper {len(tilt_upper)}, lower {len(tilt_lower)}")

    # water oxygen atoms
    water_oxy = [a for a in atoms if a["resname"] == "TIP3" and a["name"] == "OH2"]
    print(f"water O: {len(water_oxy)}")

    protein_heavy_record = [atom_record(a) for a in protein_heavy]
    water_oxy_record = [atom_record(a) for a in water_oxy]
    upper_leaf_record = [atom_record(a) for a in upper]
    lower_leaf_record = [atom_record(a) for a in lower]

    identity_sha = hashlib.sha256(
        json.dumps(trajectory_identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    ).hexdigest()

    record = {
        "schema_version": "1.0",
        "approval_status": "draft_not_for_execution",
        "system_id": "8kct_nirogacestat_native",
        "trajectory_atom_identity_sha256": identity_sha,
        "trajectory_atom_source": str(STEP5.relative_to(ROOT)),
        "trajectory_topology_source": str(topology_source),
        "source_records": {
            "trajectory_atom_source": {
                "path": str(STEP5.relative_to(ROOT)),
                "sha256": hashlib.sha256(STEP5.read_bytes()).hexdigest(),
            },
            "trajectory_topology": {
                "path": str(topology_source),
                "sha256": hashlib.sha256(topology_source.read_bytes()).hexdigest(),
            },
        },
        "box_geometry_required": "orthorhombic",
        "frozen_atom_groups": {
            "upper_leaflet_phosphate_atoms": upper_leaf_record,
            "lower_leaflet_phosphate_atoms": lower_leaf_record,
            "protein_tilt_upper_anchor_atoms": tilt_upper,
            "protein_tilt_lower_anchor_atoms": tilt_lower,
            "protein_heavy_atoms": protein_heavy_record,
            "water_oxygen_atoms": water_oxy_record,
        },
        "metric_settings": {
            "phosphate_density_bandwidth_nm": 0.05,
            "phosphate_density_grid_nm": [-4.0, 4.0, 161],
            "leaflet_hysteresis_nm": 0.2,
            "hydrophobic_core_half_thickness_nm": round(core_half / 10.0, 4),
            "protein_xy_exclusion_nm": round(r99 / 10.0, 4),
            "water_cluster_cutoff_nm": 0.35,
            "orthorhombic_angle_tolerance_deg": 0.01,
            "predeclaration_source": (
                "Phosphate slab, hydrophobic-core half thickness and protein XY exclusion "
                "were computed from the constructed geometry (step5_input.pdb) before any "
                "production trajectory was inspected. Bandwidth 0.05 nm and water cluster "
                "cutoff 0.35 nm follow standard practice (see references/evidence matrix)."
            ),
            "predeclared_build_geometry": {
                "phosphate_count": len(phosphates),
                "lower_leaflet_count": len(lower),
                "upper_leaflet_count": len(upper),
                "lower_leaflet_mean_z_angstrom": round(lower_mean, 4),
                "upper_leaflet_mean_z_angstrom": round(upper_mean, 4),
                "membrane_center_z_angstrom": round(membrane_center_z, 4),
                "hydrophobic_core_half_thickness_angstrom": round(core_half, 4),
                "protein_heavy_xy_r95_angstrom": round(r95, 4),
                "protein_heavy_xy_r99_angstrom": round(r99, 4),
            },
        },
        "qc_gates": {
            "maximum_cumulative_leaflet_flip_events": 0,
            "water_defect_largest_cluster_threshold": 10,
            "water_defect_persistence_frames": 5,
            "maximum_absolute_scd_adjacent_block_change": None,
            "maximum_absolute_scd_first_last_change": None,
            "predeclaration_note": (
                "Leaflet flip count 0 and water-defect thresholds are frozen before production "
                "review; SCD gates are intentionally null until a validated gorder route and "
                "an independently justified threshold record are versioned (external NO-GO)."
            ),
        },
        "external_metrics": {
            "protein_aware_area_per_lipid": {
                "status": "not_available",
                "reason": ("Hard pre-production NO-GO until this record is replaced by a "
                           "source-hashed validated APL@Voro v3.3 or FATSLiM route and exact "
                           "rep01-rep03 outputs."),
                "activation_requirements": [
                    "version-pinned source archive and version capture",
                    "exact POPC/protein atom mapping and build-level validation",
                    "exact command record and hash-bound rep01-rep03 outputs",
                ],
            },
            "popc_deuterium_order_parameters": {
                "status": "not_available",
                "reason": ("Hard pre-production NO-GO until this record is replaced by a "
                           "source-hashed validated gorder route with exact CHARMM36 POPC "
                           "mappings and fixed-block rep01-rep03 outputs."),
                "activation_requirements": [
                    "version-pinned gorder source archive and version capture",
                    "explicit CHARMM36 POPC carbon-hydrogen mapping and build validation",
                    "independently justified S_CD drift thresholds",
                    "exact five-block command record and hash-bound rep01-rep03 outputs",
                ],
            },
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256(OUT.read_bytes()).hexdigest()
    OUT.with_suffix(".json.sha256").write_text(f"{digest}  {OUT.name}\n", encoding="ascii")
    print(f"WROTE {OUT} sha256={digest}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--trajectory-topology",
        type=Path,
        default=STEP5,
        help="Frozen production TPR for formal records; defaults to step5 PDB for local contract tests",
    )
    args = parser.parse_args()
    sys.exit(main(args.trajectory_topology))
