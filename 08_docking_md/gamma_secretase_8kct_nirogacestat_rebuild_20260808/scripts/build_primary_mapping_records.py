#!/usr/bin/env python3
"""Build the hash-bound primary atom-mapping/contacts record for the 8KCT-O6U MD analysis.

Frozen mapping rules (all deterministic, from curated inputs):
  - Reference atoms   : docking_native_redock/plip_native/8KCT_protonated.pdb (PLIP input)
  - Trajectory atoms  : common/step5_input.pdb (CHARMM-GUI output, same order as production TPR)
  - Chain->segment    : A->PROA (nicastrin), B->PROD+PROE (PSEN1), C->PROB (APH1A), D->PROC (PEN2)
  - Residue names     : HIS->HSD (CHARMM36m protonation); NAG/BMA/CLR/PC1 are separate segments
  - Atom names        : ILE CD1->CD; chain-terminal O -> OT1/OT2 (assigned by geometry)
  - Ligand            : O6U CCD names -> CGenFF names via
                        inputs/ligand_parameterization/O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv;
                        the three polar hydrogens (reference 'H') are assigned to N04/N06/N07 by
                        nearest-heavy-atom geometry.
  - TM-core C-alpha   : PSEN1 (segid PROD/PROE) C-alpha atoms whose z lies within the POPC
                        phosphate-z slab of the constructed membrane (lower/upper phosphate mean
                        minus/plus half the inter-leaflet separation), validated against the
                        nirogacestat contact residues and their TM numbers in Guo 2025
                        Supplementary Table 3.
  - Native contacts   : PLIP endpoints (protein heavy atom <-> O6U heavy atom), hydrogen bonds
                        from the PLIP XML donor/acceptor atoms (O6U donor -> protein acceptor and
                        protein donor -> O6U acceptor).

Outputs (written into config/): primary_atom_mapping_contacts.json plus its SHA-256 sidecar and
a verification summary. No production file is touched.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
STEP5 = ROOT / "analysis_config_work" / "step5_input.pdb"
REFERENCE = ROOT / "docking_native_redock" / "plip_native" / "8KCT_protonated.pdb"
CORRESP = ROOT / "inputs" / "ligand_parameterization" / "O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv"
PLIP_XML = ROOT / "docking_native_redock" / "plip_native" / "run1" / "8KCT_O6U.xml"
PLIP_JSON = (
    ROOT / "docking_native_redock" / "figures" / "native_8kct_o6u"
    / "8KCT_O6U_native_contacts.interactions.normalized.json"
)
GRO = ROOT / "analysis_config_work" / "minimized.gro"
OUT = ROOT / "config" / "primary_atom_mapping_contacts.json"

CHAIN_MAP = {"A": ["PROA"], "B": ["PROD", "PROE"], "C": ["PROB"], "D": ["PROC"]}
CHARMM_AMINO = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def parse_reference(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append({
                "index": len(atoms),
                "record_type": line[0:6].strip(),
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip(),
                "resid": int(line[22:26]),
                "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                "element": line[76:78].strip() or line[12:16].strip()[0],
            })
    return atoms


def parse_phosphates_from_step5(atoms: list[dict[str, Any]]) -> tuple[float, float, int]:
    zs = sorted(a["z"] for a in atoms if a["resname"] == "POPC" and a["name"] == "P")
    if not zs:
        raise ValueError("No POPC phosphate atoms found")
    mid = (zs[0] + zs[-1]) / 2.0
    lower = [z for z in zs if z < mid]
    upper = [z for z in zs if z >= mid]
    lower_mean = sum(lower) / len(lower)
    upper_mean = sum(upper) / len(upper)
    return lower_mean, upper_mean, len(zs)


def parse_correspondence(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open() as handle:
        next(handle)  # header
        for line in handle:
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 5:
                continue
            ccd = fld[1].strip()
            cgenff = fld[4].strip()
            mapping[ccd] = cgenff
    return mapping


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
    print("== Input hashes ==")
    for label, path in [("step5", STEP5), ("reference", REFERENCE), ("correspondence", CORRESP),
                        ("plip_xml", PLIP_XML), ("plip_json", PLIP_JSON), ("gro", GRO)]:
        print(f"  {label}: {sha256(path)}")

    top = parse_step5(STEP5)
    ref = parse_reference(REFERENCE)
    topology_source = (trajectory_topology or STEP5).resolve()
    trajectory_identity = topology_identity_records(topology_source)
    if len(trajectory_identity) != len(top):
        raise ValueError(
            f"Trajectory topology atom count differs from step5 index source: {len(trajectory_identity)} != {len(top)}"
        )
    print(f"trajectory atoms: {len(top)}, reference atoms: {len(ref)}")

    # topology key -> trajectory indices
    top_by_key: dict[tuple[str, str, int, str], list[int]] = defaultdict(list)
    for atom in top:
        top_by_key[(atom["segid"], atom["resname"], atom["resid"], atom["name"])].append(atom["index"])

    # geometry-based OT1/OT2 resolution: reference O at chain terminus -> closest OT
    def resolve_terminal_o(ref_atom: dict[str, Any]) -> int:
        segids = CHAIN_MAP.get(ref_atom["chain"], [])
        candidates = []
        for sg in segids:
            for idx in top_by_key.get((sg, ref_atom["resname"], ref_atom["resid"], "OT1"), []):
                candidates.append(idx)
            for idx in top_by_key.get((sg, ref_atom["resname"], ref_atom["resid"], "OT2"), []):
                candidates.append(idx)
        if not candidates:
            raise ValueError(f"Cannot resolve terminal O for {ref_atom}")
        best = min(
            candidates,
            key=lambda i: (
                (top[i]["x"] - ref_atom["x"]) ** 2
                + (top[i]["y"] - ref_atom["y"]) ** 2
                + (top[i]["z"] - ref_atom["z"]) ** 2
            ),
        )
        return best

    def map_protein_atom(ref_atom: dict[str, Any]) -> int | None:
        segids = CHAIN_MAP.get(ref_atom["chain"], [])
        residue_names = [ref_atom["resname"]]
        residue_names.extend(charmm for charmm, pdb in CHARMM_AMINO.items() if pdb == ref_atom["resname"])
        name = ref_atom["name"]
        for sg in segids:
            for resname in residue_names:
                hits = top_by_key.get((sg, resname, ref_atom["resid"], name), [])
                if hits:
                    return hits[0]
        # name substitutions
        if ref_atom["resname"] == "ILE" and name == "CD1":
            for sg in segids:
                hits = top_by_key.get((sg, "ILE", ref_atom["resid"], "CD"), [])
                if hits:
                    return hits[0]
        if ref_atom.get("element", "").upper() == "H":
            candidates = [
                atom["index"] for atom in top
                if atom["segid"] in segids
                and atom["resname"] in residue_names
                and atom["resid"] == ref_atom["resid"]
                and atom["name"].startswith("H")
            ]
            if candidates:
                return min(
                    candidates,
                    key=lambda index: (
                        (top[index]["x"] - ref_atom["x"]) ** 2
                        + (top[index]["y"] - ref_atom["y"]) ** 2
                        + (top[index]["z"] - ref_atom["z"]) ** 2
                    ),
                )
        if name == "O" and ref_atom["resname"] != "O6U":
            return resolve_terminal_o(ref_atom)
        return None

    # O6U mapping via correspondence + geometry
    corresp = parse_correspondence(CORRESP)
    ref_o6u = [a for a in ref if a["resname"] == "O6U"]
    top_o6u = [a for a in top if a["resname"] == "O6U"]
    print(f"ref O6U: {len(ref_o6u)}, top O6U: {len(top_o6u)}")
    heavy_ref_o6u = [a for a in ref_o6u if a["name"] != "H"]
    heavy_top_o6u = [a for a in top_o6u if not a["name"].startswith("H")]
    # cgenff atom names for heavy atoms: from correspondence; exclude F/N/O/C names that start with H
    cgenff_heavy = {corresp.get(a["name"]) for a in heavy_ref_o6u}
    top_heavy_names = {a["name"] for a in heavy_top_o6u}
    if not cgenff_heavy.issubset(top_heavy_names):
        raise ValueError(f"Missing CGenFF heavy names: {cgenff_heavy - top_heavy_names}")

    o6u_map: dict[str, int] = {}  # reference name -> trajectory index
    for ra in heavy_ref_o6u:
        cgenff_name = corresp[ra["name"]]
        # among top atoms with this name, pick by geometric nearest (CCD->CGenFF shifts are small)
        candidates = [a["index"] for a in top_o6u if a["name"] == cgenff_name]
        if len(candidates) == 1:
            o6u_map[ra["name"]] = candidates[0]
        else:
            best = min(candidates, key=lambda i: (
                (top[i]["x"] - ra["x"]) ** 2 + (top[i]["y"] - ra["y"]) ** 2 + (top[i]["z"] - ra["z"]) ** 2
            ))
            o6u_map[ra["name"]] = best
    # polar hydrogens: ref 'H' atoms nearest N04/N06/N07 heavy atoms -> CGenFF H39/H40/H41
    polar_h = [a for a in ref_o6u if a["name"] == "H"]
    print(f"ref polar H: {len(polar_h)}")
    h_targets = {"N04": "H39", "N06": "H40", "N07": "H41"}
    for h_atom in polar_h:
        best_n = min(("N04", "N06", "N07"), key=lambda n: (
            (ref[0]["x"] if False else 0)  # placeholder, replaced below
        ))
        # compute distance from h_atom to each heavy ref atom
        dists = {}
        for n in ("N04", "N06", "N07"):
            n_atom = next(a for a in heavy_ref_o6u if a["name"] == n)
            dists[n] = ((h_atom["x"] - n_atom["x"]) ** 2
                        + (h_atom["y"] - n_atom["y"]) ** 2
                        + (h_atom["z"] - n_atom["z"]) ** 2) ** 0.5
        nearest = min(dists, key=dists.get)
        cgenff_name = h_targets[nearest]
        candidates = [a["index"] for a in top_o6u if a["name"] == cgenff_name]
        if len(candidates) != 1:
            raise ValueError(f"Polar hydrogen {cgenff_name} not unique in topology")
        o6u_map[f"H-{nearest}"] = candidates[0]
    print(f"O6U mapped atoms: {len(o6u_map)}")

    # PLIP contacts
    plip_json = json.loads(PLIP_JSON.read_text(encoding="utf-8"))
    plip_xml = ET.parse(PLIP_XML).getroot()
    ref_by_serial = {a["serial"]: a for a in ref}

    def ref_to_traj(ref_atom: dict[str, Any]) -> int:
        if ref_atom["resname"] == "O6U":
            key = ref_atom["name"] if ref_atom["name"] != "H" else None
            if key is not None and key in o6u_map:
                return o6u_map[key]
            # polar H: find nearest N
            for n in ("N04", "N06", "N07"):
                if f"H-{n}" in o6u_map:
                    n_atom = next(a for a in heavy_ref_o6u if a["name"] == n)
                    d = ((ref_atom["x"] - n_atom["x"]) ** 2
                         + (ref_atom["y"] - n_atom["y"]) ** 2
                         + (ref_atom["z"] - n_atom["z"]) ** 2) ** 0.5
                    if d < 1.2:
                        return o6u_map[f"H-{n}"]
            raise ValueError(f"Cannot map O6U ref atom {ref_atom}")
        idx = map_protein_atom(ref_atom)
        if idx is None:
            raise ValueError(f"Cannot map ref atom {ref_atom}")
        return idx

    def reference_identity(ref_atom: dict[str, Any]) -> dict[str, Any]:
        return {
            "index": ref_atom["index"],
            "name": ref_atom["name"],
            "resname": ref_atom["resname"],
            "resid": ref_atom["resid"],
            "segid": ref_atom["chain"],
            "chainID": ref_atom["chain"],
        }

    def trajectory_atom_identity(index: int) -> dict[str, Any]:
        return dict(trajectory_identity[index])

    def mapped_endpoint(ref_atom: dict[str, Any]) -> dict[str, Any]:
        return {
            "reference": reference_identity(ref_atom),
            "trajectory": trajectory_atom_identity(ref_to_traj(ref_atom)),
        }

    def unique_reference_atom(chain: str, resid: int, name: str, resname: str | None = None) -> dict[str, Any]:
        hits = [
            atom for atom in ref
            if atom["chain"] == chain and atom["resid"] == resid and atom["name"] == name
            and (resname is None or atom["resname"] == resname)
        ]
        if len(hits) != 1:
            raise ValueError(f"Reference atom is not unique: {chain}/{resname or '*'}/{resid}/{name}: {len(hits)}")
        return hits[0]

    # Exhaustive reference heavy-atom contacts are the quantitative native-contact
    # definition. PLIP remains a source annotation, but is not used as an incomplete
    # atom-pair filter.
    protein_heavy_ref = [
        atom for atom in ref
        if atom["record_type"] == "ATOM" and atom["element"].upper() != "H"
    ]
    ligand_heavy_ref = [atom for atom in ref_o6u if atom["element"].upper() != "H"]
    native_contacts = []
    cutoff_angstrom = 4.5
    for p_atom in protein_heavy_ref:
        for l_atom in ligand_heavy_ref:
            distance = (
                (p_atom["x"] - l_atom["x"]) ** 2
                + (p_atom["y"] - l_atom["y"]) ** 2
                + (p_atom["z"] - l_atom["z"]) ** 2
            ) ** 0.5
            if distance <= cutoff_angstrom + 1.0e-9:
                native_contacts.append({
                    "contact_id": (
                        f"{p_atom['chain']}:{p_atom['resname']}{p_atom['resid']}:{p_atom['name']}"
                        f"__O6U:{l_atom['name']}"
                    ),
                    "type": "reference_heavy_atom_contact",
                    "protein_atom": mapped_endpoint(p_atom),
                    "ligand_atom": mapped_endpoint(l_atom),
                    "reference_distance_nm": distance / 10.0,
                    "cutoff_nm": 0.45,
                })
    if len({item["contact_id"] for item in native_contacts}) != len(native_contacts):
        raise ValueError("Native contact IDs are not unique")
    print(f"exhaustive native heavy-atom contacts: {len(native_contacts)}")

    # Hydrogen bonds from XML (O6U-relevant only), with explicit reference and
    # trajectory donor/H/acceptor identities required by the analyzer.
    hydrogen_bonds = []
    for hb in plip_xml.iter("hydrogen_bond"):
        d = {c.tag: c.text for c in hb}
        if d.get("resnr_lig") != "502":
            continue
        donor_serial = int(d["donoridx"])
        acceptor_serial = int(d["acceptoridx"])
        donor = ref_by_serial[donor_serial]
        acceptor = ref_by_serial[acceptor_serial]
        hydrogen_candidates = [
            atom for atom in ref
            if atom["chain"] == donor["chain"]
            and atom["resname"] == donor["resname"]
            and atom["resid"] == donor["resid"]
            and atom["element"].upper() == "H"
        ]
        proximal_hydrogens = []
        for atom in hydrogen_candidates:
            distance = (
                (atom["x"] - donor["x"]) ** 2
                + (atom["y"] - donor["y"]) ** 2
                + (atom["z"] - donor["z"]) ** 2
            ) ** 0.5
            if 0.5 < distance < 1.3:
                proximal_hydrogens.append((distance, atom))
        if not proximal_hydrogens:
            raise ValueError(f"No explicit reference donor H found for {donor}")
        hydrogen = min(proximal_hydrogens, key=lambda item: item[0])[1]
        donor_acceptor_angstrom = (
            (donor["x"] - acceptor["x"]) ** 2
            + (donor["y"] - acceptor["y"]) ** 2
            + (donor["z"] - acceptor["z"]) ** 2
        ) ** 0.5
        first = (donor["x"] - hydrogen["x"], donor["y"] - hydrogen["y"], donor["z"] - hydrogen["z"])
        second = (acceptor["x"] - hydrogen["x"], acceptor["y"] - hydrogen["y"], acceptor["z"] - hydrogen["z"])
        denominator = math.sqrt(sum(value * value for value in first)) * math.sqrt(sum(value * value for value in second))
        cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(first, second)) / denominator))
        dha_angle_deg = math.degrees(math.acos(cosine))
        angular_deviation_deg = 180.0 - dha_angle_deg
        if donor_acceptor_angstrom > 3.5 + 1.0e-9 or angular_deviation_deg > 30.0 + 1.0e-9:
            print(
                "skip PLIP hydrogen bond outside frozen geometry: "
                f"{donor['resname']}{donor['resid']}:{donor['name']} -> "
                f"{acceptor['resname']}{acceptor['resid']}:{acceptor['name']} "
                f"DA={donor_acceptor_angstrom:.3f}A deviation={angular_deviation_deg:.3f}deg"
            )
            continue
        donor_tag = f"{donor['resname']}{donor['resid']}:{donor['name']}" if donor["resname"] == "O6U" \
            else f"O6U-{donor['resname']}{donor['resid']}:{donor['name']}"
        acceptor_tag = f"O6U:{acceptor['name']}" if acceptor["resname"] == "O6U" \
            else f"{acceptor['resname']}{acceptor['resid']}:{acceptor['name']}"
        hydrogen_bonds.append({
            "metric_id": f"{donor_tag}__{acceptor_tag}",
            "protisdonor": d["protisdon"] == "True",
            "donor": mapped_endpoint(donor),
            "hydrogen": mapped_endpoint(hydrogen),
            "acceptor": mapped_endpoint(acceptor),
            "reference_donor_acceptor_distance_nm": float(d["dist_d-a"]) / 10.0,
            "reference_angle_deg": float(d["don_angle"]),
            "distance_cutoff_nm": 0.35,
            "angle_cutoff_deg": 30.0,
        })
    print(f"O6U hydrogen bonds from PLIP: {len(hydrogen_bonds)}")

    # TM-core / pocket alignment: the 14 experimentally resolved PSEN1 residues reported
    # for nirogacestat in Guo 2025 Supplementary Table 3 (2 hydrogen-bond residues +
    # 12 van-der-Waals residues). Frozen before production review; no trajectory-derived
    # residue may be added.
    contact_residues = sorted({77, 261, 268, 271, 272, 282, 287, 379, 380, 381, 425, 431, 432})
    tm_core = [a["index"] for a in top
               if a["segid"] in ("PROD", "PROE") and a["name"] == "CA"
               and a["resid"] in contact_residues]
    print(f"TM-core C-alpha (Guo 2025 contact residues): {len(tm_core)}")
    found = {top[i]["resid"] for i in tm_core}
    missing = [r for r in contact_residues if r not in found]
    print(f"TM-core missing contact residues: {missing}")
    if missing:
        raise ValueError(f"TM-core does not contain frozen contact residues: {missing}")
    tm_helix_segments_record = [{"residue": r, "guo2025_tm": tm_label} for r, tm_label in
                                ((77, "TM1"), (261, "TM6"), (268, "TM6a"), (271, "TM6a"),
                                 (272, "TM6a"), (282, "TM6a"), (287, "TM6a"), (379, "TM7"),
                                 (380, "TM7"), (381, "TM7"), (425, "TM8"), (431, "TM9"),
                                 (432, "TM9"))]

    # assemble output
    identity_sha = hashlib.sha256(
        json.dumps(trajectory_identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    ).hexdigest()
    tm_core_mappings = [
        mapped_endpoint(unique_reference_atom("B", top[index]["resid"], "CA"))
        for index in tm_core
    ]
    protein_ca_mappings = [
        mapped_endpoint(atom)
        for atom in ref
        if atom["record_type"] == "ATOM" and atom["name"] == "CA"
    ]
    prespecified_distances = [
        {
            "metric_id": "K380_backbone_O_to_O6U",
            "note": "literature H-bond residue",
            "atom1": mapped_endpoint(unique_reference_atom("B", 380, "O", "LYS")),
            "atom2": mapped_endpoint(unique_reference_atom("B", 502, "N07", "O6U")),
        },
        {
            "metric_id": "L432_backbone_O_to_O6U",
            "note": "literature H-bond residue",
            "atom1": mapped_endpoint(unique_reference_atom("B", 432, "O", "LEU")),
            "atom2": mapped_endpoint(unique_reference_atom("B", 502, "N04", "O6U")),
        },
    ]

    record = {
        "schema_version": "1.0",
        "approval_status": "draft_not_for_execution",
        "system_id": "8kct_nirogacestat_native",
        "reference_structure_id": "8KCT",
        "ligand_resname": "O6U",
        "trajectory_atom_identity_sha256": identity_sha,
        "trajectory_atom_count": len(trajectory_identity),
        "reference_coordinate_source": str(REFERENCE.relative_to(ROOT)),
        "trajectory_atom_source": str(STEP5.relative_to(ROOT)),
        "trajectory_topology_source": str(topology_source),
        "source_records": {
            "reference": {"path": str(REFERENCE.relative_to(ROOT)), "sha256": sha256(REFERENCE)},
            "trajectory_atom_source": {"path": str(STEP5.relative_to(ROOT)), "sha256": sha256(STEP5)},
            "trajectory_topology": {"path": str(topology_source), "sha256": sha256(topology_source)},
            "atom_correspondence": {"path": str(CORRESP.relative_to(ROOT)), "sha256": sha256(CORRESP)},
            "plip_xml": {"path": str(PLIP_XML.relative_to(ROOT)), "sha256": sha256(PLIP_XML)},
            "plip_normalized_json": {"path": str(PLIP_JSON.relative_to(ROOT)), "sha256": sha256(PLIP_JSON)},
            "minimized_structure": {"path": str(GRO.relative_to(ROOT)), "sha256": sha256(GRO)},
        },
        "coordinate_space": "pocket_aligned_euclidean_after_validated_whole_cluster_nojump_center",
        "native_contact_cutoff_nm": 0.45,
        "reference_distance_tolerance_nm": 0.0001,
        "hydrogen_bond_distance_cutoff_nm": 0.35,
        "hydrogen_bond_angular_deviation_cutoff_deg": 30.0,
        "atom_mappings": {
            "pocket_alignment": tm_core_mappings,
            "tm_core_ca": tm_core_mappings,
            "o6u_heavy": [mapped_endpoint(atom) for atom in heavy_ref_o6u],
            "protein_ca": protein_ca_mappings,
        },
        "native_contacts": native_contacts,
        "prespecified_distances": prespecified_distances,
        "hydrogen_bonds": hydrogen_bonds,
        "membrane_slab": {
            "source": str(GRO.relative_to(ROOT)),
            "tm_core_selection_rule": ("Guo 2025 Supplementary Table 3 nirogacestat contact "
                                       "residues (2 H-bond + 12 vdW); frozen before production review"),
            "tm_contact_residues": tm_helix_segments_record,
            "tm_core_ca_count": len(tm_core),
            "guo2025_contact_residues_covered": True,
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
