#!/usr/bin/env python3
"""Select a deterministic O6U QM conformer set from the validated CREST ensemble.

The selector uses RDKit torsion-fingerprint-distance (TFD) clustering with one
prespecified cutoff.  Every Butina cluster contributes its centroid, and the
CREST global minimum is retained even when it is not a cluster centroid.  The
deposited 8KCT heavy-atom geometry is emitted separately as an experimental
mapping/reference structure; it is not silently promoted to a hydrogen-complete
QM starting geometry.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints
from rdkit.ML.Cluster import Butina

from prepare_o6u_crest_input import EXPECTED_ATOMS, load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble


FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def comment_energy(comment: str) -> float:
    values = FLOAT_RE.findall(comment)
    if not values:
        raise ValueError(f"CREST frame comment contains no numeric energy: {comment!r}")
    value = float(values[0])
    if not math.isfinite(value):
        raise ValueError(f"Non-finite CREST frame energy: {comment!r}")
    return value


def xyz_text(elements: list[str], coordinates: np.ndarray, comment: str) -> str:
    rows = [str(len(elements)), comment]
    rows.extend(
        f"{element:<2s} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for element, (x, y, z) in zip(elements, coordinates, strict=True)
    )
    return "\n".join(rows) + "\n"


def native_heavy_reference(path: Path, source: Chem.Mol) -> tuple[list[str], np.ndarray, list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    heavy = [row for row in rows if row["hydrogen_or_heavy"] == "heavy"]
    heavy.sort(key=lambda row: int(row["sdf_zero_based_index"]))
    expected = [atom for atom in source.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy) != len(expected):
        raise ValueError(f"Native correspondence has {len(heavy)} heavy atoms, expected {len(expected)}")
    elements: list[str] = []
    names: list[str] = []
    coordinates: list[list[float]] = []
    for row, atom in zip(heavy, expected, strict=True):
        index = int(row["sdf_zero_based_index"])
        if index != atom.GetIdx() or row["element"] != atom.GetSymbol():
            raise ValueError("Native 8KCT correspondence is not in immutable SDF heavy-atom order")
        if row["native_8kct_present"].lower() != "yes":
            raise ValueError(f"Native coordinate missing for {row['ccd_atom_id']}")
        xyz = [float(row[key]) for key in ("native_x_angstrom", "native_y_angstrom", "native_z_angstrom")]
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Non-finite native coordinate for {row['ccd_atom_id']}")
        elements.append(row["element"])
        names.append(row["ccd_atom_id"])
        coordinates.append(xyz)
    return elements, np.asarray(coordinates, dtype=float), names


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--ensemble", required=True, type=Path)
    parser.add_argument("--native-correspondence", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tfd-cutoff", type=float, default=0.20)
    args = parser.parse_args()

    source_path = args.source_sdf.resolve()
    ensemble_path = args.ensemble.resolve()
    correspondence_path = args.native_correspondence.resolve()
    output_dir = args.output_dir.resolve()
    for path in (source_path, ensemble_path, correspondence_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    if not 0.0 < args.tfd_cutoff < 1.0:
        raise SystemExit("TFD cutoff must be strictly between zero and one")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {output_dir}")

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    frames = read_xyz_ensemble(ensemble_path)
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if len(frames) < 2:
        raise SystemExit("At least two validated CREST frames are required for clustering")

    ensemble_mol = Chem.Mol(source)
    ensemble_mol.RemoveAllConformers()
    energies: list[float] = []
    for frame_number, frame in enumerate(frames, start=1):
        if frame["elements"] != expected_elements:
            raise SystemExit(f"Frame {frame_number} element order differs from immutable source SDF")
        conformer = Chem.Conformer(EXPECTED_ATOMS)
        conformer.Set3D(True)
        for atom_index, (x, y, z) in enumerate(frame["coordinates"]):
            conformer.SetAtomPosition(atom_index, (float(x), float(y), float(z)))
        ensemble_mol.AddConformer(conformer, assignId=True)
        energies.append(comment_energy(str(frame["comment"])))

    distances = TorsionFingerprints.GetTFDMatrix(ensemble_mol)
    clusters = Butina.ClusterData(
        distances,
        len(frames),
        args.tfd_cutoff,
        isDistData=True,
        reordering=True,
    )
    if not clusters or sorted(index for cluster in clusters for index in cluster) != list(range(len(frames))):
        raise SystemExit("TFD clustering did not partition every CREST conformer exactly once")

    global_minimum = min(range(len(frames)), key=lambda index: (energies[index], index))
    roles: dict[int, list[str]] = {cluster[0]: [f"tfd_cluster_{number:02d}_centroid"] for number, cluster in enumerate(clusters, 1)}
    roles.setdefault(global_minimum, []).append("crest_global_minimum")
    selected_indices = sorted(roles, key=lambda index: (energies[index], index))

    output_dir.mkdir(parents=True, exist_ok=True)
    selected_dir = output_dir / "selected_xyz"
    selected_dir.mkdir()
    selected_manifest = []
    combined_blocks = []
    for selection_number, index in enumerate(selected_indices, start=1):
        frame = frames[index]
        comment = (
            f"O6U CREST frame={index + 1}; energy_raw={energies[index]:.12f}; "
            f"roles={','.join(roles[index])}; source_ensemble_sha256={sha256(ensemble_path)}"
        )
        text = xyz_text(expected_elements, frame["coordinates"], comment)
        output = selected_dir / f"o6u_qm_rep_{selection_number:02d}_frame_{index + 1:04d}.xyz"
        output.write_text(text, encoding="utf-8", newline="\n")
        combined_blocks.append(text)
        selected_manifest.append(
            {
                "selection_number": selection_number,
                "crest_frame_1based": index + 1,
                "energy_raw": energies[index],
                "roles": roles[index],
                "output_xyz": str(output),
                "output_xyz_sha256": sha256(output),
            }
        )

    combined = output_dir / "O6U_QM_REPRESENTATIVES.xyz"
    combined.write_text("".join(combined_blocks), encoding="utf-8", newline="\n")

    native_elements, native_coordinates, native_names = native_heavy_reference(correspondence_path, source)
    native = output_dir / "O6U_8KCT_NATIVE_HEAVY_REFERENCE.xyz"
    native.write_text(
        xyz_text(
            native_elements,
            native_coordinates,
            "8KCT deposited O6U heavy atoms only; mapping/reference, not a QM start; atom_names=" + ",".join(native_names),
        ),
        encoding="utf-8",
        newline="\n",
    )

    cluster_membership = []
    for cluster_number, cluster in enumerate(clusters, start=1):
        for within_cluster_order, index in enumerate(cluster, start=1):
            cluster_membership.append(
                {
                    "crest_frame_1based": index + 1,
                    "energy_raw": energies[index],
                    "tfd_cluster": cluster_number,
                    "within_cluster_order": within_cluster_order,
                    "is_centroid": within_cluster_order == 1,
                    "is_global_minimum": index == global_minimum,
                    "selected": index in roles,
                }
            )
    membership_path = output_dir / "O6U_TFD_CLUSTER_MEMBERSHIP.tsv"
    with membership_path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(cluster_membership[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(cluster_membership, key=lambda row: int(row["crest_frame_1based"])))

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_qm_representative_selection",
        "status": "pass",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "source_sdf": str(source_path),
        "source_sdf_sha256": sha256(source_path),
        "crest_ensemble": str(ensemble_path),
        "crest_ensemble_sha256": sha256(ensemble_path),
        "native_correspondence": str(correspondence_path),
        "native_correspondence_sha256": sha256(correspondence_path),
        "method": {
            "distance": "RDKit torsion fingerprint deviation (TFD)",
            "clustering": "Butina on condensed TFD distance matrix with reordering=True",
            "tfd_cutoff": args.tfd_cutoff,
            "selection_rule": "every cluster centroid plus CREST global minimum; duplicate roles retained on one frame",
            "post_qm_duplicate_rule": "collapse only if heavy-atom RMSD < 0.1 A AND absolute energy difference < 0.1 kcal/mol",
        },
        "ensemble_frame_count": len(frames),
        "cluster_count": len(clusters),
        "cluster_sizes": [len(cluster) for cluster in clusters],
        "cluster_centroids_1based": [cluster[0] + 1 for cluster in clusters],
        "crest_global_minimum_frame_1based": global_minimum + 1,
        "selected_count": len(selected_manifest),
        "selected": selected_manifest,
        "combined_xyz": str(combined),
        "combined_xyz_sha256": sha256(combined),
        "cluster_membership_tsv": str(membership_path),
        "cluster_membership_tsv_sha256": sha256(membership_path),
        "native_heavy_reference_xyz": str(native),
        "native_heavy_reference_xyz_sha256": sha256(native),
        "release_boundary": "Representative starts only; no QM target or force-field parameter is approved by this report.",
    }
    report_path = output_dir / "O6U_QM_REPRESENTATIVE_SELECTION.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "clusters": len(clusters), "selected": len(selected_manifest)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
