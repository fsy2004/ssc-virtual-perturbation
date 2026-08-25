#!/usr/bin/env python3
"""Fail-closed identity and stereochemistry audit for an O6U CREST ensemble."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem

from prepare_o6u_crest_input import (
    EXPECTED_ATOMS,
    EXPECTED_CHIRAL_CENTERS,
    load_single_sdf,
    validate_identity,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_xyz_ensemble(path: Path) -> list[dict[str, object]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    frames = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ atom count at line {cursor + 1}") from exc
        if atom_count != EXPECTED_ATOMS:
            raise ValueError(f"Frame {len(frames) + 1} has {atom_count} atoms, expected {EXPECTED_ATOMS}")
        if cursor + atom_count + 1 >= len(lines):
            raise ValueError(f"Truncated XYZ frame {len(frames) + 1}")
        comment = lines[cursor + 1]
        elements = []
        coordinates = []
        for offset in range(atom_count):
            fields = lines[cursor + 2 + offset].split()
            if len(fields) < 4:
                raise ValueError(f"Malformed XYZ coordinate at frame {len(frames) + 1}, atom {offset + 1}")
            try:
                xyz = [float(fields[1]), float(fields[2]), float(fields[3])]
            except ValueError as exc:
                raise ValueError(f"Non-numeric XYZ coordinate at frame {len(frames) + 1}, atom {offset + 1}") from exc
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError(f"Non-finite coordinate at frame {len(frames) + 1}, atom {offset + 1}")
            elements.append(fields[0])
            coordinates.append(xyz)
        frames.append({"comment": comment, "elements": elements, "coordinates": np.asarray(coordinates)})
        cursor += atom_count + 2
    if not frames:
        raise ValueError("CREST ensemble contains no frames")
    return frames


def signed_chiral_volume(coordinates: np.ndarray, center: int, neighbors: list[int]) -> float:
    if len(neighbors) != 4:
        raise ValueError(f"Chiral center {center} does not have four explicit neighbours")
    points = coordinates[neighbors]
    return float(np.linalg.det(np.stack((points[0] - points[3], points[1] - points[3], points[2] - points[3]))))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--ensemble", required=True, type=Path)
    parser.add_argument("--crest-log", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    source = args.source_sdf.resolve()
    ensemble = args.ensemble.resolve()
    report_path = args.report.resolve()
    if not source.is_file() or not ensemble.is_file():
        raise SystemExit("Source SDF and CREST ensemble must both exist")
    if report_path.exists():
        raise SystemExit("Refusing to overwrite an existing report")

    mol = load_single_sdf(source)
    identity = validate_identity(mol)
    expected_elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    reference = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    bonds = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds())
    reference_lengths = {(i, j): float(np.linalg.norm(reference[i] - reference[j])) for i, j in bonds}
    chiral_neighbors = {
        center: sorted(neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(center).GetNeighbors())
        for center in EXPECTED_CHIRAL_CENTERS
    }
    reference_signs = {}
    for center, neighbors in chiral_neighbors.items():
        volume = signed_chiral_volume(reference, center, neighbors)
        if abs(volume) <= 1.0e-4:
            raise SystemExit(f"Reference chiral center {center} is geometrically degenerate")
        reference_signs[center] = math.copysign(1.0, volume)

    frames = read_xyz_ensemble(ensemble)
    frame_summaries = []
    for frame_index, frame in enumerate(frames, start=1):
        if frame["elements"] != expected_elements:
            raise SystemExit(f"Frame {frame_index} element sequence differs from the immutable SDF")
        coordinates = frame["coordinates"]
        bond_ratios = []
        for i, j in bonds:
            observed = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            ratio = observed / reference_lengths[(i, j)]
            bond_ratios.append(ratio)
            if not 0.70 <= ratio <= 1.35:
                raise SystemExit(
                    f"Frame {frame_index} bond {i}-{j} ratio {ratio:.4f} falls outside the topology guard"
                )
        minimum_abs_chiral_volume = math.inf
        for center, neighbors in chiral_neighbors.items():
            volume = signed_chiral_volume(coordinates, center, neighbors)
            minimum_abs_chiral_volume = min(minimum_abs_chiral_volume, abs(volume))
            if abs(volume) <= 1.0e-4 or math.copysign(1.0, volume) != reference_signs[center]:
                raise SystemExit(f"Frame {frame_index} inverts or degenerates O6U chiral center {center}")
        frame_summaries.append(
            {
                "frame": frame_index,
                "comment": frame["comment"],
                "minimum_bond_length_ratio_to_source": min(bond_ratios),
                "maximum_bond_length_ratio_to_source": max(bond_ratios),
                "minimum_absolute_signed_chiral_volume_angstrom3": minimum_abs_chiral_volume,
            }
        )

    crest_log = args.crest_log.resolve() if args.crest_log else None
    if crest_log is not None and (not crest_log.is_file() or crest_log.stat().st_size == 0):
        raise SystemExit("Declared CREST log is missing or empty")
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_crest_ensemble_validation",
        "status": "pass",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sdf": str(source),
        "source_sdf_sha256": sha256(source),
        "ensemble": str(ensemble),
        "ensemble_sha256": sha256(ensemble),
        "crest_log": str(crest_log) if crest_log else None,
        "crest_log_sha256": sha256(crest_log) if crest_log else None,
        "identity": identity,
        "frame_count": len(frames),
        "frame_summaries": frame_summaries,
        "release_boundary": "Conformer-start identity only; no force-field parameter is approved by this report.",
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "frame_count": len(frames)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
