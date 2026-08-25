#!/usr/bin/env python3
"""Prepare immutable, membrane-centered 300-frame endpoint-energy inputs."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


REPLICAS = ("rep01", "rep02", "rep03")
TARGETS_NS = tuple(200.5 + index for index in range(300))
BLOCKS_NS = ((200.0, 260.0), (260.0, 320.0), (320.0, 380.0), (380.0, 440.0), (440.0, 500.0))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def select_midpoint_frames(times_ns: Iterable[float]) -> list[dict[str, Any]]:
    times = [float(value) for value in times_ns]
    if not times or any(not math.isfinite(value) for value in times):
        raise ValueError("trajectory times must be finite and nonempty")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("trajectory times must be strictly increasing")
    selected: list[dict[str, Any]] = []
    for stratum, target in enumerate(TARGETS_NS):
        insertion = bisect.bisect_left(times, target)
        candidates = []
        if insertion < len(times):
            candidates.append(insertion)
        if insertion > 0:
            candidates.append(insertion - 1)
        if not candidates:
            raise ValueError(f"no frame is available for target {target}")
        source_index = min(candidates, key=lambda index: (abs(times[index] - target), times[index]))
        block_index = min(int((target - 200.0) // 60.0), 4)
        selected.append({
            "stratum_index_zero_based": stratum,
            "target_time_ns": target,
            "source_index_zero_based": source_index,
            "source_time_ns": times[source_index],
            "absolute_offset_ns": abs(times[source_index] - target),
            "block_index": block_index,
            "block_start_ns": BLOCKS_NS[block_index][0],
            "block_end_ns": BLOCKS_NS[block_index][1],
        })
    indices = [row["source_index_zero_based"] for row in selected]
    if len(set(indices)) != 300:
        raise ValueError("midpoint selection does not yield 300 unique source frames")
    if [sum(row["block_index"] == index for row in selected) for index in range(5)] != [60] * 5:
        raise ValueError("fixed block membership is not 60 frames per block")
    return selected


def render_frame_index(selected: list[dict[str, Any]]) -> str:
    if len(selected) != 300:
        raise ValueError("frame index requires exactly 300 selected frames")
    values = [str(int(row["source_index_zero_based"])) for row in selected]
    return "[ endpoint_midpoint_frames_zero_based ]\n" + "\n".join(values) + "\n"


def preparation_output_names(replica: str) -> dict[str, str]:
    if replica not in REPLICAS:
        raise ValueError(f"unknown replica: {replica}")
    return {
        "structure": f"{replica}_endpoint_structure.gro",
        "reference": f"{replica}_endpoint_complex_reference.pdb",
        "trajectory": f"{replica}_endpoint_300frames_midplane0.xtc",
        "canary_trajectory": f"{replica}_endpoint_canary_3frames_midplane0.xtc",
    }


def parse_ndx(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    current: str | None = None
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith((";", "#")):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            if not current or current in groups:
                raise ValueError(f"{path}:{number}: empty or duplicate group")
            groups[current] = []
            continue
        if current is None:
            raise ValueError(f"{path}:{number}: atom indices precede a group")
        groups[current].extend(int(value) for value in line.split())
    return groups


def validate_endpoint_groups(
    groups: dict[str, list[int]], residue_counts: dict[str, int]
) -> dict[str, Any]:
    required = {"Receptor", "Ligand_O6U", "Complex"}
    if set(groups) != required:
        raise ValueError(f"endpoint index groups must be exactly {sorted(required)}")
    normalized: dict[str, list[int]] = {}
    for name, atoms in groups.items():
        if not atoms or any(not isinstance(value, int) or value <= 0 for value in atoms):
            raise ValueError(f"{name}: atom indices must be positive integers")
        if len(atoms) != len(set(atoms)):
            raise ValueError(f"{name}: duplicate atom indices")
        normalized[name] = sorted(atoms)
    receptor = set(normalized["Receptor"])
    ligand = set(normalized["Ligand_O6U"])
    complex_atoms = set(normalized["Complex"])
    if len(ligand) != 76:
        raise ValueError("Ligand_O6U must contain exactly 76 atoms")
    if receptor & ligand:
        raise ValueError("Receptor and Ligand_O6U groups overlap")
    if receptor | ligand != complex_atoms:
        raise ValueError("Complex must be the exact receptor/ligand union")
    expected = {"O6U": 1, "NAG": 18, "BMA": 3, "CLR": 3}
    for name, count in expected.items():
        if int(residue_counts.get(name, 0)) != count:
            raise ValueError(f"{name}: expected {count} retained residues")
    structural_lipids = int(residue_counts.get("PC1", 0)) + int(residue_counts.get("DSPC", 0))
    if structural_lipids != 2:
        raise ValueError("PC1/DSPC: expected exactly two retained structural lipids")
    return {
        "receptor_atom_count": len(receptor),
        "ligand_atom_count": len(ligand),
        "complex_atom_count": len(complex_atoms),
        "structural_lipid_residue_count": structural_lipids,
        "residue_counts": {key: int(value) for key, value in residue_counts.items()},
    }


def render_atom_index(groups: dict[str, list[int]]) -> str:
    lines: list[str] = []
    for name in ("Receptor", "Ligand_O6U", "Complex"):
        lines.append(f"[ {name} ]")
        atoms = groups[name]
        for start in range(0, len(atoms), 15):
            lines.append(" ".join(str(value) for value in atoms[start:start + 15]))
    return "\n".join(lines) + "\n"


def compare_selected_distances(
    source: list[tuple[float, float]],
    derived: list[tuple[float, float]],
    tolerance_nm: float,
) -> dict[str, Any]:
    if tolerance_nm <= 0:
        raise ValueError("tolerance must be positive")
    source_times = [row[0] for row in source]
    derived_times = [row[0] for row in derived]
    if source_times != derived_times:
        raise ValueError("source and derived distance times differ")
    if not source:
        raise ValueError("distance series are empty")
    differences = []
    for (_, left), (_, right) in zip(source, derived):
        if not (math.isfinite(left) and math.isfinite(right) and left >= 0 and right >= 0):
            raise ValueError("distance series contain non-finite or negative values")
        differences.append(abs(left - right))
    maximum = max(differences)
    return {
        "metric": "selected_frame_minimum_image_protein_O6U_heavy_atom_distance",
        "frame_count": len(source),
        "tolerance_nm": tolerance_nm,
        "maximum_absolute_difference_nm": maximum,
        "status": "pass" if maximum <= tolerance_nm else "fail",
    }


def write_new_json(path: Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _mapping_indices(entries: list[dict[str, Any]], atom_count: int, label: str) -> list[int]:
    values = [int(entry.get("index", -1)) for entry in entries]
    if len(values) < 2 or len(values) != len(set(values)) or any(value < 0 or value >= atom_count for value in values):
        raise ValueError(f"{label}: frozen atom mapping is invalid")
    return values


def _minimum_distance_nm(first: Any, second: Any, box: Any) -> float:
    import numpy as np
    from MDAnalysis.lib.distances import distance_array

    matrix = distance_array(first.positions, second.positions, box=box)
    value = float(np.min(matrix)) / 10.0
    if not math.isfinite(value) or value < 0:
        raise ValueError("minimum distance is invalid")
    return value


def execute_preparation(args: argparse.Namespace) -> dict[str, Any]:
    import MDAnalysis as mda
    import numpy as np

    from analyze_membrane_qc_mdanalysis import leaflet_relative_z_nm

    if args.replica not in REPLICAS:
        raise ValueError(f"replica must be one of {REPLICAS}")
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {args.output_dir}")
    gate = _load_json(args.all_three_gate)
    if gate.get("status") != "pass" or gate.get("eligible_replicas") != list(REPLICAS):
        raise ValueError("all-three eligibility gate is not passing")
    if gate.get("production_runners_active") not in (False, 0, []):
        raise ValueError("endpoint preprocessing cannot overlap production")
    pbc = _load_json(args.pbc_report)
    if pbc.get("status") != "pass" or float(pbc.get("maximum_absolute_difference_nm", math.inf)) > 0.01:
        raise ValueError("PBC invariance report is not passing at 0.01 nm")
    provenance = _load_json(args.trajectory_provenance)
    source_hash = file_sha256(args.source_trajectory)
    expected_source_hash = provenance.get("outputs", {}).get("center_and_rebox", {}).get("sha256")
    if expected_source_hash != source_hash:
        raise ValueError("centered/reboxed trajectory hash differs from provenance")
    if file_sha256(args.tpr) != args.expected_tpr_sha256:
        raise ValueError("production TPR hash mismatch")

    mapping = _load_json(args.membrane_mapping)
    frozen = mapping.get("frozen_atom_groups", {})
    universe = mda.Universe(str(args.tpr), str(args.source_trajectory))
    upper_indices = _mapping_indices(frozen.get("upper_leaflet_phosphate_atoms", []), len(universe.atoms), "upper leaflet")
    lower_indices = _mapping_indices(frozen.get("lower_leaflet_phosphate_atoms", []), len(universe.atoms), "lower leaflet")

    receptor = universe.select_atoms("protein or resname NAG BMA CLR PC1 DSPC")
    ligand = universe.select_atoms("resname O6U")
    complex_indices = np.unique(np.concatenate([receptor.indices, ligand.indices]))
    groups = {
        "Receptor": [int(value) + 1 for value in receptor.indices],
        "Ligand_O6U": [int(value) + 1 for value in ligand.indices],
        "Complex": [int(value) + 1 for value in complex_indices],
    }
    residue_counts = Counter(str(name) for name in universe.residues.resnames)
    group_report = validate_endpoint_groups(groups, dict(residue_counts))
    times_ns = [float(ts.time) / 1000.0 for ts in universe.trajectory]
    selected = select_midpoint_frames(times_ns)

    args.output_dir.mkdir(parents=True)
    frame_index_path = args.output_dir / "endpoint_frames.ndx"
    frame_index_path.write_text(render_frame_index(selected), encoding="utf-8")
    endpoint_index_path = args.output_dir / "endpoint_groups.ndx"
    endpoint_index_path.write_text(render_atom_index(groups), encoding="utf-8")
    output_names = preparation_output_names(args.replica)
    structure_path = args.output_dir / output_names["structure"]
    reference_path = args.output_dir / output_names["reference"]
    trajectory_path = args.output_dir / output_names["trajectory"]
    canary_trajectory_path = args.output_dir / output_names["canary_trajectory"]
    frame_map_path = args.output_dir / "frame_map.csv"
    source_distances: list[tuple[float, float]] = []
    shifts: list[float] = []
    protein_heavy = universe.select_atoms("protein and not name H*")
    ligand_heavy = universe.select_atoms("resname O6U and not name H*")
    canary_output_indices = {0, 150, 299}
    with (
        mda.Writer(str(trajectory_path), n_atoms=len(universe.atoms)) as writer,
        mda.Writer(str(canary_trajectory_path), n_atoms=len(universe.atoms)) as canary_writer,
    ):
        for output_index, row in enumerate(selected):
            universe.trajectory[int(row["source_index_zero_based"])]
            box_z = float(universe.dimensions[2])
            if not math.isfinite(box_z) or box_z <= 0:
                raise ValueError("trajectory frame has invalid box z length")
            source_distances.append((float(row["source_time_ns"]), _minimum_distance_nm(protein_heavy, ligand_heavy, universe.dimensions)))
            _, _, midplane = leaflet_relative_z_nm(
                universe.atoms[upper_indices].positions[:, 2],
                universe.atoms[lower_indices].positions[:, 2],
                box_z,
            )
            universe.atoms.positions[:, 2] -= midplane
            shifts.append(float(midplane))
            writer.write(universe.atoms)
            if output_index in canary_output_indices:
                canary_writer.write(universe.atoms)
            if output_index == 0:
                universe.atoms.write(str(structure_path))
                universe.atoms[complex_indices].write(str(reference_path))

    derived = mda.Universe(str(args.tpr), str(trajectory_path))
    derived_protein = derived.select_atoms("protein and not name H*")
    derived_ligand = derived.select_atoms("resname O6U and not name H*")
    derived_distances = []
    for row, ts in zip(selected, derived.trajectory):
        derived_distances.append((float(row["source_time_ns"]), _minimum_distance_nm(derived_protein, derived_ligand, ts.dimensions)))
    invariance = compare_selected_distances(source_distances, derived_distances, 0.01)
    if invariance["status"] != "pass":
        raise ValueError("derived midpoint trajectory failed selected-frame distance invariance")

    with frame_map_path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(selected[0]) + ["output_index_zero_based", "midplane_shift_angstrom"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for output_index, (row, shift) in enumerate(zip(selected, shifts)):
            writer.writerow({**row, "output_index_zero_based": output_index, "midplane_shift_angstrom": shift})

    outputs = {}
    for path in (
        frame_index_path,
        endpoint_index_path,
        structure_path,
        reference_path,
        trajectory_path,
        canary_trajectory_path,
        frame_map_path,
    ):
        outputs[path.name] = {"bytes": path.stat().st_size, "sha256": file_sha256(path)}
    manifest = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_preparation",
        "status": "pass",
        "replica": args.replica,
        "frame_count": 300,
        "window_ns": [200.0, 500.0],
        "target_times_ns": [200.5, 499.5],
        "fixed_block_counts": [60, 60, 60, 60, 60],
        "canary_output_indices_zero_based": [0, 150, 299],
        "canary_target_times_ns": [200.5, 350.5, 499.5],
        "membrane_normal_axis": "z",
        "membrane_midplane_z_angstrom": 0.0,
        "coordinate_operation": "per_frame_global_z_translation_only_after_validated_PBC_chain",
        "groups": group_report,
        "pbc_invariance": invariance,
        "source": {
            "tpr": {"path": str(args.tpr), "sha256": file_sha256(args.tpr)},
            "trajectory": {"path": str(args.source_trajectory), "sha256": source_hash},
            "trajectory_provenance": {"path": str(args.trajectory_provenance), "sha256": file_sha256(args.trajectory_provenance)},
            "pbc_report": {"path": str(args.pbc_report), "sha256": file_sha256(args.pbc_report)},
            "membrane_mapping": {"path": str(args.membrane_mapping), "sha256": file_sha256(args.membrane_mapping)},
            "all_three_gate": {"path": str(args.all_three_gate), "sha256": file_sha256(args.all_three_gate)},
        },
        "outputs": outputs,
        "raw_inputs_immutable": True,
    }
    manifest_path = args.output_dir / "PREPARATION_MANIFEST.json"
    write_new_json(manifest_path, manifest)
    digest = file_sha256(manifest_path)
    manifest_path.with_suffix(".json.sha256").write_text(f"{digest}  {manifest_path.name}\n", encoding="ascii")
    return {"status": "pass", "replica": args.replica, "manifest": str(manifest_path), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica", required=True, choices=REPLICAS)
    parser.add_argument("--tpr", type=Path, required=True)
    parser.add_argument("--expected-tpr-sha256", required=True)
    parser.add_argument("--source-trajectory", type=Path, required=True)
    parser.add_argument("--trajectory-provenance", type=Path, required=True)
    parser.add_argument("--pbc-report", type=Path, required=True)
    parser.add_argument("--membrane-mapping", type=Path, required=True)
    parser.add_argument("--all-three-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = execute_preparation(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
