#!/usr/bin/env python3
"""Unsmoothed membrane/cell QC for the full centered 0-500 ns trajectories.

This script deliberately does not report Lx*Ly/N as area per lipid.  The box
lateral area is cell QC only.  Protein-aware APL and POPC order parameters are
external pre-production NO-GO requirements until separately validated atom
mappings/tools exist; their absence can never yield a complete primary bundle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance

from primary_postprocessing_common import (
    PRIMARY_WINDOW_NS,
    REALIZATION_IDS,
    MEMBRANE_STATIONARITY_METRICS,
    ContractError,
    atomic_write_csv,
    atomic_write_json,
    block_diagnostics,
    check_mdanalysis_version,
    has_placeholder,
    load_json,
    primary_window_mask,
    require,
    resolve_record,
    robust_first_difference,
    sha256_file,
    stationarity_diagnostics,
    validate_primary_manifest,
    validate_time_axis,
    verify_atom,
)


RAW_FIELDS = [
    "system_id",
    "realization_id",
    "frame_index_zero_based",
    "time_ns",
    "in_primary_window_200_500_ns",
    "phosphate_peak_thickness_nm",
    "upper_phosphate_peak_z_relative_nm",
    "lower_phosphate_peak_z_relative_nm",
    "upper_leaflet_mismatch_count",
    "lower_leaflet_mismatch_count",
    "new_leaflet_flip_events",
    "cumulative_leaflet_flip_events",
    "cell_lateral_area_nm2_not_apl",
    "box_z_vector_length_nm",
    "cell_volume_nm3",
    "protein_tilt_deg",
    "outside_protein_hydrophobic_core_water_count",
    "outside_protein_hydrophobic_core_largest_water_cluster",
]


def topology_identity_sha256(universe: mda.Universe) -> str:
    records = []
    for atom in universe.atoms:
        record = {"index": int(atom.index), "name": str(atom.name), "resname": str(atom.resname), "resid": int(atom.resid)}
        for field in ("segid", "chainID"):
            try:
                value = getattr(atom, field)
            except Exception:
                continue
            if value not in (None, ""):
                record[field] = str(value)
        records.append(record)
    return hashlib.sha256(json.dumps(records, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")).hexdigest()


def atom_is_hydrogen(atom: Any) -> bool:
    """Identify hydrogens by chemical identity, independent of HMR mass."""
    try:
        element = str(atom.element).strip().upper()
    except Exception:
        element = ""
    normalized_name = str(atom.name).strip().upper().lstrip("0123456789")
    return element == "H" or normalized_name.startswith("H")


def frozen_group(universe: mda.Universe, entries: Sequence[Mapping[str, Any]], label: str, minimum: int = 1) -> np.ndarray:
    require(isinstance(entries, Sequence) and len(entries) >= minimum, f"{label} requires at least {minimum} explicit atoms")
    indices = []
    for number, expected in enumerate(entries):
        index = int(expected.get("index", -1))
        require(0 <= index < len(universe.atoms), f"{label}[{number}] index is out of range")
        verify_atom(universe.atoms[index], expected, f"{label}[{number}]")
        indices.append(index)
    require(len(indices) == len(set(indices)), f"{label} contains duplicate atom indices")
    return np.asarray(indices, dtype=np.int64)


def _circular_center_angstrom(z_angstrom: np.ndarray, box_z_angstrom: float, label: str) -> float:
    z = np.mod(np.asarray(z_angstrom, dtype=np.float64), box_z_angstrom)
    require(len(z) >= 2 and box_z_angstrom > 0.0, f"{label} circular center requires at least two atoms and a positive box length")
    angles = 2.0 * np.pi * z / box_z_angstrom
    vector = np.mean(np.exp(1j * angles))
    require(abs(vector) > 1e-8, f"{label} circular center is undefined")
    return float((np.angle(vector) % (2.0 * np.pi)) * box_z_angstrom / (2.0 * np.pi))


def _relative_to_midplane_nm(z_angstrom: np.ndarray, midplane_angstrom: float, box_z_angstrom: float) -> np.ndarray:
    z = np.mod(np.asarray(z_angstrom, dtype=np.float64), box_z_angstrom)
    center = midplane_angstrom % box_z_angstrom
    relative = (z - center + box_z_angstrom / 2.0) % box_z_angstrom - box_z_angstrom / 2.0
    return relative / 10.0


def leaflet_relative_z_nm(
    upper_z_angstrom: np.ndarray,
    lower_z_angstrom: np.ndarray,
    box_z_angstrom: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return leaflet positions relative to one PBC-safe phosphate midplane.

    Each frozen leaflet is circularly centered independently.  The midplane is
    then the midpoint on the minimum-image path from the frozen lower center to
    the frozen upper center.  This avoids the undefined combined circular mean
    of two nearly opposing leaflets and supplies the same physical midplane to
    both phosphate and water-defect calculations.
    """

    upper_center = _circular_center_angstrom(upper_z_angstrom, box_z_angstrom, "Upper phosphate leaflet")
    lower_center = _circular_center_angstrom(lower_z_angstrom, box_z_angstrom, "Lower phosphate leaflet")
    separation = (upper_center - lower_center + box_z_angstrom / 2.0) % box_z_angstrom - box_z_angstrom / 2.0
    ambiguity_tolerance_angstrom = max(1e-6, box_z_angstrom * 1e-8)
    require(
        ambiguity_tolerance_angstrom < separation < box_z_angstrom / 2.0 - ambiguity_tolerance_angstrom,
        "Frozen upper/lower phosphate centers are inverted, coincident, or half-box ambiguous",
    )
    midplane = (lower_center + separation / 2.0) % box_z_angstrom
    upper_relative = _relative_to_midplane_nm(upper_z_angstrom, midplane, box_z_angstrom)
    lower_relative = _relative_to_midplane_nm(lower_z_angstrom, midplane, box_z_angstrom)
    return upper_relative, lower_relative, float(midplane)


def gaussian_density(values_nm: np.ndarray, grid_nm: np.ndarray, bandwidth_nm: float) -> np.ndarray:
    values = np.asarray(values_nm, dtype=np.float64)
    require(len(values) >= 2 and bandwidth_nm > 0.0, "Density peak requires at least two atoms and positive bandwidth")
    scaled = (grid_nm[:, None] - values[None, :]) / bandwidth_nm
    density = np.mean(np.exp(-0.5 * scaled * scaled), axis=1) / (bandwidth_nm * math.sqrt(2.0 * math.pi))
    require(np.all(np.isfinite(density)), "Phosphate density contains NaN or infinity")
    return density


def _largest_periodic_cluster_3d(points_nm: np.ndarray, box_nm: np.ndarray, cutoff_nm: float) -> int:
    points = np.asarray(points_nm, dtype=np.float64)
    box = np.asarray(box_nm, dtype=np.float64)
    require(points.ndim == 2 and points.shape[1] == 3, "Water-cluster coordinates must be N x 3")
    require(box.shape == (3,) and np.all(np.isfinite(box)) and np.all(box > 0.0), "Water-cluster box is invalid")
    count = len(points)
    if count == 0:
        return 0
    parents = list(range(count))

    def find(value: int) -> int:
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value

    def union(first: int, second: int) -> None:
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    pairs = capped_distance(
        points,
        points,
        max_cutoff=float(cutoff_nm),
        min_cutoff=0.0,
        box=np.asarray([box[0], box[1], box[2], 90.0, 90.0, 90.0], dtype=np.float32),
        return_distances=False,
    )
    for first, second in np.asarray(pairs, dtype=np.int64):
        if first < second:
            union(int(first), int(second))
    sizes: dict[int, int] = {}
    for index in range(count):
        root = find(index)
        sizes[root] = sizes.get(root, 0) + 1
    return max(sizes.values())


def _outside_protein_water(
    water_positions_angstrom: np.ndarray,
    water_relative_z_nm: np.ndarray,
    protein_positions_angstrom: np.ndarray,
    protein_relative_z_nm: np.ndarray,
    lx_nm: float,
    ly_nm: float,
    lz_nm: float,
    core_half_thickness_nm: float,
    protein_exclusion_nm: float,
    cluster_cutoff_nm: float,
) -> tuple[int, int]:
    in_core = np.abs(water_relative_z_nm) <= core_half_thickness_nm
    water_positions = np.asarray(water_positions_angstrom, dtype=np.float64)
    protein_positions = np.asarray(protein_positions_angstrom, dtype=np.float64)
    protein_relative = np.asarray(protein_relative_z_nm, dtype=np.float64)
    require(len(water_positions) == len(water_relative_z_nm), "Water relative-z mapping differs")
    require(len(protein_positions) == len(protein_relative), "Protein relative-z mapping differs")
    candidates = np.column_stack((water_positions[in_core, :2] / 10.0, np.asarray(water_relative_z_nm, dtype=np.float64)[in_core]))
    if len(candidates) == 0:
        return 0, 0
    protein_xyz = np.column_stack((protein_positions[:, :2] / 10.0, protein_relative))
    close_pairs = capped_distance(
        candidates,
        protein_xyz,
        max_cutoff=float(protein_exclusion_nm),
        min_cutoff=0.0,
        box=np.asarray([lx_nm, ly_nm, lz_nm, 90.0, 90.0, 90.0], dtype=np.float32),
        return_distances=False,
    )
    excluded = np.zeros(len(candidates), dtype=bool)
    if len(close_pairs):
        excluded[np.asarray(close_pairs, dtype=np.int64)[:, 0]] = True
    outside_array = candidates[~excluded]
    if len(outside_array) == 0:
        return 0, 0
    return len(outside_array), _largest_periodic_cluster_3d(outside_array, np.asarray([lx_nm, ly_nm, lz_nm]), cluster_cutoff_nm)


def _validate_mapping(mapping: Mapping[str, Any], universe: mda.Universe) -> dict[str, Any]:
    require(mapping.get("schema_version") == "1.0", "Membrane mapping schema_version must be 1.0")
    require(mapping.get("approval_status") in {"approved", "synthetic_self_test"}, "Membrane mapping is not approved")
    require(not has_placeholder(mapping), "Membrane mapping contains TODO/REPLACE_ME placeholders")
    require(mapping.get("system_id") == "8kct_nirogacestat_native", "Membrane mapping system_id differs")
    require(mapping.get("box_geometry_required") == "orthorhombic", "Only the validated orthorhombic cell implementation is available")
    fingerprint = topology_identity_sha256(universe)
    require(mapping.get("trajectory_atom_identity_sha256") == fingerprint, "Trajectory atom-identity fingerprint differs from the frozen membrane mapping")
    groups = mapping.get("frozen_atom_groups", {})
    upper = frozen_group(universe, groups.get("upper_leaflet_phosphate_atoms", []), "upper_leaflet_phosphate_atoms", minimum=2)
    lower = frozen_group(universe, groups.get("lower_leaflet_phosphate_atoms", []), "lower_leaflet_phosphate_atoms", minimum=2)
    tilt_upper = frozen_group(universe, groups.get("protein_tilt_upper_anchor_atoms", []), "protein_tilt_upper_anchor_atoms", minimum=2)
    tilt_lower = frozen_group(universe, groups.get("protein_tilt_lower_anchor_atoms", []), "protein_tilt_lower_anchor_atoms", minimum=2)
    protein = frozen_group(universe, groups.get("protein_heavy_atoms", []), "protein_heavy_atoms", minimum=3)
    waters = frozen_group(universe, groups.get("water_oxygen_atoms", []), "water_oxygen_atoms", minimum=1)
    require(set(upper).isdisjoint(set(lower)), "Upper and lower phosphate groups overlap")
    require(set(protein).isdisjoint(set(waters)), "Protein and water groups overlap")
    names = np.asarray(universe.atoms.names)
    resnames = np.asarray(universe.atoms.resnames)
    phosphate_keys = {(str(universe.atoms[index].resname), str(universe.atoms[index].name)) for index in np.concatenate([upper, lower])}
    expected_phosphates = {int(index) for index in range(len(universe.atoms)) if (str(resnames[index]), str(names[index])) in phosphate_keys}
    require(set(int(index) for index in np.concatenate([upper, lower])) == expected_phosphates, "Frozen leaflet groups must cover every topology phosphate atom matching their frozen lipid resname/name identities")
    water_keys = {(str(universe.atoms[index].resname), str(universe.atoms[index].name)) for index in waters}
    expected_waters = {int(index) for index in range(len(universe.atoms)) if (str(resnames[index]), str(names[index])) in water_keys}
    require(set(int(index) for index in waters) == expected_waters, "Frozen water group must cover every topology water oxygen matching its frozen resname/name identities")
    protein_atoms = universe.select_atoms("protein")
    protein_masses = np.asarray(protein_atoms.masses, dtype=np.float64)
    require(np.all(np.isfinite(protein_masses)) and np.all(protein_masses > 0.0), "Protein atoms lack positive finite masses")
    expected_protein_heavy = {
        int(index) for index in protein_atoms.indices if not atom_is_hydrogen(universe.atoms[index])
    }
    require(set(int(index) for index in protein) == expected_protein_heavy, "Frozen protein_heavy_atoms must cover every topology protein heavy atom")
    availability = mapping.get("external_metrics", {})
    require(set(availability) == {"protein_aware_area_per_lipid", "popc_deuterium_order_parameters"}, "External membrane metric keys differ")
    for metric in ("protein_aware_area_per_lipid", "popc_deuterium_order_parameters"):
        record = availability.get(metric)
        require(isinstance(record, Mapping) and record.get("status") in {"not_available", "validated"}, f"{metric} status must be not_available or validated")
        if record.get("status") == "not_available":
            require(isinstance(record.get("reason"), str) and len(record["reason"]) >= 20, f"{metric} requires a concrete unavailability reason")
        else:
            require(record.get("frozen_tool_route_before_production") is True, f"{metric} tool route was not frozen before production")
            tool = record.get("tool", {})
            require(isinstance(tool, Mapping) and isinstance(tool.get("name"), str) and isinstance(tool.get("version"), str), f"{metric} tool identity is incomplete")
            if metric == "protein_aware_area_per_lipid":
                require(tool.get("name") in {"APL@Voro", "FATSLiM"}, "Protein-aware APL tool must be APL@Voro or FATSLiM")
                if tool.get("name") == "APL@Voro":
                    require(tool.get("version") == "3.3", "APL@Voro must be version 3.3")
                require(record.get("output_schema") == "per_saved_frame_protein_aware_popc_apl_v1", "Protein-aware APL output schema differs")
            else:
                require(tool.get("name") == "gorder" and tool.get("version"), "POPC S_CD tool must be a version-pinned gorder")
                require(record.get("output_schema") == "five_fixed_60ns_blocks_charmm36_popc_scd_v1", "POPC S_CD output schema differs")
            required_records = ("source_code_record", "version_capture_record", "command_record", "atom_mapping_record", "validation_report_record")
            require(all(isinstance(record.get(key), Mapping) for key in required_records), f"{metric} source-hashed validation records are incomplete")
            outputs = record.get("per_realization_outputs")
            require(isinstance(outputs, list) and [item.get("realization_id") for item in outputs] == list(REALIZATION_IDS), f"{metric} outputs must cover rep01-rep03 in order")
    gates = mapping.get("qc_gates", {})
    require(isinstance(gates.get("maximum_cumulative_leaflet_flip_events"), int) and int(gates["maximum_cumulative_leaflet_flip_events"]) >= 0, "qc_gates.maximum_cumulative_leaflet_flip_events is not frozen")
    for key in ("water_defect_largest_cluster_threshold", "water_defect_persistence_frames"):
        require(isinstance(gates.get(key), int) and int(gates[key]) >= 1, f"qc_gates.{key} must be a frozen positive integer")
    scd_validated = availability["popc_deuterium_order_parameters"].get("status") == "validated"
    for key in ("maximum_absolute_scd_adjacent_block_change", "maximum_absolute_scd_first_last_change"):
        value = gates.get(key)
        if scd_validated:
            parsed = float(value)
            require(math.isfinite(parsed) and 0.0 < parsed <= 1.5, f"qc_gates.{key} must be frozen in (0,1.5]")
        else:
            require(value is None, f"qc_gates.{key} must remain null while the gorder route is unavailable")
    settings = mapping.get("metric_settings", {})
    for key in ("phosphate_density_bandwidth_nm", "leaflet_hysteresis_nm", "hydrophobic_core_half_thickness_nm", "protein_xy_exclusion_nm", "water_cluster_cutoff_nm"):
        require(float(settings.get(key)) > 0.0, f"metric_settings.{key} must be positive")
    grid = np.asarray(settings.get("phosphate_density_grid_nm"), dtype=np.float64)
    require(grid.shape == (3,), "phosphate_density_grid_nm must be [minimum, maximum, points]")
    minimum, maximum, points_float = grid
    points = int(points_float)
    require(points == points_float and points >= 101 and minimum < 0.0 < maximum, "Phosphate density grid is invalid")
    require(points % 2 == 1 and math.isclose(abs(float(minimum)), abs(float(maximum)), rel_tol=0.0, abs_tol=1e-12), "Phosphate density grid must be odd and symmetric about zero")
    angle_tolerance = float(settings.get("orthorhombic_angle_tolerance_deg"))
    require(0.0 < angle_tolerance <= 0.1, "orthorhombic_angle_tolerance_deg must be frozen in (0, 0.1]")
    return {
        "upper": upper,
        "lower": lower,
        "tilt_upper": tilt_upper,
        "tilt_lower": tilt_lower,
        "protein": protein,
        "waters": waters,
        "availability": availability,
        "gates": gates,
        "settings": settings,
        "density_grid_nm": np.linspace(float(minimum), float(maximum), points),
        "topology_identity_sha256": fingerprint,
    }


def _read_external_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None, f"External membrane CSV lacks a header: {path}")
        rows = list(reader)
    require(rows, f"External membrane CSV is empty: {path}")
    return list(reader.fieldnames), rows


def _resolve_external_common(metric: str, record: Mapping[str, Any], mapping_base: Path) -> dict[str, Any]:
    resolved: dict[str, Any] = {}
    for key in ("source_code_record", "version_capture_record", "command_record", "atom_mapping_record", "validation_report_record"):
        path = resolve_record(mapping_base, record[key], f"external_metrics.{metric}.{key}")
        resolved[key] = {"path": str(path), "sha256": sha256_file(path)}
    return resolved


def _ingest_external_metrics(
    mapping: Mapping[str, Any],
    mapping_base: Path,
    realization_id: str,
    times_ns: np.ndarray,
    trajectory_sha256: str,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    results: dict[str, Any] = {}
    for metric, record in mapping["availability"].items():
        if record["status"] == "not_available":
            results[metric] = {"status": "not_available", "reason": record["reason"]}
            continue
        common = _resolve_external_common(metric, record, mapping_base)
        atom_mapping_payload = load_json(Path(common["atom_mapping_record"]["path"]))
        require(atom_mapping_payload.get("schema_version") == "1.0" and atom_mapping_payload.get("system_id") == "8kct_nirogacestat_native", f"{metric} atom-mapping schema/system differs")
        require(atom_mapping_payload.get("trajectory_atom_identity_sha256") == mapping["topology_identity_sha256"], f"{metric} atom mapping is bound to another topology")
        expected_scd_profile_keys: set[tuple[str, str]] | None = None
        if metric == "protein_aware_area_per_lipid":
            require(atom_mapping_payload.get("metric") == "protein_aware_area_per_lipid" and atom_mapping_payload.get("protein_footprint_included") is True, "APL mapping must explicitly include the protein footprint")
            require(atom_mapping_payload.get("lipid_resnames") == ["POPC"], "APL mapping must explicitly identify bulk POPC")
            for key in ("protein_atom_indices_sha256", "popc_atom_indices_sha256"):
                digest = str(atom_mapping_payload.get(key, ""))
                require(len(digest) == 64 and set(digest) <= set("0123456789abcdef"), f"APL mapping lacks {key}")
        else:
            require(atom_mapping_payload.get("metric") == "popc_deuterium_order_parameters", "gorder mapping metric differs")
            require(atom_mapping_payload.get("force_field_family") == "CHARMM36" and atom_mapping_payload.get("lipid_resname") == "POPC", "gorder mapping must identify CHARMM36 POPC")
            require(atom_mapping_payload.get("unsaturated_chain_geometry_explicit") is True, "gorder mapping must explicitly encode the unsaturated chain geometry")
            entries = atom_mapping_payload.get("carbon_hydrogen_mappings")
            require(isinstance(entries, list) and entries, "gorder carbon-hydrogen mapping is empty")
            expected_scd_profile_keys = set()
            for number, entry in enumerate(entries):
                require(isinstance(entry, Mapping), f"gorder mapping row {number} is invalid")
                key = (str(entry.get("chain_id", "")), str(entry.get("carbon_id", "")))
                require(all(key) and key not in expected_scd_profile_keys, f"gorder mapping row {number} has an empty/duplicate profile key")
                require(isinstance(entry.get("carbon_atom_name"), str) and entry["carbon_atom_name"], f"gorder mapping row {number} lacks a carbon atom name")
                hydrogens = entry.get("hydrogen_atom_names")
                require(isinstance(hydrogens, list) and hydrogens and all(isinstance(name, str) and name for name in hydrogens), f"gorder mapping row {number} lacks explicit hydrogen names")
                expected_scd_profile_keys.add(key)
        outputs = {item["realization_id"]: item["output"] for item in record["per_realization_outputs"]}
        output_path = resolve_record(mapping_base, outputs[realization_id], f"external_metrics.{metric}.{realization_id}.output")
        header, rows = _read_external_csv(output_path)
        tool_version = str(record["tool"]["version"])
        mapping_sha = common["atom_mapping_record"]["sha256"]
        required_metadata = {"realization_id", "source_trajectory_sha256", "atom_mapping_sha256", "tool_version"}
        require(required_metadata.issubset(header), f"{metric} output lacks provenance columns")
        for row in rows:
            require(row["realization_id"] == realization_id, f"{metric} output realization differs")
            require(row["source_trajectory_sha256"] == trajectory_sha256, f"{metric} output is bound to another trajectory")
            require(row["atom_mapping_sha256"] == mapping_sha, f"{metric} output is bound to another atom mapping")
            require(row["tool_version"] == tool_version, f"{metric} output tool version differs")
        if metric == "protein_aware_area_per_lipid":
            require({"time_ns", "protein_aware_popc_area_per_lipid_nm2"}.issubset(header), "APL output columns differ")
            external_times = np.asarray([float(row["time_ns"]) for row in rows], dtype=np.float64)
            values = np.asarray([float(row["protein_aware_popc_area_per_lipid_nm2"]) for row in rows], dtype=np.float64)
            require(np.all(np.isfinite(values)) and np.all(values > 0.0), "Protein-aware APL contains nonpositive/nonfinite values")
            require(len(external_times) == len(times_ns) and np.allclose(external_times, times_ns, rtol=0.0, atol=float(manifest["time_contract"]["endpoint_tolerance_ns"])), "Protein-aware APL must cover every saved 0-500 ns frame")
            mask = primary_window_mask(external_times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
            stationarity = stationarity_diagnostics(
                values[mask],
                external_times[mask],
                manifest["diagnostics"]["stationarity"],
                float(manifest["acceptance_gates"]["stationarity_scale_floors"]["membrane"]["protein_aware_area_per_lipid_nm2"]),
            )
            results[metric] = {
                "status": "validated",
                "qc_status": "pass" if stationarity["status"] == "pass" else "fail",
                "tool": dict(record["tool"]),
                "output": {"path": str(output_path), "sha256": sha256_file(output_path), "rows": len(rows)},
                "source_hashed_records": common,
                "primary_mean_nm2": float(np.mean(values[mask])),
                "stationarity": stationarity,
            }
        else:
            required = {"block_index_zero_based", "block_start_ns", "block_end_ns", "chain_id", "carbon_id", "s_cd"}
            require(required.issubset(header), "gorder S_CD output columns differ")
            profile_by_block: dict[int, dict[tuple[str, str], float]] = {}
            for row in rows:
                block = int(row["block_index_zero_based"])
                require(0 <= block < 5, "gorder S_CD block index is outside 0-4")
                require(math.isclose(float(row["block_start_ns"]), 200.0 + 60.0 * block, abs_tol=1e-9) and math.isclose(float(row["block_end_ns"]), 260.0 + 60.0 * block, abs_tol=1e-9), "gorder S_CD blocks must be the fixed five 60 ns intervals")
                key = (row["chain_id"], row["carbon_id"])
                require(key not in profile_by_block.setdefault(block, {}), "Duplicate gorder S_CD profile row")
                value = float(row["s_cd"])
                require(math.isfinite(value) and -0.5 <= value <= 1.0, "gorder S_CD value is outside the physical order-tensor range [-0.5,1]")
                profile_by_block[block][key] = value
            require(set(profile_by_block) == set(range(5)), "gorder S_CD output lacks a fixed block")
            profile_keys = set(profile_by_block[0])
            require(profile_keys and all(set(profile_by_block[index]) == profile_keys for index in range(1, 5)), "gorder S_CD atom/carbon profile differs among blocks")
            require(profile_keys == expected_scd_profile_keys, "gorder S_CD output does not exactly cover the frozen CHARMM36 carbon-hydrogen mapping")
            maximum_adjacent = max(abs(profile_by_block[index + 1][key] - profile_by_block[index][key]) for index in range(4) for key in profile_keys)
            maximum_first_last = max(abs(profile_by_block[4][key] - profile_by_block[0][key]) for key in profile_keys)
            gates = mapping["gates"]
            scd_pass = maximum_adjacent <= float(gates["maximum_absolute_scd_adjacent_block_change"]) and maximum_first_last <= float(gates["maximum_absolute_scd_first_last_change"])
            results[metric] = {
                "status": "validated",
                "qc_status": "pass" if scd_pass else "fail",
                "tool": dict(record["tool"]),
                "output": {"path": str(output_path), "sha256": sha256_file(output_path), "rows": len(rows)},
                "source_hashed_records": common,
                "profile_key_count": len(profile_keys),
                "maximum_absolute_adjacent_block_change": float(maximum_adjacent),
                "maximum_absolute_first_last_change": float(maximum_first_last),
                "gates": {
                    "maximum_absolute_scd_adjacent_block_change": float(gates["maximum_absolute_scd_adjacent_block_change"]),
                    "maximum_absolute_scd_first_last_change": float(gates["maximum_absolute_scd_first_last_change"]),
                },
            }
    return results


def _process_realization(
    manifest: Mapping[str, Any],
    manifest_base: Path,
    realization: Mapping[str, Any],
    mapping_record: Mapping[str, Any],
    mapping_base: Path,
    output_directory: Path,
) -> tuple[dict[str, Any], np.ndarray]:
    realization_id = str(realization["realization_id"])
    topology_path = resolve_record(manifest_base, realization["topology"], f"{realization_id}.topology")
    trajectory_path = resolve_record(manifest_base, realization["centered_system_trajectory"], f"{realization_id}.centered_system_trajectory")
    try:
        universe = mda.Universe(str(topology_path), str(trajectory_path))
    except Exception as exc:
        raise ContractError(f"MDAnalysis could not open {realization_id} for membrane QC: {exc}") from exc
    mapping = _validate_mapping(mapping_record, universe)
    times_ns = np.asarray([float(ts.time) / 1000.0 for ts in universe.trajectory], dtype=np.float64)
    tolerance_ns = float(manifest["time_contract"]["endpoint_tolerance_ns"])
    saved_step_ns = validate_time_axis(times_ns, tolerance_ns)
    window_mask = primary_window_mask(times_ns, tolerance_ns)
    trajectory_sha256 = sha256_file(trajectory_path)
    external_metrics = _ingest_external_metrics(mapping, mapping_base, realization_id, times_ns, trajectory_sha256, manifest)

    grid_nm = mapping["density_grid_nm"]
    bandwidth_nm = float(mapping["settings"]["phosphate_density_bandwidth_nm"])
    hysteresis_nm = float(mapping["settings"]["leaflet_hysteresis_nm"])
    core_half_thickness_nm = float(mapping["settings"]["hydrophobic_core_half_thickness_nm"])
    protein_exclusion_nm = float(mapping["settings"]["protein_xy_exclusion_nm"])
    cluster_cutoff_nm = float(mapping["settings"]["water_cluster_cutoff_nm"])
    angle_tolerance = float(mapping["settings"]["orthorhombic_angle_tolerance_deg"])
    upper_state = np.ones(len(mapping["upper"]), dtype=np.int8)
    lower_state = -np.ones(len(mapping["lower"]), dtype=np.int8)
    cumulative_flips = 0
    leaflet_events: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    primary_upper_density = np.zeros(len(grid_nm), dtype=np.float64)
    primary_lower_density = np.zeros(len(grid_nm), dtype=np.float64)
    primary_density_frames = 0

    for frame_index, ts in enumerate(universe.trajectory):
        coordinates = np.asarray(universe.atoms.positions, dtype=np.float64)
        require(np.all(np.isfinite(coordinates)), f"{realization_id} frame {frame_index} contains nonfinite coordinates")
        dimensions = np.asarray(ts.dimensions, dtype=np.float64)
        require(dimensions.shape == (6,) and np.all(np.isfinite(dimensions)) and np.all(dimensions[:3] > 0.0), f"{realization_id} frame {frame_index} has an invalid box")
        require(np.all(np.abs(dimensions[3:] - 90.0) <= angle_tolerance), f"{realization_id} frame {frame_index} is not orthorhombic; protein-aware triclinic water-defect logic is not validated")
        vectors = np.asarray(ts.triclinic_dimensions, dtype=np.float64)
        require(vectors.shape == (3, 3) and np.all(np.isfinite(vectors)), f"{realization_id} frame {frame_index} lacks finite triclinic vectors")
        a, b, c = vectors
        lateral_area_nm2 = float(np.linalg.norm(np.cross(a, b)) / 100.0)
        c_length_nm = float(np.linalg.norm(c) / 10.0)
        volume_nm3 = float(abs(np.linalg.det(vectors)) / 1000.0)
        require(lateral_area_nm2 > 0.0 and c_length_nm > 0.0 and volume_nm3 > 0.0, f"{realization_id} frame {frame_index} has a nonpositive cell metric")
        lx_nm, ly_nm = float(dimensions[0] / 10.0), float(dimensions[1] / 10.0)

        upper_relative, lower_relative, phosphate_midplane_angstrom = leaflet_relative_z_nm(
            coordinates[mapping["upper"], 2],
            coordinates[mapping["lower"], 2],
            float(dimensions[2]),
        )
        upper_density = gaussian_density(upper_relative, grid_nm, bandwidth_nm)
        lower_density = gaussian_density(lower_relative, grid_nm, bandwidth_nm)
        upper_peak_index = int(np.argmax(upper_density))
        lower_peak_index = int(np.argmax(lower_density))
        require(0 < upper_peak_index < len(grid_nm) - 1 and 0 < lower_peak_index < len(grid_nm) - 1, f"{realization_id} frame {frame_index} phosphate-density peak reached a frozen grid boundary")
        upper_peak = float(grid_nm[upper_peak_index])
        lower_peak = float(grid_nm[lower_peak_index])
        thickness_nm = upper_peak - lower_peak
        require(thickness_nm > 0.0, f"{realization_id} frame {frame_index} has inverted/nonpositive phosphate peak separation")
        if frame_index == 0:
            require(np.all(upper_relative > hysteresis_nm) and np.all(lower_relative < -hysteresis_nm), "Frozen leaflet identities do not match the first centered frame")

        new_flips = 0
        for group_name, indices, values, states, expected_state in (
            ("upper", mapping["upper"], upper_relative, upper_state, 1),
            ("lower", mapping["lower"], lower_relative, lower_state, -1),
        ):
            for local_index, value in enumerate(values):
                new_state = int(states[local_index])
                if value > hysteresis_nm:
                    new_state = 1
                elif value < -hysteresis_nm:
                    new_state = -1
                if new_state != int(states[local_index]):
                    new_flips += 1
                    cumulative_flips += 1
                    atom_index = int(indices[local_index])
                    leaflet_events.append({
                        "system_id": manifest["system_id"],
                        "realization_id": realization_id,
                        "frame_index_zero_based": frame_index,
                        "time_ns": float(times_ns[frame_index]),
                        "frozen_leaflet": group_name,
                        "trajectory_atom_index_zero_based": atom_index,
                        "old_state": int(states[local_index]),
                        "new_state": new_state,
                        "relative_z_nm": float(value),
                        "event_retained": True,
                    })
                    states[local_index] = new_state
        upper_mismatch = int(np.count_nonzero(upper_state != 1))
        lower_mismatch = int(np.count_nonzero(lower_state != -1))

        normal = np.cross(a, b)
        normal /= np.linalg.norm(normal)
        upper_anchor = np.mean(coordinates[mapping["tilt_upper"]], axis=0)
        lower_anchor = np.mean(coordinates[mapping["tilt_lower"]], axis=0)
        axis = upper_anchor - lower_anchor
        axis -= dimensions[:3] * np.rint(axis / dimensions[:3])
        require(np.linalg.norm(axis) > 0.0, f"{realization_id} frame {frame_index} has a zero-length protein tilt axis")
        cosine = float(np.clip(abs(np.dot(axis, normal)) / np.linalg.norm(axis), 0.0, 1.0))
        tilt_deg = float(np.degrees(np.arccos(cosine)))

        water_relative = _relative_to_midplane_nm(
            coordinates[mapping["waters"], 2],
            phosphate_midplane_angstrom,
            float(dimensions[2]),
        )
        protein_relative = _relative_to_midplane_nm(
            coordinates[mapping["protein"], 2],
            phosphate_midplane_angstrom,
            float(dimensions[2]),
        )
        water_count, largest_cluster = _outside_protein_water(
            coordinates[mapping["waters"]],
            water_relative,
            coordinates[mapping["protein"]],
            protein_relative,
            lx_nm,
            ly_nm,
            float(dimensions[2] / 10.0),
            core_half_thickness_nm,
            protein_exclusion_nm,
            cluster_cutoff_nm,
        )
        row = {
            "system_id": manifest["system_id"],
            "realization_id": realization_id,
            "frame_index_zero_based": frame_index,
            "time_ns": float(times_ns[frame_index]),
            "in_primary_window_200_500_ns": int(bool(window_mask[frame_index])),
            "phosphate_peak_thickness_nm": thickness_nm,
            "upper_phosphate_peak_z_relative_nm": upper_peak,
            "lower_phosphate_peak_z_relative_nm": lower_peak,
            "upper_leaflet_mismatch_count": upper_mismatch,
            "lower_leaflet_mismatch_count": lower_mismatch,
            "new_leaflet_flip_events": new_flips,
            "cumulative_leaflet_flip_events": cumulative_flips,
            "cell_lateral_area_nm2_not_apl": lateral_area_nm2,
            "box_z_vector_length_nm": c_length_nm,
            "cell_volume_nm3": volume_nm3,
            "protein_tilt_deg": tilt_deg,
            "outside_protein_hydrophobic_core_water_count": water_count,
            "outside_protein_hydrophobic_core_largest_water_cluster": largest_cluster,
        }
        require(all(math.isfinite(float(value)) for key, value in row.items() if key not in {"system_id", "realization_id"}), f"{realization_id} frame {frame_index} produced a nonfinite membrane metric")
        raw_rows.append(row)
        if window_mask[frame_index]:
            primary_upper_density += upper_density
            primary_lower_density += lower_density
            primary_density_frames += 1

    require(len(raw_rows) == len(times_ns), f"{realization_id} lost frames during membrane QC")
    require(primary_density_frames == int(np.count_nonzero(window_mask)), f"{realization_id} primary density frame accounting differs")
    primary_upper_density /= primary_density_frames
    primary_lower_density /= primary_density_frames
    profile_rows = [
        {
            "system_id": manifest["system_id"],
            "realization_id": realization_id,
            "window_start_ns": PRIMARY_WINDOW_NS[0],
            "window_end_ns": PRIMARY_WINDOW_NS[1],
            "z_relative_nm": float(position),
            "upper_leaflet_phosphate_density_arbitrary": float(upper_density_value),
            "lower_leaflet_phosphate_density_arbitrary": float(lower_density_value),
        }
        for position, upper_density_value, lower_density_value in zip(grid_nm, primary_upper_density, primary_lower_density, strict=True)
    ]

    realization_directory = output_directory / realization_id
    raw_path = realization_directory / "membrane_raw_unsmoothed.csv"
    atomic_write_csv(raw_path, RAW_FIELDS, raw_rows)
    profile_path = realization_directory / "phosphate_density_profile_200_500ns.csv"
    atomic_write_csv(profile_path, list(profile_rows[0]), profile_rows)
    events_path = realization_directory / "leaflet_events.csv"
    event_fields = ["system_id", "realization_id", "frame_index_zero_based", "time_ns", "frozen_leaflet", "trajectory_atom_index_zero_based", "old_state", "new_state", "relative_z_nm", "event_retained"]
    atomic_write_csv(events_path, event_fields, leaflet_events)

    scalar_metrics = [field for field in RAW_FIELDS if field not in {"system_id", "realization_id", "frame_index_zero_based", "time_ns", "in_primary_window_200_500_ns"}]
    diagnostics_payload: dict[str, Any] = {}
    block_rows: list[dict[str, Any]] = []
    review_rows: list[dict[str, Any]] = []
    for metric in scalar_metrics:
        values = np.asarray([float(row[metric]) for row in raw_rows], dtype=np.float64)
        diagnostics_payload[metric] = {}
        for window_name, mask in (("full_0_500_ns", np.ones(len(times_ns), dtype=bool)), ("primary_200_500_ns", window_mask)):
            diagnostic, blocks = block_diagnostics(values[mask], times_ns[mask], manifest["diagnostics"])
            if window_name == "primary_200_500_ns" and metric in MEMBRANE_STATIONARITY_METRICS:
                diagnostic["stationarity"] = stationarity_diagnostics(
                    values[mask],
                    times_ns[mask],
                    manifest["diagnostics"]["stationarity"],
                    float(manifest["acceptance_gates"]["stationarity_scale_floors"]["membrane"][metric]),
                )
            diagnostics_payload[metric][window_name] = diagnostic
            for block in blocks:
                block_rows.append({"system_id": manifest["system_id"], "realization_id": realization_id, "metric": metric, "window": window_name, **block})
        first_difference_summary, flags = robust_first_difference(values, times_ns, float(manifest["diagnostics"]["robust_first_difference_z_threshold"]))
        diagnostics_payload[metric]["first_difference_review"] = first_difference_summary
        for flag in flags:
            review_rows.append({"system_id": manifest["system_id"], "realization_id": realization_id, "metric": metric, **flag})
    block_path = realization_directory / "membrane_block_summaries.csv"
    block_fields = ["system_id", "realization_id", "metric", "window", "block_index_zero_based", "start_time_ns", "end_time_ns", "frame_count", "mean", "median", "minimum", "maximum"]
    atomic_write_csv(block_path, block_fields, block_rows)
    review_path = realization_directory / "membrane_first_difference_review_flags.csv"
    review_fields = ["system_id", "realization_id", "metric", "row_index_zero_based", "time_before_ns", "time_after_ns", "value_before", "value_after", "first_difference", "median_first_difference", "mad_first_difference", "robust_z", "method", "review_required", "point_retained"]
    atomic_write_csv(review_path, review_fields, review_rows)
    availability_path = realization_directory / "membrane_metric_availability.json"
    atomic_write_json(availability_path, {
        "schema_version": "1.0",
        "system_id": manifest["system_id"],
        "realization_id": realization_id,
        "cell_lateral_area": {"status": "available", "label": "periodic cell lateral area; not area per lipid", "units": "nm^2"},
        **external_metrics,
    })

    mandatory_metrics = [
        "phosphate_peak_thickness_nm",
        "cell_lateral_area_nm2_not_apl",
        "box_z_vector_length_nm",
        "cell_volume_nm3",
        "protein_tilt_deg",
    ]
    sampling_failures = []
    stationarity_failures = []
    for metric in mandatory_metrics:
        primary_diagnostic = diagnostics_payload[metric]["primary_200_500_ns"]
        if primary_diagnostic["status"] != "pass":
            sampling_failures.append(metric)
        if primary_diagnostic["stationarity"]["status"] != "pass":
            stationarity_failures.append(metric)
    external_unavailable = [metric for metric, record in external_metrics.items() if record.get("status") != "validated"]
    external_qc_failures = [metric for metric, record in external_metrics.items() if record.get("status") == "validated" and record.get("qc_status") != "pass"]
    leaflet_failed = cumulative_flips > int(mapping["gates"]["maximum_cumulative_leaflet_flip_events"])
    cluster_threshold = int(mapping["gates"]["water_defect_largest_cluster_threshold"])
    persistence_required = int(mapping["gates"]["water_defect_persistence_frames"])
    consecutive = 0
    maximum_consecutive = 0
    for row in raw_rows:
        if int(row["outside_protein_hydrophobic_core_largest_water_cluster"]) >= cluster_threshold:
            consecutive += 1
            maximum_consecutive = max(maximum_consecutive, consecutive)
        else:
            consecutive = 0
    water_defect_failed = cluster_threshold > 0 and maximum_consecutive >= persistence_required
    qc_failures = []
    if leaflet_failed:
        qc_failures.append("leaflet_flip_gate")
    if water_defect_failed:
        qc_failures.append("persistent_outside_protein_water_defect_gate")
    technical_status = "fail" if qc_failures else "pass"
    sampling_status = "inconclusive" if sampling_failures or stationarity_failures or external_qc_failures else "pass"
    summary = {
        "schema_version": "1.0",
        "analysis": "membrane_qc_mdanalysis",
        "system_id": manifest["system_id"],
        "realization_id": realization_id,
        "construction_count": 1,
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "input_frame_count": len(times_ns),
        "primary_window_frame_count": int(np.count_nonzero(window_mask)),
        "saved_step_ns": saved_step_ns,
        "first_time_ns": float(times_ns[0]),
        "last_time_ns": float(times_ns[-1]),
        "technical_status": technical_status,
        "sampling_status": sampling_status,
        "qc_failures": qc_failures,
        "sampling_failures": sampling_failures,
        "stationarity_status": "pass" if not stationarity_failures else "fail",
        "stationarity_failures": stationarity_failures,
        "external_membrane_metric_qc_failures": external_qc_failures,
        "cumulative_leaflet_flip_events": cumulative_flips,
        "maximum_consecutive_water_defect_frames": maximum_consecutive,
        "first_difference_review_flag_count": len(review_rows),
        "first_difference_flags_remove_points": False,
        "cell_lateral_area_is_area_per_lipid": False,
        "preproduction_status": "pass" if not external_unavailable and not external_qc_failures else "blocked_external_membrane_metrics",
        "preproduction_no_go_requirements": sorted(set(external_unavailable + external_qc_failures)),
        "mdanalysis_version": mda.__version__,
        "trajectory_atom_identity_sha256": mapping["topology_identity_sha256"],
        "input_files": {
            "topology": {"path": str(topology_path), "sha256": sha256_file(topology_path)},
            "centered_system_trajectory": {"path": str(trajectory_path), "sha256": sha256_file(trajectory_path)},
        },
        "data_handling": {"raw_frames_retained": True, "smoothing": False, "frame_deletion": False, "interpolation": False, "realization_pooling": False},
        "diagnostics": diagnostics_payload,
        "outputs": {},
    }
    summary_path = realization_directory / "membrane_summary.json"
    for path in (raw_path, profile_path, events_path, block_path, review_path, availability_path):
        summary["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    atomic_write_json(summary_path, summary)
    return summary, times_ns


def run(manifest_path: Path, output_root: Path, allow_synthetic: bool = False) -> dict[str, Any]:
    check_mdanalysis_version(mda.__version__)
    manifest, manifest_base = validate_primary_manifest(manifest_path, allow_synthetic=allow_synthetic)
    mapping_path = resolve_record(manifest_base, manifest["mapping_records"]["membrane"], "mapping_records.membrane")
    mapping_record = load_json(mapping_path)
    output_directory = output_root.resolve() / "membrane_qc"
    require(not output_directory.exists(), f"Refusing to overwrite an existing membrane output directory: {output_directory}")
    output_directory.mkdir(parents=True)
    summaries = []
    shared_times: np.ndarray | None = None
    for realization in manifest["realizations"]:
        summary, times = _process_realization(manifest, manifest_base, realization, mapping_record, mapping_path.parent, output_directory)
        if shared_times is None:
            shared_times = times
        else:
            require(np.allclose(times, shared_times, rtol=0.0, atol=float(manifest["time_contract"]["endpoint_tolerance_ns"])), "Saved frame times differ among realizations")
        summaries.append(summary)
    require([summary["realization_id"] for summary in summaries] == list(REALIZATION_IDS), "Membrane QC lost or reordered a realization")
    technical_status = "pass" if all(summary["technical_status"] == "pass" for summary in summaries) else "fail"
    sampling_status = "pass" if all(summary["sampling_status"] == "pass" for summary in summaries) else "inconclusive"
    preproduction_no_go = sorted({metric for summary in summaries for metric in summary["preproduction_no_go_requirements"]})
    overall_status = "pass" if technical_status == "pass" and sampling_status == "pass" and not preproduction_no_go else "inconclusive"
    complete = {
        "schema_version": "1.0",
        "status": overall_status,
        "technical_status": technical_status,
        "sampling_status": sampling_status,
        "system_id": manifest["system_id"],
        "construction_count": 1,
        "realization_ids": list(REALIZATION_IDS),
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "extension_or_recovery_window": False,
        "mdanalysis_version": mda.__version__,
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path.resolve())},
        "mapping_record": {"path": str(mapping_path), "sha256": sha256_file(mapping_path)},
        "cell_lateral_area_is_area_per_lipid": False,
        "preproduction_status": "pass" if not preproduction_no_go else "blocked_external_membrane_metrics",
        "preproduction_no_go_requirements": preproduction_no_go,
        "realization_summaries": [
            {"realization_id": item["realization_id"], "technical_status": item["technical_status"], "sampling_status": item["sampling_status"], "stationarity_status": item["stationarity_status"], "qc_failures": item["qc_failures"], "preproduction_no_go_requirements": item["preproduction_no_go_requirements"]}
            for item in summaries
        ],
    }
    atomic_write_json(output_directory / "COMPLETE.json", complete)
    return complete


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True, help="Approved primary postprocessing manifest")
    parser.add_argument("--output-root", type=Path, required=True, help="New output root; existing membrane_qc is never overwritten")
    parser.add_argument("--allow-synthetic", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = run(args.manifest, args.output_root, allow_synthetic=args.allow_synthetic)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
