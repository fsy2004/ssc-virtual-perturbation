#!/usr/bin/env python3
"""Independently validate the three-realization 500 ns primary MD outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from gmx_energy_qc import REQUIRED_TERM_KEYS
from primary_postprocessing_common import (
    PRIMARY_WINDOW_NS,
    REALIZATION_IDS,
    ContractError,
    atomic_write_json,
    continuous_true_events,
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
)


STRUCTURAL_REQUIRED = (
    "pocket_aligned_o6u_heavy_rmsd_nm",
    "pocket_aligned_o6u_com_displacement_nm",
    "tm_core_ca_rmsd_nm",
    "protein_ca_rmsd_nm",
    "native_contact_fraction",
)

MEMBRANE_REQUIRED = (
    "phosphate_peak_thickness_nm",
    "cell_lateral_area_nm2_not_apl",
    "box_z_vector_length_nm",
    "cell_volume_nm3",
    "protein_tilt_deg",
    "outside_protein_hydrophobic_core_water_count",
    "outside_protein_hydrophobic_core_largest_water_cluster",
)


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file() and path.stat().st_size > 0, f"Missing/empty CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None, f"CSV has no header: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def finite_column(rows: Iterable[Mapping[str, str]], name: str) -> np.ndarray:
    values = []
    for number, row in enumerate(rows):
        require(name in row, f"CSV lacks {name}")
        try:
            value = float(row[name])
        except (TypeError, ValueError) as exc:
            raise ContractError(f"{name} row {number} is nonnumeric") from exc
        require(math.isfinite(value), f"{name} row {number} is NaN or infinity")
        values.append(value)
    return np.asarray(values, dtype=np.float64)


def verify_external_membrane_metric(
    record: Mapping[str, Any],
    metric: str,
    realization_id: str,
    expected_times: np.ndarray,
    expected_trajectory_sha256: str,
    expected_topology_identity_sha256: str,
    manifest: Mapping[str, Any],
) -> list[str]:
    if record.get("status") == "not_available":
        require(isinstance(record.get("reason"), str) and len(record["reason"]) >= 20, f"{metric} lacks a concrete NO-GO reason")
        return [f"preproduction_no_go:{metric}:external_build_specific_validation_required"]
    require(record.get("status") == "validated", f"Unexpected external membrane metric status: {metric}")
    require(record.get("qc_status") in {"pass", "fail"}, f"{metric} lacks a QC status")
    output = record.get("output", {})
    output_path = Path(str(output.get("path", "")))
    require(output_path.is_file() and output.get("sha256") == sha256_file(output_path), f"{metric} output is missing or changed: {realization_id}")
    source_records = record.get("source_hashed_records", {})
    expected_source_keys = {"source_code_record", "version_capture_record", "command_record", "atom_mapping_record", "validation_report_record"}
    require(set(source_records) == expected_source_keys, f"{metric} source-hashed records differ: {realization_id}")
    for label, source in source_records.items():
        path = Path(str(source.get("path", "")))
        require(path.is_file() and source.get("sha256") == sha256_file(path), f"{metric} source record is missing or changed: {realization_id}/{label}")
    header, rows = read_csv(output_path)
    provenance_columns = {"realization_id", "source_trajectory_sha256", "atom_mapping_sha256", "tool_version"}
    require(provenance_columns.issubset(header), f"{metric} output lacks provenance columns")
    expected_mapping_sha = source_records["atom_mapping_record"]["sha256"]
    atom_mapping_payload = load_json(Path(str(source_records["atom_mapping_record"]["path"])))
    require(atom_mapping_payload.get("schema_version") == "1.0" and atom_mapping_payload.get("system_id") == "8kct_nirogacestat_native", f"{metric} atom mapping schema/system differs")
    require(atom_mapping_payload.get("trajectory_atom_identity_sha256") == expected_topology_identity_sha256, f"{metric} atom mapping is bound to another topology")
    expected_tool_version = str(record.get("tool", {}).get("version", ""))
    for row in rows:
        require(row["realization_id"] == realization_id, f"{metric} output realization differs")
        require(row["source_trajectory_sha256"] == expected_trajectory_sha256, f"{metric} output is bound to another trajectory")
        require(row["atom_mapping_sha256"] == expected_mapping_sha, f"{metric} output is bound to another atom mapping")
        require(row["tool_version"] == expected_tool_version, f"{metric} output tool version differs")
    if metric == "protein_aware_area_per_lipid":
        require(atom_mapping_payload.get("metric") == "protein_aware_area_per_lipid" and atom_mapping_payload.get("protein_footprint_included") is True and atom_mapping_payload.get("lipid_resnames") == ["POPC"], "Protein-aware APL atom mapping differs")
        for key in ("protein_atom_indices_sha256", "popc_atom_indices_sha256"):
            digest = str(atom_mapping_payload.get(key, ""))
            require(len(digest) == 64 and set(digest) <= set("0123456789abcdef"), f"Protein-aware APL mapping lacks {key}")
        require({"time_ns", "protein_aware_popc_area_per_lipid_nm2"}.issubset(header), "Protein-aware APL output columns differ")
        times = finite_column(rows, "time_ns")
        require(len(times) == len(expected_times) and np.array_equal(times, expected_times), f"Protein-aware APL does not cover every saved frame: {realization_id}")
        values = finite_column(rows, "protein_aware_popc_area_per_lipid_nm2")
        require(np.all(values > 0.0), f"Protein-aware APL contains nonpositive values: {realization_id}")
        mask = primary_window_mask(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
        recomputed = stationarity_diagnostics(
            values[mask],
            times[mask],
            manifest["diagnostics"]["stationarity"],
            float(manifest["acceptance_gates"]["stationarity_scale_floors"]["membrane"]["protein_aware_area_per_lipid_nm2"]),
        )
        verify_stationarity_report(record.get("stationarity", {}), recomputed, f"{realization_id}/membrane/protein_aware_area_per_lipid")
        expected_qc = "pass" if recomputed["status"] == "pass" else "fail"
    else:
        require(atom_mapping_payload.get("metric") == "popc_deuterium_order_parameters" and atom_mapping_payload.get("force_field_family") == "CHARMM36" and atom_mapping_payload.get("lipid_resname") == "POPC" and atom_mapping_payload.get("unsaturated_chain_geometry_explicit") is True, "gorder atom mapping differs")
        mapping_entries = atom_mapping_payload.get("carbon_hydrogen_mappings")
        require(isinstance(mapping_entries, list) and mapping_entries, "gorder atom mapping is empty")
        for number, entry in enumerate(mapping_entries):
            require(isinstance(entry.get("carbon_atom_name"), str) and entry["carbon_atom_name"], f"gorder mapping row {number} lacks a carbon atom name")
            hydrogens = entry.get("hydrogen_atom_names")
            require(isinstance(hydrogens, list) and hydrogens and all(isinstance(name, str) and name for name in hydrogens), f"gorder mapping row {number} lacks explicit hydrogen names")
        expected_profile_keys = {(str(entry.get("chain_id", "")), str(entry.get("carbon_id", ""))) for entry in mapping_entries}
        require(len(expected_profile_keys) == len(mapping_entries) and all(all(key) for key in expected_profile_keys), "gorder atom mapping has empty/duplicate profile keys")
        required = {"block_index_zero_based", "block_start_ns", "block_end_ns", "chain_id", "carbon_id", "s_cd"}
        require(required.issubset(header), "gorder S_CD output columns differ")
        profiles: dict[int, dict[tuple[str, str], float]] = {}
        for row in rows:
            block = int(row["block_index_zero_based"])
            require(0 <= block < 5, f"gorder block index is invalid: {realization_id}")
            require(math.isclose(float(row["block_start_ns"]), 200.0 + 60.0 * block, rel_tol=0.0, abs_tol=1e-9) and math.isclose(float(row["block_end_ns"]), 260.0 + 60.0 * block, rel_tol=0.0, abs_tol=1e-9), f"gorder fixed block boundaries differ: {realization_id}/{block}")
            key = (row["chain_id"], row["carbon_id"])
            require(key not in profiles.setdefault(block, {}), f"Duplicate gorder row: {realization_id}/{block}/{key}")
            value = float(row["s_cd"])
            require(math.isfinite(value) and -0.5 <= value <= 1.0, f"gorder S_CD is outside [-0.5,1]: {realization_id}")
            profiles[block][key] = value
        require(set(profiles) == set(range(5)), f"gorder S_CD lacks a fixed block: {realization_id}")
        keys = set(profiles[0])
        require(keys and all(set(profiles[index]) == keys for index in range(1, 5)), f"gorder profile keys differ among blocks: {realization_id}")
        require(keys == expected_profile_keys, f"gorder output does not exactly cover the frozen mapping: {realization_id}")
        adjacent = max(abs(profiles[index + 1][key] - profiles[index][key]) for index in range(4) for key in keys)
        first_last = max(abs(profiles[4][key] - profiles[0][key]) for key in keys)
        gates = record.get("gates", {})
        require(math.isclose(float(record.get("maximum_absolute_adjacent_block_change")), adjacent, rel_tol=0.0, abs_tol=1e-12), f"gorder adjacent-block result differs: {realization_id}")
        require(math.isclose(float(record.get("maximum_absolute_first_last_change")), first_last, rel_tol=0.0, abs_tol=1e-12), f"gorder first-last result differs: {realization_id}")
        expected_qc = "pass" if adjacent <= float(gates["maximum_absolute_scd_adjacent_block_change"]) and first_last <= float(gates["maximum_absolute_scd_first_last_change"]) else "fail"
    require(record.get("qc_status") == expected_qc, f"{metric} QC status differs: {realization_id}")
    return [] if expected_qc == "pass" else [f"{realization_id}:membrane:{metric}:external_metric_qc_failed"]


def verify_summary_hashes(summary: Mapping[str, Any], directory: Path) -> None:
    outputs = summary.get("outputs")
    require(isinstance(outputs, Mapping) and outputs, f"Output hashes are missing in {directory}")
    for filename, record in outputs.items():
        path = directory / filename
        require(path.is_file(), f"Hashed output is missing: {path}")
        require(int(record.get("bytes", -1)) == path.stat().st_size, f"Byte count differs: {path}")
        require(record.get("sha256") == sha256_file(path), f"SHA-256 differs: {path}")


def verify_stationarity_report(reported: Mapping[str, Any], recomputed: Mapping[str, Any], label: str) -> None:
    for key in ("method", "window_ns", "fixed_time_blocks", "minimum_frames_per_block", "failed_gates", "status", "all_input_points_retained"):
        require(reported.get(key) == recomputed.get(key), f"Stationarity {key} differs: {label}")
    for key in ("metric_scale_floor", "raw_median", "raw_mad", "robust_scale", "linear_slope_per_ns", "linear_change_over_window"):
        require(math.isclose(float(reported.get(key)), float(recomputed.get(key)), rel_tol=1e-12, abs_tol=1e-12), f"Stationarity {key} differs: {label}")
    for group in ("observed", "limits"):
        require(set(reported.get(group, {})) == set(recomputed.get(group, {})), f"Stationarity {group} keys differ: {label}")
        for key, value in recomputed[group].items():
            require(math.isclose(float(reported[group][key]), float(value), rel_tol=1e-12, abs_tol=1e-12), f"Stationarity {group}.{key} differs: {label}")
    reported_blocks = reported.get("blocks")
    recomputed_blocks = recomputed.get("blocks")
    require(isinstance(reported_blocks, list) and len(reported_blocks) == len(recomputed_blocks), f"Stationarity block count differs: {label}")
    for observed, expected in zip(reported_blocks, recomputed_blocks, strict=True):
        for key in ("block_index_zero_based", "frame_count"):
            require(int(observed[key]) == int(expected[key]), f"Stationarity block {key} differs: {label}")
        for key in ("start_time_ns", "end_time_ns", "median", "mean"):
            require(math.isclose(float(observed[key]), float(expected[key]), rel_tol=1e-12, abs_tol=1e-12), f"Stationarity block {key} differs: {label}")


def verify_complete(path: Path, manifest_path: Path, allowed_statuses: set[str]) -> dict[str, Any]:
    complete = load_json(path)
    require(complete.get("schema_version") == "1.0", f"Unexpected COMPLETE schema: {path}")
    require(complete.get("status") in allowed_statuses, f"Unexpected COMPLETE status: {path}")
    require(complete.get("system_id") == "8kct_nirogacestat_native", f"Wrong system: {path}")
    require(complete.get("realization_ids") == list(REALIZATION_IDS), f"COMPLETE does not contain rep01-rep03: {path}")
    require(float(complete.get("production_duration_ns")) == 500.0, f"COMPLETE duration is not 500 ns: {path}")
    require(complete.get("primary_analysis_window_ns") == list(PRIMARY_WINDOW_NS), f"COMPLETE window is not 200-500 ns: {path}")
    require(complete.get("extension_or_recovery_window") is False, f"COMPLETE contains recovery/extension logic: {path}")
    require(complete.get("manifest", {}).get("sha256") == sha256_file(manifest_path), f"Manifest hash differs in {path}")
    return complete


def validate_structural(output_root: Path, manifest_path: Path, manifest: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, np.ndarray], list[str], list[dict[str, str]]]:
    root = output_root / "structural_analysis"
    complete = verify_complete(root / "COMPLETE.json", manifest_path, {"pass", "inconclusive"})
    time_arrays: dict[str, np.ndarray] = {}
    reasons: list[str] = []
    all_flags: list[dict[str, str]] = []
    for realization_id in REALIZATION_IDS:
        directory = root / realization_id
        summary = load_json(directory / "structural_summary.json")
        require(summary.get("realization_id") == realization_id, f"Structural summary ID differs: {realization_id}")
        require(summary.get("construction_count") == 1, "Structural summary construction_count must be 1")
        require(summary.get("data_handling") == {"raw_frames_retained": True, "smoothing": False, "frame_deletion": False, "interpolation": False, "realization_pooling": False}, f"Structural data-handling contract differs: {realization_id}")
        require(summary.get("first_difference_flags_remove_points") is False, f"Structural flags remove points: {realization_id}")
        verify_summary_hashes(summary, directory)
        header, rows = read_csv(directory / "structural_raw_unsmoothed.csv")
        require(rows, f"Structural raw CSV is empty: {realization_id}")
        require(all(name in header for name in STRUCTURAL_REQUIRED), f"Structural raw CSV lacks primary metrics: {realization_id}")
        require(len(rows) == int(summary["input_frame_count"]), f"Structural row count differs: {realization_id}")
        indices = finite_column(rows, "frame_index_zero_based")
        require(np.array_equal(indices, np.arange(len(rows), dtype=np.float64)), f"Structural frame indices are not complete: {realization_id}")
        times = finite_column(rows, "time_ns")
        validate_time_axis(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
        window = primary_window_mask(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
        flags = finite_column(rows, "in_primary_window_200_500_ns").astype(np.int64)
        require(np.array_equal(flags, window.astype(np.int64)), f"Structural window flags differ: {realization_id}")
        require(int(np.count_nonzero(window)) == int(summary["primary_window_frame_count"]), f"Structural primary-window row count differs: {realization_id}")
        for metric in STRUCTURAL_REQUIRED:
            values = finite_column(rows, metric)
            if metric == "native_contact_fraction":
                require(np.all((0.0 <= values) & (values <= 1.0)), f"Native-contact fraction is outside [0,1]: {realization_id}")
            else:
                require(np.all(values >= 0.0), f"Negative structural distance/RMSD: {realization_id}/{metric}")
            diagnostic = summary.get("diagnostics", {}).get(metric, {}).get("primary_200_500_ns", {})
            require(diagnostic.get("minimum_complete_blocks", 0) >= 5, f"Structural minimum block gate was weakened: {realization_id}/{metric}")
            require(float(diagnostic.get("block_tau_multiplier", 0.0)) >= 10.0, f"Structural block multiplier was weakened: {realization_id}/{metric}")
            if diagnostic.get("status") != "pass":
                reasons.append(f"{realization_id}:structural:{metric}:insufficient_sampling")
            recomputed_stationarity = stationarity_diagnostics(
                values[window],
                times[window],
                manifest["diagnostics"]["stationarity"],
                float(manifest["acceptance_gates"]["stationarity_scale_floors"]["structural"][metric]),
            )
            reported_stationarity = diagnostic.get("stationarity", {})
            verify_stationarity_report(reported_stationarity, recomputed_stationarity, f"{realization_id}/structural/{metric}")
            if recomputed_stationarity["status"] != "pass":
                reasons.append(f"{realization_id}:structural:{metric}:nonstationary")
        _, observed_flags = read_csv(directory / "structural_first_difference_review_flags.csv")
        expected_flags: list[dict[str, Any]] = []
        for metric in STRUCTURAL_REQUIRED:
            values = finite_column(rows, metric)
            _, flags = robust_first_difference(
                values,
                times,
                float(manifest["diagnostics"]["robust_first_difference_z_threshold"]),
            )
            expected_flags.extend({"metric": metric, **flag} for flag in flags)
        require(len(observed_flags) == len(expected_flags), f"Structural review-flag count differs: {realization_id}")
        require(len(observed_flags) == int(summary.get("first_difference_review_flag_count", -1)), f"Structural summary review-flag count differs: {realization_id}")
        for observed, expected in zip(observed_flags, expected_flags, strict=True):
            require(observed.get("metric") == expected["metric"], f"Structural review metric differs: {realization_id}")
            require(int(observed["row_index_zero_based"]) == int(expected["row_index_zero_based"]), f"Structural review row differs: {realization_id}")
            require(math.isclose(float(observed["time_after_ns"]), float(expected["time_after_ns"]), rel_tol=0.0, abs_tol=1e-12), f"Structural review time differs: {realization_id}")
            require(observed.get("point_retained") == "True", f"Structural review flag deleted a point: {realization_id}")
            all_flags.append(observed)
        pose_gates = manifest["acceptance_gates"]["native_pose"]
        ligand_rmsd = finite_column(rows, "pocket_aligned_o6u_heavy_rmsd_nm")
        ligand_com = finite_column(rows, "pocket_aligned_o6u_com_displacement_nm")
        native_fraction = finite_column(rows, "native_contact_fraction")
        joint_pose_pass = (
            (ligand_rmsd[window] <= float(pose_gates["maximum_pocket_aligned_o6u_heavy_rmsd_nm"]))
            & (ligand_com[window] <= float(pose_gates["maximum_o6u_com_displacement_nm"]))
            & (native_fraction[window] >= float(pose_gates["minimum_native_contact_fraction"]))
        )
        joint_pose_fraction = float(np.mean(joint_pose_pass))
        reported_pose = summary.get("native_pose_acceptance", {})
        require(math.isclose(float(reported_pose.get("fraction_of_primary_frames_meeting_all_pose_gates", -1.0)), joint_pose_fraction, rel_tol=0.0, abs_tol=1e-12), f"Native-pose gate fraction differs: {realization_id}")
        event_definitions = (
            ("o6u_heavy_rmsd_egress", ligand_rmsd > float(pose_gates["maximum_pocket_aligned_o6u_heavy_rmsd_nm"])),
            ("o6u_com_displacement_egress", ligand_com > float(pose_gates["maximum_o6u_com_displacement_nm"])),
            ("native_contact_loss", native_fraction < float(pose_gates["minimum_native_contact_fraction"])),
        )
        expected_events = []
        for event_type, condition in event_definitions:
            expected_events.extend(continuous_true_events(condition, times, float(pose_gates["minimum_continuous_event_duration_ns"]), event_type))
        expected_events.sort(key=lambda item: (float(item["start_time_ns"]), str(item["event_type"])))
        _, event_rows = read_csv(directory / "native_pose_continuous_events.csv")
        require(len(event_rows) == len(expected_events), f"Continuous ligand-event count differs: {realization_id}")
        for observed, expected in zip(event_rows, expected_events, strict=True):
            require(observed["event_type"] == expected["event_type"], f"Continuous ligand-event type differs: {realization_id}")
            for key in ("start_frame_index_zero_based", "end_frame_index_zero_based", "frame_count"):
                require(int(observed[key]) == int(expected[key]), f"Continuous ligand-event {key} differs: {realization_id}")
            for key in ("start_time_ns", "end_time_ns", "continuous_duration_ns", "minimum_duration_ns"):
                require(math.isclose(float(observed[key]), float(expected[key]), rel_tol=0.0, abs_tol=1e-12), f"Continuous ligand-event {key} differs: {realization_id}")
        require(int(reported_pose.get("qualifying_continuous_event_count", -1)) == len(expected_events), f"Reported ligand-event count differs: {realization_id}")
        pose_pass = joint_pose_fraction >= float(pose_gates["minimum_fraction_of_primary_frames_meeting_all_pose_gates"]) and not expected_events
        require(summary.get("scientific_status") == ("pass" if pose_pass else "fail"), f"Structural scientific status differs: {realization_id}")
        if not pose_pass:
            reasons.append(f"{realization_id}:structural:native_pose_scientific_failure_no_rerun")
        contact_columns = [name for name in header if name.startswith("contact__") and name.endswith("__present")]
        require(contact_columns, f"No explicit native-contact columns: {realization_id}")
        for name in contact_columns:
            values = finite_column(rows, name)
            require(np.all(np.isin(values, [0.0, 1.0])), f"Contact presence is not binary: {realization_id}/{name}")
        hbond_columns = [name for name in header if name.startswith("hbond__") and name.endswith("__present")]
        require(hbond_columns, f"No endpoint-resolved hydrogen-bond columns: {realization_id}")
        for name in hbond_columns:
            require(np.all(np.isin(finite_column(rows, name), [0.0, 1.0])), f"Hydrogen-bond presence is not binary: {realization_id}/{name}")
        _, occupancy = read_csv(directory / "native_contact_occupancy.csv")
        require(occupancy, f"Native-contact occupancy is empty: {realization_id}")
        require(np.all((0.0 <= finite_column(occupancy, "occupancy")) & (finite_column(occupancy, "occupancy") <= 1.0)), f"Contact occupancy is outside [0,1]: {realization_id}")
        _, rmsf = read_csv(directory / "protein_ca_rmsf_200_500ns.csv")
        require(rmsf and np.all(finite_column(rmsf, "rmsf_nm") >= 0.0), f"Protein RMSF is invalid: {realization_id}")
        time_arrays[realization_id] = times
        if summary.get("technical_status") != "pass" or summary.get("sampling_status") != "pass" or summary.get("scientific_status") != "pass":
            reasons.append(f"{realization_id}:structural_summary_not_pass")
    reference_times = time_arrays["rep01"]
    for realization_id in REALIZATION_IDS[1:]:
        require(np.array_equal(time_arrays[realization_id], reference_times), "Structural saved times differ among realizations")
    return complete, time_arrays, sorted(set(reasons)), all_flags


def validate_membrane(output_root: Path, manifest_path: Path, manifest: Mapping[str, Any], structural_times: Mapping[str, np.ndarray]) -> tuple[dict[str, Any], list[str], list[dict[str, str]]]:
    root = output_root / "membrane_qc"
    complete = verify_complete(root / "COMPLETE.json", manifest_path, {"pass", "inconclusive"})
    require(complete.get("cell_lateral_area_is_area_per_lipid") is False, "Cell lateral area was mislabeled as area per lipid")
    require(complete.get("preproduction_status") in {"pass", "blocked_external_membrane_metrics"}, "Membrane pre-production status differs")
    reasons: list[str] = []
    all_flags: list[dict[str, str]] = []
    observed_no_go: set[str] = set()
    for realization_id in REALIZATION_IDS:
        directory = root / realization_id
        summary = load_json(directory / "membrane_summary.json")
        require(summary.get("realization_id") == realization_id, f"Membrane summary ID differs: {realization_id}")
        require(summary.get("cell_lateral_area_is_area_per_lipid") is False, f"Cell area mislabeled as APL: {realization_id}")
        require(summary.get("first_difference_flags_remove_points") is False, f"Membrane flags remove points: {realization_id}")
        require(summary.get("data_handling") == {"raw_frames_retained": True, "smoothing": False, "frame_deletion": False, "interpolation": False, "realization_pooling": False}, f"Membrane data-handling contract differs: {realization_id}")
        verify_summary_hashes(summary, directory)
        header, rows = read_csv(directory / "membrane_raw_unsmoothed.csv")
        require(rows and all(name in header for name in MEMBRANE_REQUIRED), f"Membrane raw CSV lacks required metrics: {realization_id}")
        require(len(rows) == int(summary["input_frame_count"]), f"Membrane row count differs: {realization_id}")
        times = finite_column(rows, "time_ns")
        require(np.array_equal(times, structural_times[realization_id]), f"Structural and membrane frame accounting differs: {realization_id}")
        for metric in ("phosphate_peak_thickness_nm", "cell_lateral_area_nm2_not_apl", "box_z_vector_length_nm", "cell_volume_nm3"):
            require(np.all(finite_column(rows, metric) > 0.0), f"Nonpositive membrane/cell metric: {realization_id}/{metric}")
        tilt = finite_column(rows, "protein_tilt_deg")
        require(np.all((0.0 <= tilt) & (tilt <= 90.0)), f"Protein tilt is outside [0,90]: {realization_id}")
        for metric in ("outside_protein_hydrophobic_core_water_count", "outside_protein_hydrophobic_core_largest_water_cluster"):
            values = finite_column(rows, metric)
            require(np.all(values >= 0.0) and np.all(values == np.floor(values)), f"Water-defect counts are invalid: {realization_id}/{metric}")
        availability = load_json(directory / "membrane_metric_availability.json")
        require(availability.get("cell_lateral_area", {}).get("label") == "periodic cell lateral area; not area per lipid", "Cell-area label differs")
        for metric in ("protein_aware_area_per_lipid", "popc_deuterium_order_parameters"):
            metric_reasons = verify_external_membrane_metric(
                availability.get(metric, {}),
                metric,
                realization_id,
                times,
                str(summary.get("input_files", {}).get("centered_system_trajectory", {}).get("sha256", "")),
                str(summary.get("trajectory_atom_identity_sha256", "")),
                manifest,
            )
            reasons.extend(metric_reasons)
            if metric_reasons:
                observed_no_go.add(metric)
        for metric in ("phosphate_peak_thickness_nm", "cell_lateral_area_nm2_not_apl", "box_z_vector_length_nm", "protein_tilt_deg"):
            diagnostic = summary.get("diagnostics", {}).get(metric, {}).get("primary_200_500_ns", {})
            require(diagnostic.get("minimum_complete_blocks", 0) >= 5 and float(diagnostic.get("block_tau_multiplier", 0.0)) >= 10.0, f"Membrane block gate was weakened: {realization_id}/{metric}")
            if diagnostic.get("status") != "pass":
                reasons.append(f"{realization_id}:membrane:{metric}:insufficient_sampling")
            values = finite_column(rows, metric)
            primary = primary_window_mask(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
            recomputed_stationarity = stationarity_diagnostics(
                values[primary],
                times[primary],
                manifest["diagnostics"]["stationarity"],
                float(manifest["acceptance_gates"]["stationarity_scale_floors"]["membrane"][metric]),
            )
            verify_stationarity_report(diagnostic.get("stationarity", {}), recomputed_stationarity, f"{realization_id}/membrane/{metric}")
            if recomputed_stationarity["status"] != "pass":
                reasons.append(f"{realization_id}:membrane:{metric}:nonstationary")
        _, flags = read_csv(directory / "membrane_first_difference_review_flags.csv")
        require(len(flags) == int(summary["first_difference_review_flag_count"]), f"Membrane review-flag count differs: {realization_id}")
        for row in flags:
            require(row.get("point_retained") == "True", f"Membrane review flag deleted a point: {realization_id}")
            all_flags.append(row)
        if summary.get("technical_status") != "pass" or summary.get("sampling_status") != "pass" or summary.get("qc_failures"):
            reasons.append(f"{realization_id}:membrane_qc_not_pass")
    require(complete.get("preproduction_no_go_requirements") == sorted(observed_no_go), "Membrane COMPLETE external-metric NO-GO list differs")
    require(complete.get("preproduction_status") == ("blocked_external_membrane_metrics" if observed_no_go else "pass"), "Membrane COMPLETE external-metric status differs")
    return complete, sorted(set(reasons)), all_flags


def _load_review_dispositions(
    path: Path,
    structural_complete: Path,
    energy_complete: Path,
    membrane_complete: Path,
    structural_flags: list[dict[str, str]],
    energy_flags: list[dict[str, str]],
    membrane_flags: list[dict[str, str]],
) -> list[str]:
    record = load_json(path)
    require(record.get("schema_version") == "3.0", "Review dispositions schema_version must be 3.0")
    require(record.get("approval_status") == "approved", "Review dispositions are not approved")
    require(not has_placeholder(record), "Review dispositions contain TODO/REPLACE_ME placeholders")
    require(record.get("policy_frozen_before_production") is True, "Spike adjudication policy was not frozen before production")
    resolve_record(path.parent, record.get("policy_source_record", {}), "review_dispositions.policy_source_record")
    source_hashes = record.get("source_complete_sha256", {})
    require(source_hashes == {
        "structural": sha256_file(structural_complete),
        "energy": sha256_file(energy_complete),
        "membrane": sha256_file(membrane_complete),
    }, "Review dispositions are bound to different component outputs")
    require(record.get("rules") == {
        "all_flags_covered_exactly_once": True,
        "source_rows_retained": True,
        "unresolved_or_blocking_disposition_blocks_claim": True,
        "every_adjudication_requires_source_hashed_evidence": True,
    }, "Spike adjudication rules differ")
    allowed = record.get("allowed_dispositions", {})
    require(allowed.get("structural") == [
        "finite_coordinate_fluctuation_retained",
        "pbc_or_coordinate_artifact_raw_point_retained_and_analysis_blocked",
        "structural_disruption_or_corruption_analysis_blocked",
    ], "Structural disposition vocabulary differs")
    require(allowed.get("energy") == [
        "finite_physical_or_barostat_fluctuation_retained",
        "pbc_or_output_artifact_raw_point_retained_and_analysis_blocked",
        "constraint_or_corruption_failure_analysis_blocked",
    ], "Energy disposition vocabulary differs")
    require(allowed.get("membrane") == [
        "finite_physical_fluctuation_retained",
        "pbc_or_coordinate_artifact_raw_point_retained_and_analysis_blocked",
        "structural_disruption_or_corruption_analysis_blocked",
    ], "Membrane disposition vocabulary differs")
    reasons: list[str] = []

    def audit(kind: str, flags: list[dict[str, str]], entries: Any, key_field: str) -> None:
        require(isinstance(entries, list), f"{kind} dispositions must be a list")
        observed_keys = [(row["realization_id"], row[key_field], int(row["row_index_zero_based"]), float(row["time_after_ns"])) for row in flags]
        disposition_keys = [(row["realization_id"], row[key_field], int(row["row_index_zero_based"]), float(row["time_after_ns"])) for row in entries]
        require(len(observed_keys) == len(set(observed_keys)), f"{kind} review flags contain duplicate keys")
        require(len(disposition_keys) == len(set(disposition_keys)), f"{kind} dispositions contain duplicate keys")
        require(set(observed_keys) == set(disposition_keys), f"{kind} dispositions do not cover every flag exactly once")
        for entry in entries:
            disposition = str(entry.get("disposition", ""))
            require(disposition in set(allowed[kind]), f"{kind} disposition is not allowed")
            evidence = entry.get("evidence_records")
            require(isinstance(evidence, list) and evidence, f"{kind} disposition lacks source-hashed evidence")
            for number, evidence_record in enumerate(evidence):
                resolve_record(path.parent, evidence_record, f"{kind}_flags.evidence_records[{number}]")
            require(isinstance(entry.get("reviewer"), str) and entry["reviewer"].strip(), f"{kind} disposition lacks a reviewer")
            reviewed_at = str(entry.get("reviewed_at_utc", ""))
            require(reviewed_at.endswith("Z") and "T" in reviewed_at, f"{kind} disposition lacks an ISO-8601 UTC timestamp")
            if disposition.endswith("analysis_blocked"):
                reasons.append(f"{entry['realization_id']}:{kind}:{entry[key_field]}:analysis_blocked")

    audit("structural", structural_flags, record.get("structural_flags"), "metric")
    audit("energy", energy_flags, record.get("energy_flags"), "term_key")
    audit("membrane", membrane_flags, record.get("membrane_flags"), "metric")
    return reasons


def validate_energy(output_root: Path, manifest_path: Path, manifest: Mapping[str, Any]) -> tuple[dict[str, Any], list[str], list[dict[str, str]]]:
    root = output_root / "energy_qc"
    complete_path = root / "COMPLETE.json"
    complete = verify_complete(complete_path, manifest_path, {"pass", "pass_with_review_flags", "inconclusive"})
    require(complete.get("exact_required_term_keys") == list(REQUIRED_TERM_KEYS), "Energy required-term list differs")
    require(complete.get("flagged_points_retained") is True, "Energy flags were allowed to remove points")
    if manifest.get("approval_status") == "approved_for_server_execution":
        version = complete.get("gromacs_version", {})
        expected_version = str(manifest["energy_execution"]["required_gromacs_version"])
        require(version.get("expected") == expected_version and version.get("observed") == expected_version, "Executed GROMACS version differs from the frozen manifest")
        version_log = Path(str(version.get("log", "")))
        require(version_log.is_file() and version.get("log_sha256") == sha256_file(version_log), "GROMACS version log is missing or changed")
    reasons: list[str] = []
    all_flags: list[dict[str, str]] = []
    time_arrays = {}
    for realization_id in REALIZATION_IDS:
        directory = root / realization_id
        summary = load_json(directory / "energy_summary.json")
        require(summary.get("realization_id") == realization_id, f"Energy summary ID differs: {realization_id}")
        require(summary.get("technical_status") == "pass", f"Energy technical status failed: {realization_id}")
        require(summary.get("flagged_points_retained") is True and summary.get("smoothing") is False and summary.get("frame_deletion") is False and summary.get("interpolation") is False, f"Energy data handling differs: {realization_id}")
        verify_summary_hashes(summary, directory)
        term_metadata = summary.get("terms", {})
        require(set(term_metadata) == set(REQUIRED_TERM_KEYS), f"Energy artifact metadata differs: {realization_id}")
        raw_gmx_directory = directory / "raw_gmx_energy"
        for term in REQUIRED_TERM_KEYS:
            metadata = term_metadata[term]
            for prefix, hash_key in (("raw", "xvg_sha256"), ("gmx_energy", "log_sha256"), ("selection", "selection_sha256"), ("command", "command_sha256")):
                artifact = raw_gmx_directory / f"{prefix}__{term}.{'xvg' if prefix == 'raw' else 'log' if prefix == 'gmx_energy' else 'txt' if prefix == 'selection' else 'json'}"
                require(artifact.is_file() and metadata.get(hash_key) == sha256_file(artifact), f"Energy extraction artifact hash differs: {artifact}")
        header, rows = read_csv(directory / "energy_raw_unsmoothed.csv")
        require(rows and all(term in header for term in REQUIRED_TERM_KEYS), f"Energy raw CSV lacks exact terms: {realization_id}")
        require(len(rows) == int(summary["input_row_count"]), f"Energy row count differs: {realization_id}")
        times = finite_column(rows, "time_ns")
        validate_time_axis(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
        values = {term: finite_column(rows, term) for term in REQUIRED_TERM_KEYS}
        primary = primary_window_mask(times, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
        for term in REQUIRED_TERM_KEYS:
            reported = summary.get("primary_diagnostics", {}).get(term, {})
            require(reported.get("minimum_complete_blocks", 0) >= 5 and float(reported.get("block_tau_multiplier", 0.0)) >= 10.0, f"Energy block gate was weakened: {realization_id}/{term}")
            recomputed_stationarity = stationarity_diagnostics(
                values[term][primary],
                times[primary],
                manifest["diagnostics"]["stationarity"],
                float(manifest["acceptance_gates"]["stationarity_scale_floors"]["energy"][term]),
            )
            verify_stationarity_report(reported.get("stationarity", {}), recomputed_stationarity, f"{realization_id}/energy/{term}")
            if reported.get("status") != "pass":
                reasons.append(f"{realization_id}:energy:{term}:insufficient_sampling")
            if recomputed_stationarity["status"] != "pass":
                reasons.append(f"{realization_id}:energy:{term}:nonstationary")
        gates = manifest["acceptance_gates"]["thermodynamic_cell_qc"]
        temperature_mean = float(np.mean(values["temperature_k"][primary]))
        pressure_mean = float(np.mean(values["pressure_bar"][primary]))
        density_mean = float(np.mean(values["density_kg_m3"][primary]))
        energy_denominator = np.maximum.reduce([np.ones(len(times)), np.abs(values["total_energy_kj_mol"]), np.abs(values["potential_energy_kj_mol"]) + np.abs(values["kinetic_energy_kj_mol"])])
        energy_closure = float(np.max(np.abs(values["total_energy_kj_mol"] - values["potential_energy_kj_mol"] - values["kinetic_energy_kj_mol"]) / energy_denominator))
        pressure_closure = float(np.max(np.abs(values["pressure_bar"] - (values["pressure_xx_bar"] + values["pressure_yy_bar"] + values["pressure_zz_bar"]) / 3.0)))
        box_volume = values["box_x_nm"] * values["box_y_nm"] * values["box_z_nm"]
        volume_closure = float(np.max(np.abs(values["volume_nm3"] - box_volume) / np.maximum(1.0, np.abs(values["volume_nm3"]))))
        failures = []
        if np.any(values["temperature_k"] <= 0.0) or np.any(values["kinetic_energy_kj_mol"] < 0.0) or np.any(values["density_kg_m3"] <= 0.0) or np.any(values["volume_nm3"] <= 0.0) or any(np.any(values[key] <= 0.0) for key in ("box_x_nm", "box_y_nm", "box_z_nm")):
            failures.append("nonphysical_positive_quantity_gate")
        if abs(temperature_mean - float(gates["target_temperature_k"])) > float(gates["maximum_absolute_primary_mean_temperature_deviation_k"]):
            failures.append("primary_mean_temperature_gate")
        pressure_range = [float(value) for value in gates["approved_primary_mean_pressure_range_bar"]]
        density_range = [float(value) for value in gates["approved_primary_mean_density_range_kg_m3"]]
        if not pressure_range[0] <= pressure_mean <= pressure_range[1]:
            failures.append("primary_mean_pressure_gate")
        if not density_range[0] <= density_mean <= density_range[1]:
            failures.append("primary_mean_density_gate")
        if energy_closure > float(gates["maximum_relative_total_energy_closure_error"]):
            failures.append("total_energy_closure_gate")
        if pressure_closure > float(gates["maximum_absolute_pressure_trace_closure_bar"]):
            failures.append("pressure_trace_closure_gate")
        if volume_closure > float(gates["maximum_relative_orthorhombic_volume_closure_error"]):
            failures.append("orthorhombic_volume_closure_gate")
        reported_acceptance = summary.get("thermodynamic_acceptance", {})
        require(reported_acceptance.get("status") == ("pass" if not failures else "scientific_fail"), f"Energy scientific gate status differs: {realization_id}")
        require(summary.get("scientific_status") == ("pass" if not failures else "fail"), f"Energy scientific status differs: {realization_id}")
        for failure in failures:
            reasons.append(f"{realization_id}:energy:{failure}:scientific_failure_no_rerun")
        _, flags = read_csv(directory / "energy_first_difference_review_flags.csv")
        require(len(flags) == int(summary["review_flag_count"]), f"Energy review-flag count differs: {realization_id}")
        for row in flags:
            require(row.get("point_retained") == "True", f"Energy flag deleted a point: {realization_id}")
            all_flags.append(row)
        time_arrays[realization_id] = times
    for realization_id in REALIZATION_IDS[1:]:
        require(np.array_equal(time_arrays[realization_id], time_arrays["rep01"]), "Energy saved times differ among realizations")
    return complete, sorted(set(reasons)), all_flags


def validate(manifest_path: Path, output_root: Path, dispositions_path: Path | None = None, allow_synthetic: bool = False) -> dict[str, Any]:
    manifest_path = manifest_path.resolve()
    output_root = output_root.resolve()
    manifest, _ = validate_primary_manifest(manifest_path, allow_synthetic=allow_synthetic)
    structural_complete, structural_times, structural_reasons, structural_flags = validate_structural(output_root, manifest_path, manifest)
    membrane_complete, membrane_reasons, membrane_flags = validate_membrane(output_root, manifest_path, manifest, structural_times)
    energy_complete, energy_reasons, energy_flags = validate_energy(output_root, manifest_path, manifest)
    review_reasons: list[str] = []
    if structural_flags or energy_flags or membrane_flags:
        if dispositions_path is None:
            if structural_flags:
                review_reasons.append("structural_first_difference_flags_pending_source_hashed_adjudication")
            if energy_flags:
                review_reasons.append("energy_first_difference_flags_pending_source_hashed_adjudication")
            if membrane_flags:
                review_reasons.append("membrane_first_difference_flags_pending_source_hashed_adjudication")
        else:
            review_reasons.extend(_load_review_dispositions(
                dispositions_path.resolve(),
                output_root / "structural_analysis" / "COMPLETE.json",
                output_root / "energy_qc" / "COMPLETE.json",
                output_root / "membrane_qc" / "COMPLETE.json",
                structural_flags,
                energy_flags,
                membrane_flags,
            ))
    elif dispositions_path is not None:
        review_reasons.extend(_load_review_dispositions(
            dispositions_path.resolve(),
            output_root / "structural_analysis" / "COMPLETE.json",
            output_root / "energy_qc" / "COMPLETE.json",
            output_root / "membrane_qc" / "COMPLETE.json",
            structural_flags,
            energy_flags,
            membrane_flags,
        ))
    reasons = sorted(set(structural_reasons + membrane_reasons + energy_reasons + review_reasons))
    status = "pass" if not reasons else "inconclusive"
    report = {
        "schema_version": "1.0",
        "status": status,
        "claim_gate": "eligible_for_bounded_native_pose_compatibility_statement" if status == "pass" else "blocked_inconclusive",
        "system_id": manifest["system_id"],
        "construction_count": 1,
        "realization_ids": list(REALIZATION_IDS),
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "extension_or_recovery_window": False,
        "all_three_realizations_required": True,
        "failed_realizations_may_be_dropped": False,
        "reasons": reasons,
        "manifest_sha256": sha256_file(manifest_path),
        "component_complete_sha256": {
            "structural": sha256_file(output_root / "structural_analysis" / "COMPLETE.json"),
            "membrane": sha256_file(output_root / "membrane_qc" / "COMPLETE.json"),
            "energy": sha256_file(output_root / "energy_qc" / "COMPLETE.json"),
        },
        "component_statuses": {
            "structural": structural_complete["status"],
            "membrane": membrane_complete["status"],
            "energy": energy_complete["status"],
        },
        "prohibited": {"smoothing": False, "frame_deletion": False, "interpolation": False, "mmgbsa_mmpbsa": False, "recovery_extension": False},
    }
    atomic_write_json(output_root / "primary_postprocessing_validation.json", report)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--review-dispositions", "--energy-dispositions", dest="review_dispositions", type=Path, help="Approved source-hashed exact coverage of all structural, membrane, and energy spike flags")
    parser.add_argument("--allow-synthetic", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = validate(args.manifest, args.output_root, args.review_dispositions, allow_synthetic=args.allow_synthetic)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
