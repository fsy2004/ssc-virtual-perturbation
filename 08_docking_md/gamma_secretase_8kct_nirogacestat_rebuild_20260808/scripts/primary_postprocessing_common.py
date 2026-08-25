#!/usr/bin/env python3
"""Shared fail-closed utilities for the primary 8KCT one-system analyses."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


SYSTEM_ID = "8kct_nirogacestat_native"
STUDY_ID = "gamma_secretase_native_nirogacestat_rebuild_20260808"
REALIZATION_IDS = ("rep01", "rep02", "rep03")
EXPECTED_DURATION_NS = 500.0
PRIMARY_WINDOW_NS = (200.0, 500.0)
REQUIRED_MDANALYSIS_VERSION = (2, 10)

STRUCTURAL_STATIONARITY_METRICS = (
    "pocket_aligned_o6u_heavy_rmsd_nm",
    "pocket_aligned_o6u_com_displacement_nm",
    "tm_core_ca_rmsd_nm",
    "protein_ca_rmsd_nm",
    "native_contact_fraction",
)
MEMBRANE_STATIONARITY_METRICS = (
    "phosphate_peak_thickness_nm",
    "protein_aware_area_per_lipid_nm2",
    "cell_lateral_area_nm2_not_apl",
    "box_z_vector_length_nm",
    "cell_volume_nm3",
    "protein_tilt_deg",
)
ENERGY_STATIONARITY_METRICS = (
    "temperature_k",
    "pressure_bar",
    "pressure_xx_bar",
    "pressure_yy_bar",
    "pressure_zz_bar",
    "pressure_xy_bar",
    "pressure_xz_bar",
    "pressure_yz_bar",
    "potential_energy_kj_mol",
    "kinetic_energy_kj_mol",
    "total_energy_kj_mol",
    "density_kg_m3",
    "volume_nm3",
    "box_x_nm",
    "box_y_nm",
    "box_z_nm",
)


class ContractError(RuntimeError):
    """Raised when a frozen input or output violates the analysis contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"Cannot read JSON object {path}: {exc}") from exc
    require(isinstance(value, dict), f"Expected a JSON object: {path}")
    return value


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_jsonable(item) for item in value.tolist()]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        require(math.isfinite(value), "Refusing to serialize a nonfinite JSON value")
    return value


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(_jsonable(payload), handle, indent=2, sort_keys=True, ensure_ascii=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def atomic_write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    count = 0
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="raise")
            writer.writeheader()
            for row in rows:
                converted = {key: _jsonable(row.get(key)) for key in fieldnames}
                writer.writerow(converted)
                count += 1
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    return count


def has_placeholder(value: Any) -> bool:
    if isinstance(value, str):
        upper = value.upper()
        return "TODO" in upper or "REPLACE_ME" in upper
    if isinstance(value, Mapping):
        return any(has_placeholder(key) or has_placeholder(item) for key, item in value.items())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return any(has_placeholder(item) for item in value)
    return False


def resolve_record(base: Path, record: Mapping[str, Any], label: str) -> Path:
    require(isinstance(record, Mapping), f"{label} must be a file record")
    raw_path = record.get("path")
    expected_hash = record.get("sha256")
    require(isinstance(raw_path, str) and raw_path and not has_placeholder(raw_path), f"{label}.path is not frozen")
    require(isinstance(expected_hash, str) and len(expected_hash) == 64 and not has_placeholder(expected_hash), f"{label}.sha256 is not frozen")
    path = Path(raw_path)
    if not path.is_absolute():
        path = (base / path).resolve()
    require(path.is_file() and path.stat().st_size > 0, f"Missing or empty {label}: {path}")
    actual_hash = sha256_file(path)
    require(actual_hash.lower() == expected_hash.lower(), f"SHA-256 mismatch for {label}: {path}")
    return path


def manifest_package_root(manifest_path: Path) -> Path:
    """Return the package root used by manifest-relative file records.

    Production manifests live under ``config/`` and all recorded paths are
    package-root relative. Synthetic self-test manifests may live directly in
    their temporary package root.
    """
    resolved = manifest_path.resolve()
    return resolved.parent.parent if resolved.parent.name == "config" else resolved.parent


def validate_primary_manifest(manifest_path: Path, allow_synthetic: bool = False) -> tuple[dict[str, Any], Path]:
    manifest_path = manifest_path.resolve()
    manifest = load_json(manifest_path)
    base = manifest_package_root(manifest_path)
    require(manifest.get("schema_version") == "1.0", "Primary manifest schema_version must be 1.0")
    allowed_status = {"approved_for_server_execution"}
    if allow_synthetic:
        allowed_status.add("synthetic_self_test")
    require(manifest.get("approval_status") in allowed_status, "Primary manifest is not approved for execution")
    require(not has_placeholder(manifest), "Primary manifest still contains TODO/REPLACE_ME placeholders")
    require(manifest.get("study_id") == STUDY_ID, f"study_id must be {STUDY_ID}")
    require(manifest.get("system_id") == SYSTEM_ID, f"system_id must be {SYSTEM_ID}")
    require(manifest.get("construction_count") == 1, "Exactly one membrane construction is required")
    require(manifest.get("required_realization_ids") == list(REALIZATION_IDS), "Exactly rep01-rep03 are required in fixed order")
    require(float(manifest.get("production_duration_ns")) == EXPECTED_DURATION_NS, "Production duration must be exactly 500 ns")
    require(tuple(float(x) for x in manifest.get("primary_analysis_window_ns", [])) == PRIMARY_WINDOW_NS, "Primary analysis window must be exactly 200-500 ns")
    require(manifest.get("extension_or_recovery_window") is False, "Extension/recovery logic is prohibited")
    for key in ("smooth_frames", "delete_frames", "interpolate_frames", "drop_realizations"):
        require(manifest.get("data_handling", {}).get(key) is False, f"data_handling.{key} must be false")
    require(manifest.get("data_handling", {}).get("frame_stride") == 1, "Every saved frame must be retained (frame_stride=1)")
    require(manifest.get("inference_contract", {}).get("frame_level_hypothesis_tests") is False, "Frame-level inference is prohibited")
    require(manifest.get("inference_contract", {}).get("mmgbsa_mmpbsa") is False, "MM/GBSA and MM/PBSA are prohibited")
    require(manifest.get("required_mdanalysis_version") == "2.10.x", "required_mdanalysis_version must be 2.10.x")
    time_contract = manifest.get("time_contract", {})
    require(time_contract.get("mdanalysis_trajectory_time_unit") == "ps" and time_contract.get("reported_time_unit") == "ns", "Trajectory/report time units must be ps/ns")
    require(float(time_contract.get("required_start_ns")) == 0.0 and float(time_contract.get("required_end_ns")) == EXPECTED_DURATION_NS, "Time endpoints must be frozen at 0 and 500 ns")
    require(float(time_contract.get("endpoint_tolerance_ns")) == 0.0001, "Endpoint tolerance must be exactly 0.0001 ns")
    for key in ("strictly_increasing", "uniform_saved_step", "identical_saved_times_across_realizations"):
        require(time_contract.get(key) is True, f"time_contract.{key} must be true")
    realizations = manifest.get("realizations")
    require(isinstance(realizations, list) and len(realizations) == 3, "Manifest must contain exactly three realization records")
    require([item.get("realization_id") for item in realizations] == list(REALIZATION_IDS), "Realization order/identity must be rep01-rep03")
    require(len({item.get("velocity_seed") for item in realizations}) == 3, "Velocity seeds must be unique")
    diagnostics = manifest.get("diagnostics", {})
    require(diagnostics.get("autocorrelation_method") == "geyer_initial_positive_sequence", "Autocorrelation method is not frozen")
    require(float(diagnostics.get("max_lag_fraction")) == 0.25, "max_lag_fraction must be exactly 0.25")
    require(float(diagnostics.get("block_tau_multiplier")) >= 10.0, "Block length must be at least 10 integrated autocorrelation times")
    require(int(diagnostics.get("minimum_complete_blocks")) >= 5, "At least five complete blocks are required")
    require(float(diagnostics.get("robust_first_difference_z_threshold")) == 12.0, "Robust first-difference threshold must be exactly 12")
    stationarity = diagnostics.get("stationarity", {})
    require(stationarity.get("method") == "fixed_time_block_median_drift_and_change", "Stationarity method is not frozen")
    require(stationarity.get("window_ns") == list(PRIMARY_WINDOW_NS), "Stationarity window must be exactly 200-500 ns")
    require(int(stationarity.get("fixed_time_blocks")) == 5, "Stationarity must use five fixed 60 ns calendar-time blocks")
    require(int(stationarity.get("minimum_frames_per_block")) >= 5, "Each stationarity block must contain at least five saved frames")
    require(stationarity.get("scale_estimator") == "max_1.4826_mad_and_prespecified_metric_floor", "Stationarity scale estimator differs")
    for key in (
        "maximum_abs_normalized_linear_change",
        "maximum_abs_normalized_first_last_shift",
        "maximum_abs_normalized_adjacent_shift",
        "maximum_abs_normalized_change_point_shift",
    ):
        value = float(stationarity.get(key))
        require(math.isfinite(value) and value > 0.0, f"diagnostics.stationarity.{key} must be a frozen positive value")
    acceptance = manifest.get("acceptance_gates", {})
    allowed_acceptance_status = {"approved_and_frozen_before_production"}
    if allow_synthetic:
        allowed_acceptance_status.add("synthetic_self_test")
    require(acceptance.get("approval_status") in allowed_acceptance_status, "Scientific acceptance gates are not approved and frozen before production")
    require(acceptance.get("frozen_before_production_review") is True, "Acceptance gates must be frozen before production-trajectory review")
    require(isinstance(acceptance.get("rationale"), str) and len(acceptance["rationale"]) >= 30, "Acceptance gates require a concrete rationale")
    resolve_record(base, acceptance.get("source_record", {}), "acceptance_gates.source_record")
    scale_floors = acceptance.get("stationarity_scale_floors", {})
    expected_scale_floors = {
        "structural": set(STRUCTURAL_STATIONARITY_METRICS),
        "membrane": set(MEMBRANE_STATIONARITY_METRICS),
        "energy": set(ENERGY_STATIONARITY_METRICS),
    }
    require(set(scale_floors) == set(expected_scale_floors), "Stationarity scale-floor analysis groups differ")
    for group, expected_metrics in expected_scale_floors.items():
        observed = scale_floors.get(group)
        require(isinstance(observed, Mapping) and set(observed) == expected_metrics, f"Stationarity scale floors differ for {group}")
        for metric, raw_value in observed.items():
            value = float(raw_value)
            require(math.isfinite(value) and value > 0.0, f"Stationarity scale floor must be positive: {group}.{metric}")
    native_pose = acceptance.get("native_pose", {})
    require(native_pose.get("units") == {"distance": "nm", "fraction": "unitless"}, "Native-pose gate units differ")
    require(native_pose.get("window_ns") == list(PRIMARY_WINDOW_NS), "Native-pose gates must use 200-500 ns")
    require(native_pose.get("event_search_window_ns") == [0.0, EXPECTED_DURATION_NS], "Ligand egress/contact-loss events must be searched over the full 0-500 ns production trace")
    require(native_pose.get("continuous_event_rule") == "each_geometry_evaluated_separately_without_gap_bridging", "Continuous ligand-event rule differs")
    minimum_event_duration = float(native_pose.get("minimum_continuous_event_duration_ns"))
    require(math.isfinite(minimum_event_duration) and 0.0 < minimum_event_duration <= 5.0, "Continuous ligand-event duration must be frozen in (0,5] ns")
    require(native_pose.get("ligand_egress_or_contact_loss_is_scientific_failure") is True, "Ligand egress/contact loss must be a scientific failure")
    require(native_pose.get("failure_triggers_rerun_or_extension") is False, "A native-pose failure must not trigger rerun/extension")
    maximum_rmsd = float(native_pose.get("maximum_pocket_aligned_o6u_heavy_rmsd_nm"))
    maximum_com = float(native_pose.get("maximum_o6u_com_displacement_nm"))
    minimum_contact = float(native_pose.get("minimum_native_contact_fraction"))
    minimum_frame_fraction = float(native_pose.get("minimum_fraction_of_primary_frames_meeting_all_pose_gates"))
    require(math.isfinite(maximum_rmsd) and maximum_rmsd > 0.0, "Native-pose RMSD cutoff is invalid")
    require(math.isfinite(maximum_com) and maximum_com > 0.0, "Native-pose COM cutoff is invalid")
    require(0.0 <= minimum_contact <= 1.0 and 0.0 < minimum_frame_fraction <= 1.0, "Native-pose fraction gates are invalid")
    thermodynamic = acceptance.get("thermodynamic_cell_qc", {})
    require(thermodynamic.get("units") == {"temperature": "K", "pressure": "bar", "density": "kg/m^3", "relative_closure": "unitless"}, "Thermodynamic gate units differ")
    require(thermodynamic.get("window_ns") == list(PRIMARY_WINDOW_NS), "Thermodynamic gates must use 200-500 ns")
    target_temperature = float(thermodynamic.get("target_temperature_k"))
    temperature_tolerance = float(thermodynamic.get("maximum_absolute_primary_mean_temperature_deviation_k"))
    pressure_range = [float(value) for value in thermodynamic.get("approved_primary_mean_pressure_range_bar", [])]
    density_range = [float(value) for value in thermodynamic.get("approved_primary_mean_density_range_kg_m3", [])]
    require(math.isfinite(target_temperature) and target_temperature > 0.0 and math.isfinite(temperature_tolerance) and temperature_tolerance > 0.0, "Temperature gates are invalid")
    require(len(pressure_range) == 2 and all(math.isfinite(value) for value in pressure_range) and pressure_range[0] < pressure_range[1], "Pressure mean range is invalid")
    require(len(density_range) == 2 and all(math.isfinite(value) and value > 0.0 for value in density_range) and density_range[0] < density_range[1], "Density mean range is invalid")
    for key in ("maximum_relative_total_energy_closure_error", "maximum_absolute_pressure_trace_closure_bar", "maximum_relative_orthorhombic_volume_closure_error"):
        value = float(thermodynamic.get(key))
        require(math.isfinite(value) and value > 0.0, f"thermodynamic_cell_qc.{key} is invalid")
    require(thermodynamic.get("failure_triggers_rerun_or_extension") is False, "Thermodynamic scientific failure must not trigger automatic rerun/extension")
    return manifest, base


def check_mdanalysis_version(version: str) -> None:
    pieces = version.split(".")
    try:
        observed = (int(pieces[0]), int(pieces[1]))
    except (ValueError, IndexError) as exc:
        raise ContractError(f"Cannot parse MDAnalysis version: {version}") from exc
    require(observed == REQUIRED_MDANALYSIS_VERSION, f"MDAnalysis 2.10.x is required; observed {version}")


def validate_time_axis(times_ns: np.ndarray, endpoint_tolerance_ns: float) -> float:
    times = np.asarray(times_ns, dtype=np.float64)
    require(times.ndim == 1 and len(times) >= 2, "At least two trajectory frames are required")
    require(np.all(np.isfinite(times)), "Trajectory times contain NaN or infinity")
    differences = np.diff(times)
    require(np.all(differences > 0.0), "Trajectory times are not strictly increasing")
    median_step = float(np.median(differences))
    require(np.allclose(differences, median_step, rtol=1e-5, atol=max(1e-9, endpoint_tolerance_ns / 10.0)), "Trajectory has missing, duplicate, or nonuniform saved frames")
    require(abs(float(times[0])) <= endpoint_tolerance_ns, f"Trajectory must start at 0 ns; observed {times[0]:.9g}")
    require(abs(float(times[-1]) - EXPECTED_DURATION_NS) <= endpoint_tolerance_ns, f"Trajectory must end at 500 ns; observed {times[-1]:.9g}")
    return median_step


def primary_window_mask(times_ns: np.ndarray, tolerance_ns: float) -> np.ndarray:
    times = np.asarray(times_ns, dtype=np.float64)
    mask = (times >= PRIMARY_WINDOW_NS[0] - tolerance_ns) & (times <= PRIMARY_WINDOW_NS[1] + tolerance_ns)
    require(np.count_nonzero(mask) >= 2, "No complete 200-500 ns primary window is present")
    selected = times[mask]
    require(abs(float(selected[0]) - PRIMARY_WINDOW_NS[0]) <= tolerance_ns, "Primary window lacks the 200 ns endpoint")
    require(abs(float(selected[-1]) - PRIMARY_WINDOW_NS[1]) <= tolerance_ns, "Primary window lacks the 500 ns endpoint")
    return mask


def kabsch_transform(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile = np.asarray(mobile, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    require(mobile.shape == reference.shape and mobile.ndim == 2 and mobile.shape[1] == 3, "Kabsch arrays must be matching N x 3 matrices")
    require(mobile.shape[0] >= 3, "At least three mapped alignment atoms are required")
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_t = np.linalg.svd(covariance)
    if np.linalg.det(left @ right_t) < 0.0:
        left[:, -1] *= -1.0
    rotation = left @ right_t
    require(np.linalg.det(rotation) > 0.999999, "Kabsch mapping produced an improper rotation")
    return rotation, mobile_center, reference_center


def apply_transform(coordinates: np.ndarray, rotation: np.ndarray, mobile_center: np.ndarray, reference_center: np.ndarray) -> np.ndarray:
    return (np.asarray(coordinates, dtype=np.float64) - mobile_center) @ rotation + reference_center


def rmsd_nm(mobile_angstrom: np.ndarray, reference_angstrom: np.ndarray) -> float:
    difference = np.asarray(mobile_angstrom, dtype=np.float64) - np.asarray(reference_angstrom, dtype=np.float64)
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))) / 10.0)


def atom_identity(atom: Any) -> dict[str, Any]:
    identity: dict[str, Any] = {"index": int(atom.index), "name": str(atom.name), "resname": str(atom.resname), "resid": int(atom.resid)}
    for field in ("segid", "chainID"):
        try:
            value = getattr(atom, field)
        except (AttributeError, NoDataError):  # type: ignore[name-defined]
            continue
        if value not in (None, ""):
            identity[field] = str(value)
    return identity


def verify_atom(atom: Any, expected: Mapping[str, Any], label: str) -> None:
    require(isinstance(expected, Mapping), f"{label} expected identity must be an object")
    required = ("index", "name", "resname", "resid")
    for field in required:
        require(field in expected, f"{label} identity lacks {field}")
    observed = atom_identity(atom)
    for field, expected_value in expected.items():
        if expected_value is None:
            continue
        require(field in observed, f"{label} topology lacks identity field {field}")
        if field in ("index", "resid"):
            require(int(observed[field]) == int(expected_value), f"{label} {field} mismatch: {observed[field]} != {expected_value}")
        else:
            require(str(observed[field]) == str(expected_value), f"{label} {field} mismatch: {observed[field]} != {expected_value}")


def mapped_indices(reference_universe: Any, trajectory_universe: Any, entries: Sequence[Mapping[str, Any]], label: str) -> tuple[np.ndarray, np.ndarray]:
    require(isinstance(entries, Sequence) and len(entries) >= 1, f"{label} mapping is empty")
    reference_indices: list[int] = []
    trajectory_indices: list[int] = []
    for number, entry in enumerate(entries):
        require(isinstance(entry, Mapping), f"{label}[{number}] is not an object")
        reference_expected = entry.get("reference")
        trajectory_expected = entry.get("trajectory")
        require(isinstance(reference_expected, Mapping) and isinstance(trajectory_expected, Mapping), f"{label}[{number}] lacks explicit reference/trajectory identities")
        reference_index = int(reference_expected.get("index", -1))
        trajectory_index = int(trajectory_expected.get("index", -1))
        require(0 <= reference_index < len(reference_universe.atoms), f"{label}[{number}] reference index is out of range")
        require(0 <= trajectory_index < len(trajectory_universe.atoms), f"{label}[{number}] trajectory index is out of range")
        verify_atom(reference_universe.atoms[reference_index], reference_expected, f"{label}[{number}].reference")
        verify_atom(trajectory_universe.atoms[trajectory_index], trajectory_expected, f"{label}[{number}].trajectory")
        reference_indices.append(reference_index)
        trajectory_indices.append(trajectory_index)
    require(len(set(reference_indices)) == len(reference_indices), f"{label} contains duplicate reference indices")
    require(len(set(trajectory_indices)) == len(trajectory_indices), f"{label} contains duplicate trajectory indices")
    return np.asarray(reference_indices, dtype=np.int64), np.asarray(trajectory_indices, dtype=np.int64)


def integrated_autocorrelation_time(values: Sequence[float], max_lag_fraction: float) -> dict[str, Any]:
    series = np.asarray(values, dtype=np.float64)
    require(series.ndim == 1 and len(series) >= 2, "Autocorrelation input must be a finite one-dimensional series")
    require(np.all(np.isfinite(series)), "Autocorrelation input contains NaN or infinity")
    centered = series - np.mean(series)
    variance = float(np.dot(centered, centered) / len(centered))
    max_lag = max(1, min(len(series) - 1, int(math.floor(len(series) * float(max_lag_fraction)))))
    if variance <= np.finfo(np.float64).eps:
        return {"method": "geyer_initial_positive_sequence", "tau_frames": 0.5, "cutoff_lag_frames": 0, "max_lag_frames": max_lag, "constant_series": True}
    fft_length = 1 << (2 * len(series) - 1).bit_length()
    transformed = np.fft.rfft(centered, n=fft_length)
    acovariance = np.fft.irfft(transformed * np.conjugate(transformed), n=fft_length)[: len(series)]
    normalization = np.arange(len(series), 0, -1, dtype=np.float64)
    acovariance = acovariance / normalization
    autocorrelation = acovariance / acovariance[0]
    # Geyer's initial positive sequence uses consecutive paired normalized
    # autocovariances: (rho_0 + rho_1), (rho_2 + rho_3), ... .  Stop before
    # the first non-positive pair.  In half-integrated-time convention this is
    # tau = -1/2 + sum(positive pairs).
    tau = -0.5
    cutoff = 0
    for first_lag in range(0, max_lag, 2):
        second_lag = first_lag + 1
        if second_lag > max_lag:
            break
        pair_sum = float(autocorrelation[first_lag] + autocorrelation[second_lag])
        if not math.isfinite(pair_sum) or pair_sum <= 0.0:
            break
        tau += pair_sum
        cutoff = second_lag
    return {"method": "geyer_initial_positive_sequence", "tau_frames": float(max(0.5, tau)), "cutoff_lag_frames": cutoff, "max_lag_frames": max_lag, "constant_series": False}


def block_diagnostics(values: Sequence[float], times_ns: Sequence[float], diagnostics: Mapping[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    series = np.asarray(values, dtype=np.float64)
    times = np.asarray(times_ns, dtype=np.float64)
    require(len(series) == len(times) and len(series) >= 2, "Block inputs have inconsistent lengths")
    require(np.all(np.isfinite(series)), "Block input contains NaN or infinity")
    step_ns = float(np.median(np.diff(times)))
    autocorrelation = integrated_autocorrelation_time(series, float(diagnostics["max_lag_fraction"]))
    multiplier = float(diagnostics["block_tau_multiplier"])
    block_frames = max(1, int(math.ceil(multiplier * float(autocorrelation["tau_frames"]))))
    complete_blocks = len(series) // block_frames
    minimum_blocks = int(diagnostics["minimum_complete_blocks"])
    rows: list[dict[str, Any]] = []
    for index in range(complete_blocks):
        start = index * block_frames
        stop = start + block_frames
        current = series[start:stop]
        rows.append({
            "block_index_zero_based": index,
            "start_time_ns": float(times[start]),
            "end_time_ns": float(times[stop - 1]),
            "frame_count": block_frames,
            "mean": float(np.mean(current)),
            "median": float(np.median(current)),
            "minimum": float(np.min(current)),
            "maximum": float(np.max(current)),
        })
    summary = {
        **autocorrelation,
        "tau_ns": float(autocorrelation["tau_frames"] * step_ns),
        "saved_step_ns": step_ns,
        "block_tau_multiplier": multiplier,
        "block_frames": block_frames,
        "block_duration_ns": float(block_frames * step_ns),
        "complete_blocks": complete_blocks,
        "unused_tail_frames": int(len(series) - complete_blocks * block_frames),
        "minimum_complete_blocks": minimum_blocks,
        "status": "pass" if complete_blocks >= minimum_blocks else "insufficient_sampling",
    }
    return summary, rows


def stationarity_diagnostics(
    values: Sequence[float],
    times_ns: Sequence[float],
    stationarity_contract: Mapping[str, Any],
    metric_scale_floor: float,
) -> dict[str, Any]:
    """Evaluate prespecified robust trend and change gates on fixed time blocks.

    This is an effect-size gate, not a p-value.  Five calendar-time blocks are
    always used over 200-500 ns, so an autocorrelation estimate cannot silently
    redefine the stationarity question after the trace is seen.
    """

    series = np.asarray(values, dtype=np.float64)
    times = np.asarray(times_ns, dtype=np.float64)
    require(series.ndim == times.ndim == 1 and len(series) == len(times), "Stationarity inputs have inconsistent shapes")
    require(len(series) >= 25 and np.all(np.isfinite(series)) and np.all(np.isfinite(times)), "Stationarity inputs are incomplete or nonfinite")
    require(np.all(np.diff(times) > 0.0), "Stationarity times are not strictly increasing")
    require(abs(float(times[0]) - PRIMARY_WINDOW_NS[0]) <= 0.0001 and abs(float(times[-1]) - PRIMARY_WINDOW_NS[1]) <= 0.0001, "Stationarity input must be the exact 200-500 ns window")
    require(stationarity_contract.get("method") == "fixed_time_block_median_drift_and_change", "Stationarity method differs")
    block_count = int(stationarity_contract["fixed_time_blocks"])
    require(block_count == 5, "Exactly five fixed stationarity blocks are required")
    minimum_frames = int(stationarity_contract["minimum_frames_per_block"])
    edges = np.linspace(PRIMARY_WINDOW_NS[0], PRIMARY_WINDOW_NS[1], block_count + 1, dtype=np.float64)
    block_rows: list[dict[str, Any]] = []
    block_medians: list[float] = []
    block_centers: list[float] = []
    for index in range(block_count):
        if index + 1 == block_count:
            mask = (times >= edges[index]) & (times <= edges[index + 1])
        else:
            mask = (times >= edges[index]) & (times < edges[index + 1])
        current = series[mask]
        require(len(current) >= minimum_frames, f"Stationarity block {index} has fewer than {minimum_frames} frames")
        median = float(np.median(current))
        block_medians.append(median)
        center = float((edges[index] + edges[index + 1]) / 2.0)
        block_centers.append(center)
        block_rows.append({
            "block_index_zero_based": index,
            "start_time_ns": float(edges[index]),
            "end_time_ns": float(edges[index + 1]),
            "frame_count": int(len(current)),
            "median": median,
            "mean": float(np.mean(current)),
        })
    raw_median = float(np.median(series))
    raw_mad = float(np.median(np.abs(series - raw_median)))
    robust_scale = max(1.4826 * raw_mad, float(metric_scale_floor))
    require(math.isfinite(robust_scale) and robust_scale > 0.0, "Stationarity robust scale is invalid")
    x = np.asarray(block_centers, dtype=np.float64)
    y = np.asarray(block_medians, dtype=np.float64)
    centered_x = x - np.mean(x)
    slope = float(np.dot(centered_x, y - np.mean(y)) / np.dot(centered_x, centered_x))
    total_linear_change = slope * (PRIMARY_WINDOW_NS[1] - PRIMARY_WINDOW_NS[0])
    normalized_linear = abs(total_linear_change) / robust_scale
    normalized_first_last = abs(float(y[-1] - y[0])) / robust_scale
    normalized_adjacent = float(np.max(np.abs(np.diff(y))) / robust_scale)
    change_candidates = []
    for split in range(2, block_count - 1 + 1):
        if block_count - split < 2:
            continue
        change_candidates.append(abs(float(np.median(y[:split]) - np.median(y[split:]))) / robust_scale)
    require(change_candidates, "Stationarity change-point audit has no eligible split")
    normalized_change = float(max(change_candidates))
    observed = {
        "abs_normalized_linear_change": normalized_linear,
        "abs_normalized_first_last_shift": normalized_first_last,
        "maximum_abs_normalized_adjacent_shift": normalized_adjacent,
        "maximum_abs_normalized_change_point_shift": normalized_change,
    }
    limits = {
        "abs_normalized_linear_change": float(stationarity_contract["maximum_abs_normalized_linear_change"]),
        "abs_normalized_first_last_shift": float(stationarity_contract["maximum_abs_normalized_first_last_shift"]),
        "maximum_abs_normalized_adjacent_shift": float(stationarity_contract["maximum_abs_normalized_adjacent_shift"]),
        "maximum_abs_normalized_change_point_shift": float(stationarity_contract["maximum_abs_normalized_change_point_shift"]),
    }
    failed = [name for name, value in observed.items() if value > limits[name]]
    return {
        "method": "fixed_time_block_median_drift_and_change",
        "window_ns": list(PRIMARY_WINDOW_NS),
        "fixed_time_blocks": block_count,
        "minimum_frames_per_block": minimum_frames,
        "metric_scale_floor": float(metric_scale_floor),
        "raw_median": raw_median,
        "raw_mad": raw_mad,
        "robust_scale": robust_scale,
        "linear_slope_per_ns": slope,
        "linear_change_over_window": total_linear_change,
        "observed": observed,
        "limits": limits,
        "failed_gates": failed,
        "blocks": block_rows,
        "status": "pass" if not failed else "nonstationary",
        "all_input_points_retained": True,
    }


def continuous_true_events(
    condition: Sequence[bool],
    times_ns: Sequence[float],
    minimum_duration_ns: float,
    event_type: str,
) -> list[dict[str, Any]]:
    """Return every no-gap continuous event meeting the frozen duration."""

    flags = np.asarray(condition, dtype=bool)
    times = np.asarray(times_ns, dtype=np.float64)
    require(flags.ndim == times.ndim == 1 and len(flags) == len(times) and len(times) >= 2, "Continuous-event inputs are invalid")
    require(np.all(np.isfinite(times)) and np.all(np.diff(times) > 0.0), "Continuous-event times are invalid")
    minimum_duration = float(minimum_duration_ns)
    require(math.isfinite(minimum_duration) and minimum_duration > 0.0, "Continuous-event duration must be positive")
    events: list[dict[str, Any]] = []
    start: int | None = None
    for index in range(len(flags) + 1):
        active = bool(flags[index]) if index < len(flags) else False
        if active and start is None:
            start = index
        if not active and start is not None:
            stop = index - 1
            duration = float(times[stop] - times[start])
            if duration + 1e-12 >= minimum_duration:
                events.append({
                    "event_type": event_type,
                    "start_frame_index_zero_based": int(start),
                    "end_frame_index_zero_based": int(stop),
                    "start_time_ns": float(times[start]),
                    "end_time_ns": float(times[stop]),
                    "continuous_duration_ns": duration,
                    "frame_count": int(stop - start + 1),
                    "minimum_duration_ns": minimum_duration,
                    "gap_bridging": False,
                    "all_frames_retained": True,
                })
            start = None
    return events


def robust_first_difference(values: Sequence[float], times_ns: Sequence[float], threshold: float = 12.0) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    series = np.asarray(values, dtype=np.float64)
    times = np.asarray(times_ns, dtype=np.float64)
    require(series.ndim == times.ndim == 1 and len(series) == len(times) and len(series) >= 2, "First-difference inputs are invalid")
    require(np.all(np.isfinite(series)) and np.all(np.isfinite(times)), "First-difference inputs contain NaN or infinity")
    differences = np.diff(series)
    median = float(np.median(differences))
    mad = float(np.median(np.abs(differences - median)))
    scale = 1.4826 * mad
    flags: list[dict[str, Any]] = []
    for index, difference in enumerate(differences, start=1):
        if scale > 0.0:
            robust_z = float((difference - median) / scale)
            flagged = abs(robust_z) > threshold
            method = "median_mad"
        else:
            flagged = not math.isclose(float(difference), median, rel_tol=0.0, abs_tol=1e-15)
            robust_z = math.copysign(float("inf"), float(difference - median)) if flagged else 0.0
            method = "zero_mad_exact_difference_review"
        if flagged:
            flags.append({
                "row_index_zero_based": index,
                "time_before_ns": float(times[index - 1]),
                "time_after_ns": float(times[index]),
                "value_before": float(series[index - 1]),
                "value_after": float(series[index]),
                "first_difference": float(difference),
                "median_first_difference": median,
                "mad_first_difference": mad,
                "robust_z": "inf" if math.isinf(robust_z) and robust_z > 0 else "-inf" if math.isinf(robust_z) else robust_z,
                "method": method,
                "review_required": True,
                "point_retained": True,
            })
    return {"threshold": float(threshold), "median_first_difference": median, "mad_first_difference": mad, "robust_scale": scale, "flag_count": len(flags)}, flags


try:
    from MDAnalysis.exceptions import NoDataError
except Exception:  # pragma: no cover - MDAnalysis import is checked by callers
    class NoDataError(Exception):
        pass
