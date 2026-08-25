#!/usr/bin/env python3
"""Extract and parse exact GROMACS EDR terms without altering any data point."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from primary_postprocessing_common import (
    PRIMARY_WINDOW_NS,
    REALIZATION_IDS,
    ENERGY_STATIONARITY_METRICS,
    ContractError,
    atomic_write_csv,
    atomic_write_json,
    block_diagnostics,
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


REQUIRED_TERM_SPECS = (
    ("temperature_k", "Temperature", "K"),
    ("pressure_bar", "Pressure", "bar"),
    ("pressure_xx_bar", "Pres-XX", "bar"),
    ("pressure_yy_bar", "Pres-YY", "bar"),
    ("pressure_zz_bar", "Pres-ZZ", "bar"),
    ("pressure_xy_bar", "Pres-XY", "bar"),
    ("pressure_xz_bar", "Pres-XZ", "bar"),
    ("pressure_yz_bar", "Pres-YZ", "bar"),
    ("potential_energy_kj_mol", "Potential", "kJ/mol"),
    ("kinetic_energy_kj_mol", "Kinetic-En.", "kJ/mol"),
    ("total_energy_kj_mol", "Total-Energy", "kJ/mol"),
    ("density_kg_m3", "Density", "kg/m^3"),
    ("volume_nm3", "Volume", "nm^3"),
    ("box_x_nm", "Box-X", "nm"),
    ("box_y_nm", "Box-Y", "nm"),
    ("box_z_nm", "Box-Z", "nm"),
)
REQUIRED_TERM_KEYS = tuple(item[0] for item in REQUIRED_TERM_SPECS)


def validate_terms_record(record: Mapping[str, Any], allow_synthetic: bool = False) -> list[dict[str, Any]]:
    require(record.get("schema_version") == "1.0", "Energy-term schema_version must be 1.0")
    allowed = {"approved"}
    if allow_synthetic:
        allowed.add("synthetic_self_test")
    require(record.get("approval_status") in allowed, "Energy-term record is not approved")
    require(not has_placeholder(record), "Energy-term record contains TODO/REPLACE_ME placeholders")
    require(record.get("gromacs_time_unit") == "ps", "GROMACS energy time must be parsed as ps")
    terms = record.get("terms")
    require(isinstance(terms, list) and len(terms) == len(REQUIRED_TERM_KEYS), "Energy record must contain exactly the required terms")
    observed_specs = tuple((item.get("key"), item.get("gmx_name"), item.get("unit")) for item in terms)
    require(observed_specs == REQUIRED_TERM_SPECS, "Energy keys, exact GROMACS names, units, or order differ from the frozen specification")
    require(tuple(item[0] for item in REQUIRED_TERM_SPECS) == ENERGY_STATIONARITY_METRICS, "Every exact energy term must have a mandatory stationarity gate")
    require(all(term.get("required") is True for term in terms), "Every frozen energy term must be required")
    return terms


def parse_xvg(path: Path, expected_legend: str) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    require(path.is_file() and path.stat().st_size > 0, f"Missing/empty XVG: {path}")
    times = []
    values = []
    legends = []
    metadata_lines = 0
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(("#", "@")):
                metadata_lines += 1
                if "legend" in line and '"' in line:
                    legends.append(line.split('"', 1)[1].rsplit('"', 1)[0])
                continue
            fields = line.split()
            require(len(fields) == 2, f"{path.name}:{line_number} must contain exactly time plus one separately extracted term")
            try:
                time_value, energy_value = (float(fields[0]), float(fields[1]))
            except ValueError as exc:
                raise ContractError(f"{path.name}:{line_number} contains a nonnumeric value") from exc
            require(math.isfinite(time_value) and math.isfinite(energy_value), f"{path.name}:{line_number} contains NaN or infinity")
            times.append(time_value / 1000.0)
            values.append(energy_value)
    require(len(times) >= 2, f"{path.name} has fewer than two data rows")
    normalized_expected = re.sub(r"[-\s]+", "", expected_legend)
    matching_legends = [legend for legend in legends if re.sub(r"[-\s]+", "", legend) == normalized_expected]
    require(
        bool(matching_legends),
        f"{path.name} does not declare the selected term legend {expected_legend!r} allowing only GROMACS space/hyphen label variants; observed {legends}",
    )
    return np.asarray(times, dtype=np.float64), np.asarray(values, dtype=np.float64), {
        "metadata_lines": metadata_lines,
        "legends": legends,
        "matched_legend": matching_legends[0],
        "expected_legend": expected_legend,
        "legend_match_policy": "exact_after_space_hyphen_normalization_only",
        "data_rows": len(times),
    }


def _extract_term(gmx_executable: str, edr_path: Path, term: Mapping[str, Any], destination: Path) -> dict[str, Any]:
    key = term["key"]
    destination.mkdir(parents=True, exist_ok=True)
    selection_path = destination / f"selection__{key}.txt"
    xvg_path = destination / f"raw__{key}.xvg"
    log_path = destination / f"gmx_energy__{key}.log"
    command_path = destination / f"command__{key}.json"
    selection_text = f"{term['gmx_name']}\n0\n"
    selection_path.write_text(selection_text, encoding="utf-8", newline="\n")
    command = [gmx_executable, "energy", "-f", str(edr_path), "-o", str(xvg_path)]
    atomic_write_json(command_path, {"argv": command, "stdin_file": selection_path.name, "stdin_sha256": sha256_file(selection_path)})
    result = subprocess.run(command, input=selection_text, text=True, capture_output=True, check=False)
    log_path.write_text(
        f"returncode={result.returncode}\n--- STDOUT ---\n{result.stdout}\n--- STDERR ---\n{result.stderr}\n",
        encoding="utf-8",
        newline="\n",
    )
    require(result.returncode == 0, f"gmx energy failed for {key}; see {log_path}")
    require(xvg_path.is_file() and xvg_path.stat().st_size > 0, f"gmx energy produced no XVG for {key}")
    return {
        "selection": selection_path,
        "xvg": xvg_path,
        "log": log_path,
        "command": command_path,
    }


def _verify_gromacs_version(gmx_executable: str, expected_version: str, output_directory: Path) -> dict[str, Any]:
    require(isinstance(expected_version, str) and expected_version and not has_placeholder(expected_version), "required_gromacs_version is not frozen")
    result = subprocess.run([gmx_executable, "--version"], text=True, capture_output=True, check=False)
    log_path = output_directory / "gromacs_version.log"
    log_path.write_text(
        f"returncode={result.returncode}\n--- STDOUT ---\n{result.stdout}\n--- STDERR ---\n{result.stderr}\n",
        encoding="utf-8",
        newline="\n",
    )
    require(result.returncode == 0, f"Cannot query GROMACS version; see {log_path}")
    combined = f"{result.stdout}\n{result.stderr}"
    match = re.search(r"^\s*GROMACS version:\s*(\S+)\s*$", combined, flags=re.IGNORECASE | re.MULTILINE)
    require(match is not None, f"Cannot parse GROMACS version; see {log_path}")
    observed_version = str(match.group(1))
    require(observed_version == expected_version, f"GROMACS version differs: expected {expected_version}, observed {observed_version}")
    return {"expected": expected_version, "observed": observed_version, "log": str(log_path), "log_sha256": sha256_file(log_path)}


def _copy_existing_term(source_root: Path, term: Mapping[str, Any], destination: Path) -> dict[str, Any]:
    key = term["key"]
    destination.mkdir(parents=True, exist_ok=True)
    source_xvg = source_root / f"raw__{key}.xvg"
    source_log = source_root / f"gmx_energy__{key}.log"
    source_selection = source_root / f"selection__{key}.txt"
    source_command = source_root / f"command__{key}.json"
    for path in (source_xvg, source_log, source_selection, source_command):
        require(path.is_file() and path.stat().st_size > 0, f"Existing extraction artifact is missing: {path}")
    outputs = {
        "xvg": destination / source_xvg.name,
        "log": destination / source_log.name,
        "selection": destination / source_selection.name,
        "command": destination / source_command.name,
    }
    for label, output in outputs.items():
        source = {"xvg": source_xvg, "log": source_log, "selection": source_selection, "command": source_command}[label]
        shutil.copy2(source, output)
        require(sha256_file(source) == sha256_file(output), f"Copy hash mismatch for {source}")
    return outputs


def _parse_realization(
    manifest: Mapping[str, Any],
    manifest_base: Path,
    realization: Mapping[str, Any],
    terms: list[dict[str, Any]],
    output_directory: Path,
    mode: str,
    existing_xvg_root: Path | None,
) -> tuple[dict[str, Any], np.ndarray]:
    realization_id = str(realization["realization_id"])
    realization_directory = output_directory / realization_id
    raw_directory = realization_directory / "raw_gmx_energy"
    if mode == "extract":
        require(manifest.get("energy_execution", {}).get("server_execution_authorized") is True, "Real EDR extraction is server-gated and not authorized in the manifest")
        edr_path = resolve_record(manifest_base, realization["energy_edr"], f"{realization_id}.energy_edr")
        production_log = resolve_record(manifest_base, realization["production_log"], f"{realization_id}.production_log")
        gmx_executable = str(manifest["energy_execution"]["gmx_executable"])
        artifacts = {term["key"]: _extract_term(gmx_executable, edr_path, term, raw_directory) for term in terms}
        copied_log = realization_directory / "source_production.log"
        shutil.copy2(production_log, copied_log)
        require(sha256_file(production_log) == sha256_file(copied_log), "Production-log copy hash mismatch")
        source_records = {
            "energy_edr": {"path": str(edr_path), "sha256": sha256_file(edr_path)},
            "production_log": {"path": str(production_log), "sha256": sha256_file(production_log)},
            "copied_production_log": {"path": str(copied_log), "sha256": sha256_file(copied_log)},
        }
    else:
        require(existing_xvg_root is not None, "--existing-xvg-root is required in parse-existing mode")
        source_directory = existing_xvg_root.resolve() / realization_id
        artifacts = {term["key"]: _copy_existing_term(source_directory, term, raw_directory) for term in terms}
        source_records = {"existing_extraction_directory": str(source_directory)}

    combined: dict[str, np.ndarray] = {}
    metadata: dict[str, Any] = {}
    shared_times: np.ndarray | None = None
    endpoint_tolerance = float(manifest["time_contract"]["endpoint_tolerance_ns"])
    for term in terms:
        times, values, xvg_metadata = parse_xvg(artifacts[term["key"]]["xvg"], str(term["gmx_name"]))
        validate_time_axis(times, endpoint_tolerance)
        if shared_times is None:
            shared_times = times
        else:
            require(np.allclose(times, shared_times, rtol=0.0, atol=endpoint_tolerance), f"{realization_id} energy term time arrays differ")
        combined[term["key"]] = values
        metadata[term["key"]] = {
            **xvg_metadata,
            "gmx_name": term["gmx_name"],
            "unit": term["unit"],
            "xvg_sha256": sha256_file(artifacts[term["key"]]["xvg"]),
            "log_sha256": sha256_file(artifacts[term["key"]]["log"]),
            "selection_sha256": sha256_file(artifacts[term["key"]]["selection"]),
            "command_sha256": sha256_file(artifacts[term["key"]]["command"]),
        }
    assert shared_times is not None
    rows = []
    for index, time_ns in enumerate(shared_times):
        row: dict[str, Any] = {
            "system_id": manifest["system_id"],
            "realization_id": realization_id,
            "row_index_zero_based": index,
            "time_ns": float(time_ns),
            "in_primary_window_200_500_ns": int(200.0 - endpoint_tolerance <= time_ns <= 500.0 + endpoint_tolerance),
        }
        for term in terms:
            row[term["key"]] = float(combined[term["key"]][index])
        rows.append(row)
    raw_path = realization_directory / "energy_raw_unsmoothed.csv"
    raw_fields = ["system_id", "realization_id", "row_index_zero_based", "time_ns", "in_primary_window_200_500_ns", *REQUIRED_TERM_KEYS]
    atomic_write_csv(raw_path, raw_fields, rows)

    flag_rows = []
    first_difference_summaries = {}
    threshold = float(manifest["diagnostics"]["robust_first_difference_z_threshold"])
    for term in terms:
        summary, flags = robust_first_difference(combined[term["key"]], shared_times, threshold)
        first_difference_summaries[term["key"]] = summary
        for flag in flags:
            flag_rows.append({
                "system_id": manifest["system_id"],
                "realization_id": realization_id,
                "term_key": term["key"],
                "gmx_name": term["gmx_name"],
                "unit": term["unit"],
                **flag,
            })
    flags_path = realization_directory / "energy_first_difference_review_flags.csv"
    flag_fields = ["system_id", "realization_id", "term_key", "gmx_name", "unit", "row_index_zero_based", "time_before_ns", "time_after_ns", "value_before", "value_after", "first_difference", "median_first_difference", "mad_first_difference", "robust_z", "method", "review_required", "point_retained"]
    atomic_write_csv(flags_path, flag_fields, flag_rows)
    thermodynamic_gates = manifest["acceptance_gates"]["thermodynamic_cell_qc"]
    primary_mask = primary_window_mask(shared_times, endpoint_tolerance)
    diagnostic_payload: dict[str, Any] = {}
    block_rows: list[dict[str, Any]] = []
    sampling_failures: list[str] = []
    stationarity_failures: list[str] = []
    for term in terms:
        key = str(term["key"])
        diagnostic, blocks = block_diagnostics(combined[key][primary_mask], shared_times[primary_mask], manifest["diagnostics"])
        diagnostic["stationarity"] = stationarity_diagnostics(
            combined[key][primary_mask],
            shared_times[primary_mask],
            manifest["diagnostics"]["stationarity"],
            float(manifest["acceptance_gates"]["stationarity_scale_floors"]["energy"][key]),
        )
        diagnostic_payload[key] = diagnostic
        if diagnostic["status"] != "pass":
            sampling_failures.append(key)
        if diagnostic["stationarity"]["status"] != "pass":
            stationarity_failures.append(key)
        for block in blocks:
            block_rows.append({"system_id": manifest["system_id"], "realization_id": realization_id, "term_key": key, "window": "primary_200_500_ns", **block})
    block_path = realization_directory / "energy_block_summaries.csv"
    atomic_write_csv(
        block_path,
        ["system_id", "realization_id", "term_key", "window", "block_index_zero_based", "start_time_ns", "end_time_ns", "frame_count", "mean", "median", "minimum", "maximum"],
        block_rows,
    )
    temperature_mean = float(np.mean(combined["temperature_k"][primary_mask]))
    pressure_mean = float(np.mean(combined["pressure_bar"][primary_mask]))
    density_mean = float(np.mean(combined["density_kg_m3"][primary_mask]))
    energy_denominator = np.maximum.reduce([
        np.ones(len(shared_times), dtype=np.float64),
        np.abs(combined["total_energy_kj_mol"]),
        np.abs(combined["potential_energy_kj_mol"]) + np.abs(combined["kinetic_energy_kj_mol"]),
    ])
    energy_closure = np.abs(combined["total_energy_kj_mol"] - combined["potential_energy_kj_mol"] - combined["kinetic_energy_kj_mol"]) / energy_denominator
    pressure_trace_closure = np.abs(combined["pressure_bar"] - (combined["pressure_xx_bar"] + combined["pressure_yy_bar"] + combined["pressure_zz_bar"]) / 3.0)
    box_volume = combined["box_x_nm"] * combined["box_y_nm"] * combined["box_z_nm"]
    volume_closure = np.abs(combined["volume_nm3"] - box_volume) / np.maximum(1.0, np.abs(combined["volume_nm3"]))
    scientific_failures = []
    if np.any(combined["temperature_k"] <= 0.0) or np.any(combined["kinetic_energy_kj_mol"] < 0.0) or np.any(combined["density_kg_m3"] <= 0.0) or np.any(combined["volume_nm3"] <= 0.0) or any(np.any(combined[key] <= 0.0) for key in ("box_x_nm", "box_y_nm", "box_z_nm")):
        scientific_failures.append("nonphysical_positive_quantity_gate")
    if abs(temperature_mean - float(thermodynamic_gates["target_temperature_k"])) > float(thermodynamic_gates["maximum_absolute_primary_mean_temperature_deviation_k"]):
        scientific_failures.append("primary_mean_temperature_gate")
    pressure_range = [float(value) for value in thermodynamic_gates["approved_primary_mean_pressure_range_bar"]]
    if not pressure_range[0] <= pressure_mean <= pressure_range[1]:
        scientific_failures.append("primary_mean_pressure_gate")
    density_range = [float(value) for value in thermodynamic_gates["approved_primary_mean_density_range_kg_m3"]]
    if not density_range[0] <= density_mean <= density_range[1]:
        scientific_failures.append("primary_mean_density_gate")
    if float(np.max(energy_closure)) > float(thermodynamic_gates["maximum_relative_total_energy_closure_error"]):
        scientific_failures.append("total_energy_closure_gate")
    if float(np.max(pressure_trace_closure)) > float(thermodynamic_gates["maximum_absolute_pressure_trace_closure_bar"]):
        scientific_failures.append("pressure_trace_closure_gate")
    if float(np.max(volume_closure)) > float(thermodynamic_gates["maximum_relative_orthorhombic_volume_closure_error"]):
        scientific_failures.append("orthorhombic_volume_closure_gate")
    thermodynamic_acceptance = {
        "criteria": dict(thermodynamic_gates),
        "primary_mean_temperature_k": temperature_mean,
        "primary_mean_pressure_bar": pressure_mean,
        "primary_mean_density_kg_m3": density_mean,
        "maximum_relative_total_energy_closure_error": float(np.max(energy_closure)),
        "maximum_absolute_pressure_trace_closure_bar": float(np.max(pressure_trace_closure)),
        "maximum_relative_orthorhombic_volume_closure_error": float(np.max(volume_closure)),
        "status": "pass" if not scientific_failures else "scientific_fail",
        "failure_triggers_rerun_or_extension": False,
    }
    summary = {
        "schema_version": "1.0",
        "analysis": "gromacs_edr_exact_terms_qc",
        "system_id": manifest["system_id"],
        "realization_id": realization_id,
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "input_row_count": len(shared_times),
        "first_time_ns": float(shared_times[0]),
        "last_time_ns": float(shared_times[-1]),
        "saved_step_ns": float(np.median(np.diff(shared_times))),
        "technical_status": "pass",
        "sampling_status": "pass" if not sampling_failures and not stationarity_failures else "inconclusive",
        "sampling_failures": sampling_failures,
        "stationarity_status": "pass" if not stationarity_failures else "fail",
        "stationarity_failures": stationarity_failures,
        "scientific_status": "pass" if not scientific_failures else "fail",
        "scientific_failures": scientific_failures,
        "thermodynamic_acceptance": thermodynamic_acceptance,
        "scientific_failure_triggers_rerun_or_extension": False,
        "review_status": "pending" if flag_rows else "no_flags",
        "review_flag_count": len(flag_rows),
        "flagged_points_retained": True,
        "smoothing": False,
        "frame_deletion": False,
        "interpolation": False,
        "source_records": source_records,
        "terms": metadata,
        "first_difference_diagnostics": first_difference_summaries,
        "primary_diagnostics": diagnostic_payload,
        "outputs": {
            raw_path.name: {"bytes": raw_path.stat().st_size, "sha256": sha256_file(raw_path)},
            flags_path.name: {"bytes": flags_path.stat().st_size, "sha256": sha256_file(flags_path)},
            block_path.name: {"bytes": block_path.stat().st_size, "sha256": sha256_file(block_path)},
        },
    }
    summary_path = realization_directory / "energy_summary.json"
    atomic_write_json(summary_path, summary)
    return summary, shared_times


def run(
    manifest_path: Path,
    output_root: Path,
    mode: str,
    existing_xvg_root: Path | None = None,
    allow_synthetic: bool = False,
) -> dict[str, Any]:
    require(mode == "extract" or allow_synthetic, "parse-existing is restricted to synthetic self-tests; production must extract directly from each manifest-hashed EDR on the gated server")
    manifest, manifest_base = validate_primary_manifest(manifest_path, allow_synthetic=allow_synthetic)
    terms_path = resolve_record(manifest_base, manifest["mapping_records"]["gromacs_energy_terms"], "mapping_records.gromacs_energy_terms")
    terms_record = load_json(terms_path)
    terms = validate_terms_record(terms_record, allow_synthetic=allow_synthetic)
    output_directory = output_root.resolve() / "energy_qc"
    require(not output_directory.exists(), f"Refusing to overwrite an existing energy output directory: {output_directory}")
    output_directory.mkdir(parents=True)
    if mode == "extract":
        gromacs_version = _verify_gromacs_version(
            str(manifest["energy_execution"]["gmx_executable"]),
            str(manifest["energy_execution"]["required_gromacs_version"]),
            output_directory,
        )
    else:
        gromacs_version = {"status": "synthetic_self_test_only"}
    summaries = []
    shared_times: np.ndarray | None = None
    for realization in manifest["realizations"]:
        summary, times = _parse_realization(manifest, manifest_base, realization, terms, output_directory, mode, existing_xvg_root)
        if shared_times is None:
            shared_times = times
        else:
            require(np.allclose(times, shared_times, rtol=0.0, atol=float(manifest["time_contract"]["endpoint_tolerance_ns"])), "Energy time arrays differ among realizations")
        summaries.append(summary)
    require([summary["realization_id"] for summary in summaries] == list(REALIZATION_IDS), "Energy QC lost or reordered a realization")
    scientific_failure = any(summary["scientific_status"] != "pass" for summary in summaries)
    sampling_failure = any(summary["sampling_status"] != "pass" for summary in summaries)
    complete = {
        "schema_version": "1.0",
        "status": "inconclusive" if scientific_failure or sampling_failure else "pass_with_review_flags" if any(summary["review_flag_count"] for summary in summaries) else "pass",
        "technical_status": "pass",
        "sampling_status": "inconclusive" if sampling_failure else "pass",
        "scientific_status": "fail" if scientific_failure else "pass",
        "review_required": any(summary["review_flag_count"] for summary in summaries),
        "system_id": manifest["system_id"],
        "realization_ids": list(REALIZATION_IDS),
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "extension_or_recovery_window": False,
        "exact_required_term_keys": list(REQUIRED_TERM_KEYS),
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path.resolve())},
        "energy_terms_record": {"path": str(terms_path), "sha256": sha256_file(terms_path)},
        "gromacs_version": gromacs_version,
        "realization_summaries": [
            {"realization_id": item["realization_id"], "technical_status": item["technical_status"], "sampling_status": item["sampling_status"], "stationarity_status": item["stationarity_status"], "scientific_status": item["scientific_status"], "scientific_failures": item["scientific_failures"], "review_status": item["review_status"], "review_flag_count": item["review_flag_count"]}
            for item in summaries
        ],
        "flagged_points_retained": True,
    }
    atomic_write_json(output_directory / "COMPLETE.json", complete)
    return complete


def synthetic_xvg_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="gmx_energy_xvg_selftest_") as temporary:
        root = Path(temporary)
        path = root / "raw__pressure_bar.xvg"
        lines = [
            '# synthetic XVG for parser self-test',
            '@ s0 legend "Pressure"',
        ]
        values = np.zeros(501, dtype=np.float64)
        values[250] = 1000.0
        for index, value in enumerate(values):
            lines.append(f"{index * 1000.0:.1f} {value:.8f}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        times, parsed, metadata = parse_xvg(path, "Pressure")
        validate_time_axis(times, 1e-6)
        diagnostics, flags = robust_first_difference(parsed, times, 12.0)
        require(metadata["data_rows"] == 501, "Synthetic XVG row count differs")
        require(len(flags) == 2 and diagnostics["flag_count"] == 2, "Synthetic pressure spike was not flagged on both edges")
        require(parsed[250] == 1000.0 and len(parsed) == 501, "Synthetic flagged point was deleted or changed")
        kinetic_alias_path = root / "raw__kinetic_energy_kj_mol.xvg"
        kinetic_alias_path.write_text(
            '# synthetic XVG for GROMACS label-variant self-test\n@ s0 legend "Kinetic En."\n0.0 1.0\n1000.0 1.0\n',
            encoding="utf-8",
            newline="\n",
        )
        _, _, alias_metadata = parse_xvg(kinetic_alias_path, "Kinetic-En.")
        require(alias_metadata["matched_legend"] == "Kinetic En.", "GROMACS kinetic-energy space/hyphen label variant was not retained")
    print("SELF-TEST PASS: exact/space-hyphen GROMACS legends, finite 0-500 ns rows, robust first-difference flags, and point retention.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--mode", choices=("extract", "parse-existing"), default="extract")
    parser.add_argument("--existing-xvg-root", type=Path)
    parser.add_argument("--allow-synthetic", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--self-test", action="store_true", help="Run a synthetic XVG parser/flag-retention test")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.self_test:
            synthetic_xvg_self_test()
            return 0
        require(args.manifest is not None and args.output_root is not None, "--manifest and --output-root are required")
        report = run(args.manifest, args.output_root, args.mode, args.existing_xvg_root, allow_synthetic=args.allow_synthetic)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
