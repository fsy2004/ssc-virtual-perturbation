#!/usr/bin/env python3
"""Freeze one common implicit-membrane slab geometry for all replicas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


REPLICAS = ("rep01", "rep02", "rep03")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def freeze_geometry(replica_rows: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    if set(replica_rows) != set(REPLICAS):
        missing = [replica for replica in REPLICAS if replica not in replica_rows]
        extra = sorted(set(replica_rows) - set(REPLICAS))
        raise ValueError(f"replica inputs differ; missing={missing}, extra={extra}")
    medians: dict[str, float] = {}
    frame_counts: dict[str, int] = {}
    for replica in REPLICAS:
        window_rows = []
        for row in replica_rows[replica]:
            time_ns = float(row["time_ns"])
            if not 200.0 <= time_ns <= 500.0:
                continue
            thickness = float(row["phosphate_peak_thickness_nm"])
            upper = float(row["upper_phosphate_peak_z_relative_nm"])
            lower = float(row["lower_phosphate_peak_z_relative_nm"])
            if not all(math.isfinite(value) for value in (time_ns, thickness, upper, lower)) or thickness <= 0:
                raise ValueError(f"{replica}: membrane geometry contains invalid values")
            if upper <= 0 or lower >= 0 or upper <= lower:
                raise ValueError(f"{replica}: leaflet identities are inverted or ambiguous")
            if int(row.get("upper_leaflet_mismatch_count", 0)) != 0 or int(row.get("lower_leaflet_mismatch_count", 0)) != 0:
                raise ValueError(f"{replica}: leaflet mismatch confounds membrane thickness")
            if int(row.get("cumulative_leaflet_flip_events", 0)) != 0:
                raise ValueError(f"{replica}: leaflet flip confounds membrane thickness")
            window_rows.append(thickness)
        if len(window_rows) < 300:
            raise ValueError(f"{replica}: fixed-window membrane rows are incomplete")
        medians[replica] = float(statistics.median(window_rows))
        frame_counts[replica] = len(window_rows)
    common_nm = float(statistics.median(medians.values()))
    return {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_membrane_geometry",
        "status": "frozen_before_endpoint_energy_results",
        "rule": "across_replica_median_of_fixed_window_replica_median_phosphate_peak_thickness",
        "window_ns": [200.0, 500.0],
        "replica_frame_counts": frame_counts,
        "replica_median_thickness_nm": medians,
        "common_thickness_nm": common_nm,
        "mthick_angstrom": common_nm * 10.0,
        "mctrdz_angstrom": 0.0,
        "membrane_normal_axis": "z",
    }


def validate_preparation_manifests(manifests: dict[str, dict[str, Any]]) -> None:
    if set(manifests) != set(REPLICAS):
        raise ValueError("preparation manifests must cover rep01-rep03")
    for replica in REPLICAS:
        payload = manifests[replica]
        if (
            payload.get("status") != "pass"
            or payload.get("membrane_normal_axis") != "z"
            or float(payload.get("membrane_midplane_z_angstrom", math.inf)) != 0.0
            or int(payload.get("frame_count", 300)) != 300
        ):
            raise ValueError(f"{replica}: preparation is not z-normal, midplane-zero, 300-frame passing input")


def _parse_bindings(values: list[str], label: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        replica, separator, raw_path = value.partition("=")
        if not separator or replica not in REPLICAS or replica in result:
            raise ValueError(f"{label}: expected unique repNN=PATH bindings")
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        result[replica] = path
    if set(result) != set(REPLICAS):
        raise ValueError(f"{label}: rep01-rep03 are required")
    return result


def _load_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membrane-csv", action="append", required=True)
    parser.add_argument("--membrane-summary", action="append", required=True)
    parser.add_argument("--preparation-manifest", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix(".json.sha256").exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    csv_paths = _parse_bindings(args.membrane_csv, "membrane-csv")
    summary_paths = _parse_bindings(args.membrane_summary, "membrane-summary")
    preparation_paths = _parse_bindings(args.preparation_manifest, "preparation-manifest")
    summaries = {replica: _load_json(path) for replica, path in summary_paths.items()}
    for replica, payload in summaries.items():
        if (
            payload.get("technical_status") != "pass"
            or payload.get("sampling_status") != "pass"
            or payload.get("preproduction_status") != "pass"
        ):
            raise ValueError(f"{replica}: complete membrane QC gates are not passing")
    preparations = {replica: _load_json(path) for replica, path in preparation_paths.items()}
    validate_preparation_manifests(preparations)
    geometry = freeze_geometry({replica: _load_csv(path) for replica, path in csv_paths.items()})
    geometry["sources"] = {
        replica: {
            "membrane_raw_csv": {"path": str(csv_paths[replica]), "sha256": sha256(csv_paths[replica])},
            "membrane_summary": {"path": str(summary_paths[replica]), "sha256": sha256(summary_paths[replica])},
            "preparation_manifest": {"path": str(preparation_paths[replica]), "sha256": sha256(preparation_paths[replica])},
        }
        for replica in REPLICAS
    }
    args.output.write_text(json.dumps(geometry, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = sha256(args.output)
    args.output.with_suffix(".json.sha256").write_text(f"{digest}  {args.output.name}\n", encoding="ascii")
    print(json.dumps({"status": "pass", "output": str(args.output), "sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
