#!/usr/bin/env python3
"""Memory-safe primary structural QC runner.

This is a technical recovery runner for low-cgroup-memory nodes.  It preserves
the exact structural analysis contract and output schema from
``analyze_primary_structure_mdanalysis.py`` while executing each realization in
its own Python process.  The parent process only validates and assembles the
three per-realization reports, so large trajectory arrays, raw rows, and
diagnostic objects are released between realizations.

It does not modify raw/PBC trajectories, change scientific thresholds, drop
frames, smooth traces, or introduce a recovery/extension analysis window.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import MDAnalysis as mda
import numpy as np

import analyze_primary_structure_mdanalysis as structural
from primary_postprocessing_common import (
    PRIMARY_WINDOW_NS,
    REALIZATION_IDS,
    ContractError,
    atomic_write_json,
    check_mdanalysis_version,
    load_json,
    require,
    resolve_record,
    sha256_file,
    validate_primary_manifest,
)


def find_realization(manifest: Mapping[str, Any], realization_id: str) -> Mapping[str, Any]:
    require(realization_id in REALIZATION_IDS, f"Unknown realization_id: {realization_id}")
    for realization in manifest["realizations"]:
        if realization.get("realization_id") == realization_id:
            return realization
    raise ContractError(f"Manifest does not contain realization_id: {realization_id}")


def load_realization_summaries(output_root: Path) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for realization_id in REALIZATION_IDS:
        summary_path = output_root / "structural_analysis" / realization_id / "structural_summary.json"
        require(summary_path.is_file() and summary_path.stat().st_size > 0, f"Missing structural summary for {realization_id}: {summary_path}")
        summary = load_json(summary_path)
        require(summary.get("realization_id") == realization_id, f"Structural summary realization mismatch: {summary_path}")
        summaries.append(summary)
    return summaries


def validate_shared_time_contract(summaries: Sequence[Mapping[str, Any]], endpoint_tolerance_ns: float) -> None:
    first = summaries[0]
    for summary in summaries[1:]:
        for key in ("input_frame_count", "primary_window_frame_count"):
            require(int(summary[key]) == int(first[key]), f"Saved frame count differs across realizations for {key}")
        for key in ("saved_step_ns", "first_time_ns", "last_time_ns"):
            require(abs(float(summary[key]) - float(first[key])) <= endpoint_tolerance_ns, f"Saved time axis differs across realizations for {key}")


def assemble_complete_report(manifest_path: Path, output_root: Path, summaries: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    manifest, manifest_base = validate_primary_manifest(manifest_path)
    mapping_path = resolve_record(manifest_base, manifest["mapping_records"]["structural"], "mapping_records.structural")
    require([summary["realization_id"] for summary in summaries] == list(REALIZATION_IDS), "Structural analysis lost or reordered a realization")
    validate_shared_time_contract(summaries, float(manifest["time_contract"]["endpoint_tolerance_ns"]))
    overall_status = "pass" if all(summary["sampling_status"] == "pass" and summary["scientific_status"] == "pass" for summary in summaries) else "inconclusive"
    complete = {
        "schema_version": "1.0",
        "status": overall_status,
        "technical_status": "pass",
        "system_id": manifest["system_id"],
        "construction_count": 1,
        "realization_ids": list(REALIZATION_IDS),
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "extension_or_recovery_window": False,
        "mdanalysis_version": mda.__version__,
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path.resolve())},
        "mapping_record": {"path": str(mapping_path), "sha256": sha256_file(mapping_path)},
        "realization_summaries": [
            {
                "realization_id": item["realization_id"],
                "technical_status": item["technical_status"],
                "sampling_status": item["sampling_status"],
                "stationarity_status": item["stationarity_status"],
                "scientific_status": item["scientific_status"],
                "scientific_failures": item["scientific_failures"],
            }
            for item in summaries
        ],
        "prohibited_outputs": {"mmgbsa": False, "mmpbsa": False, "smoothed_trace": False, "deleted_frames": False, "interpolated_frames": False},
        "execution_strategy": {
            "name": "memory_safe_per_realization_subprocess",
            "scientific_contract_changed": False,
            "raw_trajectory_policy": "read_only_no_overwrite",
            "failed_prior_output_reused_as_pass_evidence": False,
        },
    }
    atomic_write_json(output_root / "structural_analysis" / "COMPLETE.json", complete)
    return complete


def drop_linux_page_cache_if_possible() -> dict[str, Any]:
    report = {"attempted": False, "status": "not_applicable"}
    drop_path = Path("/proc/sys/vm/drop_caches")
    if os.name != "posix" or not drop_path.exists():
        return report
    report["attempted"] = True
    subprocess.run(["sync"], check=False)
    try:
        drop_path.write_text("3\n", encoding="ascii")
    except OSError as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
    else:
        report["status"] = "ok"
    return report


def run_child(manifest_path: Path, output_root: Path, realization_id: str) -> dict[str, Any]:
    check_mdanalysis_version(mda.__version__)
    manifest, manifest_base = validate_primary_manifest(manifest_path)
    realization = find_realization(manifest, realization_id)
    mapping_path = resolve_record(manifest_base, manifest["mapping_records"]["structural"], "mapping_records.structural")
    reference_topology = resolve_record(manifest_base, manifest["reference"]["topology"], "reference.topology")
    reference_coordinates = resolve_record(manifest_base, manifest["reference"]["coordinates"], "reference.coordinates")
    mapping_record = load_json(mapping_path)
    output_directory = output_root.resolve() / "structural_analysis"
    require(output_directory.is_dir(), f"Parent structural output directory does not exist: {output_directory}")
    realization_directory = output_directory / realization_id
    require(not realization_directory.exists(), f"Refusing to overwrite existing structural realization directory: {realization_directory}")
    try:
        reference = mda.Universe(str(reference_topology), str(reference_coordinates))
    except Exception as exc:
        raise ContractError(f"MDAnalysis could not open the frozen 8KCT reference: {exc}") from exc
    summary, times = structural._process_realization(manifest, manifest_base, realization, reference, mapping_record, output_directory)
    require(len(times) == int(summary["input_frame_count"]), f"{realization_id} time accounting differs")
    return {
        "schema_version": "1.0",
        "status": "pass",
        "realization_id": realization_id,
        "summary_path": str((realization_directory / "structural_summary.json").resolve()),
        "summary_sha256": sha256_file(realization_directory / "structural_summary.json"),
    }


def run_parent(manifest_path: Path, output_root: Path, drop_caches_between_realizations: bool = True) -> dict[str, Any]:
    check_mdanalysis_version(mda.__version__)
    validate_primary_manifest(manifest_path)
    output_root = output_root.resolve()
    structural_directory = output_root / "structural_analysis"
    require(not structural_directory.exists(), f"Refusing to overwrite an existing structural output directory: {structural_directory}")
    structural_directory.mkdir(parents=True)
    log_directory = output_root / "logs"
    log_directory.mkdir(parents=True, exist_ok=True)
    child_reports = []
    for realization_id in REALIZATION_IDS:
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--manifest",
            str(manifest_path.resolve()),
            "--output-root",
            str(output_root),
            "--child-realization",
            realization_id,
        ]
        stdout_path = log_directory / f"structural_{realization_id}.stdout.log"
        stderr_path = log_directory / f"structural_{realization_id}.stderr.log"
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, text=True, check=False)
        child_report = {
            "realization_id": realization_id,
            "returncode": int(completed.returncode),
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "stdout_sha256": sha256_file(stdout_path),
            "stderr_sha256": sha256_file(stderr_path),
        }
        if completed.returncode != 0:
            child_report["status"] = "fail"
            atomic_write_json(output_root / "structural_memory_safe_runtime.json", {"status": "fail", "child_reports": child_reports + [child_report]})
            raise ContractError(f"Structural child failed for {realization_id}: return code {completed.returncode}; see {stderr_path}")
        child_report["status"] = "pass"
        if drop_caches_between_realizations:
            child_report["drop_caches_after_child"] = drop_linux_page_cache_if_possible()
        child_reports.append(child_report)
    summaries = load_realization_summaries(output_root)
    complete = assemble_complete_report(manifest_path, output_root, summaries)
    runtime = {
        "schema_version": "1.0",
        "status": complete["status"],
        "execution_strategy": "memory_safe_per_realization_subprocess",
        "child_reports": child_reports,
        "complete_json": {
            "path": str((structural_directory / "COMPLETE.json").resolve()),
            "sha256": sha256_file(structural_directory / "COMPLETE.json"),
        },
    }
    atomic_write_json(output_root / "structural_memory_safe_runtime.json", runtime)
    return complete


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--child-realization", choices=REALIZATION_IDS, help=argparse.SUPPRESS)
    parser.add_argument("--no-drop-caches", action="store_true", help="Do not attempt Linux page-cache drop between realization subprocesses")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.child_realization:
            report = run_child(args.manifest, args.output_root, args.child_realization)
        else:
            report = run_parent(args.manifest, args.output_root, drop_caches_between_realizations=not args.no_drop_caches)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
