#!/usr/bin/env python3
"""Parallel-memory primary structural QC runner (three realizations concurrently).

This is an execution-strategy variant of ``run_primary_structure_memory_safe.py``
for high-core / low-memory-pressure nodes.  It launches the three realization
subprocesses concurrently instead of sequentially, but each realization still runs
the identical ``analyze_primary_structure_mdanalysis`` contract and writes its own
``structural_summary.json`` under ``<output-root>/structural_analysis/<realization_id>/``.

It does NOT change any scientific parameter: the primary window (200-500 ns), full
0-500 ns visual/QC trace, frame retention (stride=1, no deletion/smoothing/
interpolation), pocket-aligned O6U RMSD, COM, native contacts, distances/H-bonds,
TM-core/protein RMSD/RMSF, stationarity and first-difference review, and the
acceptance gates are all byte-for-byte the frozen contract. ``extension_or_recovery_window``
remains False.  Only the orchestration (concurrent vs sequential subprocesses) differs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import MDAnalysis as mda

from primary_postprocessing_common import (
    REALIZATION_IDS,
    ContractError,
    atomic_write_json,
    check_mdanalysis_version,
    require,
    sha256_file,
    validate_primary_manifest,
)
from run_primary_structure_memory_safe import (
    assemble_complete_report,
    load_realization_summaries,
    run_child,
)


def spawn_child(manifest_path: Path, output_root: Path, realization_id: str, log_directory: Path):
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
    stdout = stdout_path.open("w", encoding="utf-8")
    stderr = stderr_path.open("w", encoding="utf-8")
    process = subprocess.Popen(command, stdout=stdout, stderr=stderr, text=True)
    return process, stdout, stderr, stdout_path, stderr_path


def run_parent_parallel(manifest_path: Path, output_root: Path) -> dict[str, Any]:
    check_mdanalysis_version(mda.__version__)
    validate_primary_manifest(manifest_path)
    output_root = output_root.resolve()
    structural_directory = output_root / "structural_analysis"
    require(not structural_directory.exists(), f"Refusing to overwrite an existing structural output directory: {structural_directory}")
    structural_directory.mkdir(parents=True)
    log_directory = output_root / "logs"
    log_directory.mkdir(parents=True, exist_ok=True)

    launched: list[tuple[str, Any, Any, Any, Path, Path]] = []
    for realization_id in REALIZATION_IDS:
        redirected = spawn_child(manifest_path, output_root, realization_id, log_directory)
        launched.append((realization_id, *redirected))

    child_reports: list[dict[str, Any]] = []
    for realization_id, process, stdout, stderr, stdout_path, stderr_path in launched:
        returncode = int(process.wait())
        stdout.close()
        stderr.close()
        child_report = {
            "realization_id": realization_id,
            "returncode": returncode,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "stdout_sha256": sha256_file(stdout_path),
            "stderr_sha256": sha256_file(stderr_path),
        }
        if returncode != 0:
            child_report["status"] = "fail"
            atomic_write_json(output_root / "structural_memory_safe_runtime.json", {"status": "fail", "child_reports": child_reports + [child_report]})
            raise ContractError(f"Structural child failed for {realization_id}: return code {returncode}; see {stderr_path}")
        child_report["status"] = "pass"
        child_reports.append(child_report)

    summaries = load_realization_summaries(output_root)
    complete = assemble_complete_report(manifest_path, output_root, summaries)
    runtime = {
        "schema_version": "1.0",
        "status": complete["status"],
        "execution_strategy": "parallel_per_realization_subprocess",
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
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.child_realization:
            report = run_child(args.manifest, args.output_root, args.child_realization)
        else:
            report = run_parent_parallel(args.manifest, args.output_root)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
