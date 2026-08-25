#!/usr/bin/env python3
"""Run the complete fail-fast server rehearsal before any real MD stage."""

from __future__ import annotations

import argparse
import compileall
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import warnings
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gmx", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    args.output_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    commands = [
        [sys.executable, "scripts/selftest_core_release_gates.py"],
        [sys.executable, "scripts/run_md_matrix.py", "--self-test"],
        [sys.executable, "scripts/validate_md_outputs.py", "--self-test"],
        [sys.executable, "scripts/hash_tree_manifest.py", "--self-test"],
        [sys.executable, "scripts/validate_qc_stationarity_report.py", "--self-test"],
        [sys.executable, "scripts/make_analysis_trajectories.py", "--self-test"],
        [sys.executable, "scripts/selftest_primary_postprocessing.py"],
        [sys.executable, "scripts/analyze_fel.py", "--self-test"],
        [sys.executable, "scripts/validate_analysis_outputs.py", "--self-test"],
    ]
    records = []
    overall_pass = True
    for index, command in enumerate(commands, 1):
        completed = subprocess.run(command, cwd=root, env=env, text=True, encoding="utf-8", errors="replace", capture_output=True, check=False)
        stem = f"{index:02d}_{Path(command[1]).stem}"
        stdout_path = args.output_dir / f"{stem}.stdout.log"
        stderr_path = args.output_dir / f"{stem}.stderr.log"
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        record = {
            "command": command,
            "returncode": completed.returncode,
            "stdout": {"path": str(stdout_path), "sha256": sha256(stdout_path), "bytes": stdout_path.stat().st_size},
            "stderr": {"path": str(stderr_path), "sha256": sha256(stderr_path), "bytes": stderr_path.stat().st_size},
        }
        records.append(record)
        if completed.returncode != 0:
            overall_pass = False
            break

    static_checks = {}
    if overall_pass:
        static_checks["compileall"] = compileall.compile_dir(str(root / "scripts"), quiet=1)
        try:
            json_files = sorted((root / "config").glob("*.json")) + sorted((root / "templates").glob("*.json"))
            for path in json_files:
                json.loads(path.read_text(encoding="utf-8"))
            static_checks["json_parse_count"] = len(json_files)
        except Exception as exc:
            static_checks["json_parse_error"] = str(exc)
            overall_pass = False

    versions = {}
    if overall_pass:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            import gemmi
            import h5py
            import MDAnalysis

        warning_text = [str(item.message) for item in captured]
        versions.update({"python": sys.version.split()[0], "gemmi": gemmi.__version__, "h5py": h5py.__version__, "MDAnalysis": MDAnalysis.__version__})
        static_checks["import_warnings"] = warning_text
        if warning_text or gemmi.__version__ != "0.7.3" or not MDAnalysis.__version__.startswith("2.10."):
            overall_pass = False
        gmx = subprocess.run([str(args.gmx), "--version"], text=True, encoding="utf-8", errors="replace", capture_output=True, check=False)
        version_match = re.search(r"GROMACS version:\s+(\S+)", gmx.stdout + gmx.stderr)
        versions["gromacs"] = version_match.group(1) if version_match else "unparsed"
        static_checks["gromacs_returncode"] = gmx.returncode
        if gmx.returncode != 0 or versions["gromacs"] != "2025.2":
            overall_pass = False

    report = {
        "schema_version": "1.0",
        "status": "pass" if overall_pass else "fail",
        "release_boundary": "A passing rehearsal validates code/runtime gates only; it does not approve ligand parameters, the CHARMM-GUI build, equilibration, canary MD, or production MD.",
        "captured_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python_executable": sys.executable,
        "gmx_executable": str(args.gmx.resolve()),
        "versions": versions,
        "commands": records,
        "static_checks": static_checks,
    }
    report_path = args.output_dir / "SERVER_REHEARSAL_REPORT.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "completed_commands": len(records), "report_sha256": sha256(report_path)}, sort_keys=True))
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
