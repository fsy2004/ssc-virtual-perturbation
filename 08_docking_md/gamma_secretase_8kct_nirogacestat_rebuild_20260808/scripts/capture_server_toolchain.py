#!/usr/bin/env python3
"""Capture a secret-free, hash-verifiable server toolchain record."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def execute(argv: list[str], cwd: Path, name: str) -> dict:
    completed = subprocess.run(argv, cwd=cwd, text=True, capture_output=True)
    stdout = cwd / f"{name}.stdout.txt"
    stderr = cwd / f"{name}.stderr.txt"
    stdout.write_text(completed.stdout, encoding="utf-8", newline="\n")
    stderr.write_text(completed.stderr, encoding="utf-8", newline="\n")
    record = {
        "argv": argv,
        "returncode": completed.returncode,
        "stdout": stdout.name,
        "stdout_sha256": sha256(stdout),
        "stderr": stderr.name,
        "stderr_sha256": sha256(stderr),
    }
    if completed.returncode != 0:
        raise RuntimeError(json.dumps(record, indent=2))
    return record


def module_versions() -> dict[str, str]:
    import gemmi
    import h5py
    import MDAnalysis
    import numpy
    import pandas
    import scipy

    return {
        "gemmi": gemmi.__version__,
        "h5py": h5py.__version__,
        "hdf5_runtime": h5py.version.hdf5_version,
        "MDAnalysis": MDAnalysis.__version__,
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--conda", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--nvidia-smi", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    paths = {key: Path(value).resolve() for key, value in {
        "conda": args.conda,
        "prefix": args.prefix,
        "gmx": args.gmx,
        "nvidia_smi": args.nvidia_smi,
    }.items()}
    for key in ("conda", "gmx", "nvidia_smi"):
        if not paths[key].is_file() or not os.access(paths[key], os.X_OK):
            raise SystemExit(f"Unavailable executable {key}: {paths[key]}")
    if not paths["prefix"].is_dir():
        raise SystemExit(f"Unavailable environment prefix: {paths['prefix']}")
    if Path(sys.prefix).resolve() != paths["prefix"]:
        raise SystemExit(f"Run with the recorded environment Python: {sys.prefix} != {paths['prefix']}")

    out = Path(args.output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)

    commands = [
        execute([str(paths["conda"]), "list", "--prefix", str(paths["prefix"]), "--json"], out, "conda_list_json"),
        execute([str(paths["conda"]), "list", "--prefix", str(paths["prefix"]), "--explicit"], out, "conda_explicit"),
        execute([str(paths["gmx"]), "--version"], out, "gromacs_version"),
        execute([str(paths["nvidia_smi"]), "--query-gpu=name,uuid,driver_version,memory.total,compute_cap", "--format=csv,noheader,nounits"], out, "nvidia_smi_query"),
        execute(["ldd", str(paths["gmx"])], out, "gromacs_ldd"),
        execute(["uname", "-a"], out, "uname"),
    ]
    os_release_source = Path("/etc/os-release")
    if not os_release_source.is_file():
        raise SystemExit("Missing /etc/os-release")
    os_release_copy = out / "os-release.txt"
    os_release_copy.write_bytes(os_release_source.read_bytes())
    safe_environment = {
        key: os.environ.get(key)
        for key in ["CUDA_VISIBLE_DEVICES", "GMXBIN", "GMXDATA", "GMXLDLIB", "OMP_NUM_THREADS"]
        if os.environ.get(key) is not None
    }
    executables = {
        "conda": paths["conda"],
        "gmx": paths["gmx"],
        "nvidia_smi": paths["nvidia_smi"],
        "python": Path(sys.executable).resolve(),
    }
    libgromacs = paths["gmx"].parent.parent / "lib" / "libgromacs.so.10"
    if libgromacs.is_file():
        executables["libgromacs"] = libgromacs.resolve()
    executable_records = {
        key: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for key, path in executables.items()
    }
    version_text = (out / "gromacs_version.stdout.txt").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"^GROMACS version:\s*(\S+)", version_text, re.MULTILINE)
    if version_match is None:
        raise SystemExit("Cannot parse GROMACS version for executable identity")
    linked_libraries = []
    if "libgromacs" in executable_records:
        linked_libraries.append({"role": "libgromacs", **executable_records["libgromacs"]})
    gmx_identity = {
        "schema_version": "1.0",
        "record_type": "gromacs_executable_identity",
        "resolved_path": executable_records["gmx"]["path"],
        "bytes": executable_records["gmx"]["bytes"],
        "sha256": executable_records["gmx"]["sha256"],
        "gromacs_version": version_match.group(1),
        "linked_libraries": linked_libraries,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "reviewer": "pending formal approved-build environment review",
    }
    identity_path = out / "gromacs_executable_identity.json"
    identity_path.write_text(json.dumps(gmx_identity, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    artifacts = {}
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "SERVER_TOOLCHAIN_RECORD.json":
            artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    report = {
        "schema_version": 1,
        "status": "pass",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "platform": platform.platform(),
        "python_version": sys.version,
        "python_executable": str(Path(sys.executable).resolve()),
        "environment_prefix": str(paths["prefix"]),
        "module_versions": module_versions(),
        "executables": executable_records,
        "safe_environment": safe_environment,
        "commands": commands,
        "artifacts": artifacts,
    }
    report_path = out / "SERVER_TOOLCHAIN_RECORD.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "report": str(report_path), "sha256": sha256(report_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
