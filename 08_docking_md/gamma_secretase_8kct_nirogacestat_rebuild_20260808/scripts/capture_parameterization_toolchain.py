#!/usr/bin/env python3
"""Capture a credential-free, hash-bound ligand-parameterization toolchain record."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
import pathlib
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone


PACKAGES = ("ffparam", "psi4", "openmm", "rdkit", "numpy", "scipy")


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(raw_path: str, *, required: bool = True) -> dict:
    requested = pathlib.Path(raw_path).expanduser()
    resolved = requested.resolve(strict=False)
    exists = resolved.is_file()
    if required and not exists:
        raise FileNotFoundError(f"Required file is absent: {requested}")
    return {
        "requested_path": str(requested),
        "resolved_path": str(resolved),
        "exists": exists,
        "size_bytes": resolved.stat().st_size if exists else None,
        "sha256": sha256_file(resolved) if exists else None,
    }


def executable_record(raw_path: str, version_args: list[str]) -> dict:
    found = shutil.which(raw_path) if os.sep not in raw_path else raw_path
    if not found:
        raise FileNotFoundError(f"Executable is absent: {raw_path}")
    record = file_record(found)
    command = [str(pathlib.Path(found).resolve()), *version_args]
    completed = subprocess.run(command, text=True, capture_output=True, timeout=60)
    record.update(
        {
            "version_command": command,
            "version_returncode": completed.returncode,
            "version_stdout": completed.stdout.strip(),
            "version_stderr": completed.stderr.strip(),
        }
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Version command failed: {command}")
    return record


def tree_manifest(root: pathlib.Path) -> dict:
    if not root.is_dir():
        raise FileNotFoundError(f"Source directory is absent: {root}")
    entries = []
    aggregate = hashlib.sha256()
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        relative = path.relative_to(root).as_posix()
        digest = sha256_file(path)
        size = path.stat().st_size
        entries.append({"path": relative, "size_bytes": size, "sha256": digest})
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(str(size).encode("ascii"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
    return {
        "root": str(root.resolve()),
        "file_count": len(entries),
        "tree_sha256": aggregate.hexdigest(),
        "files": entries,
    }


def memory_bytes() -> int | None:
    meminfo = pathlib.Path("/proc/meminfo")
    if not meminfo.is_file():
        return None
    for line in meminfo.read_text(encoding="utf-8").splitlines():
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) * 1024
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ffparam-archive", required=True)
    parser.add_argument("--ffparam-source-dir", required=True)
    parser.add_argument("--crest", required=True)
    parser.add_argument("--xtb", required=True)
    parser.add_argument("--gmx", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    distributions = {}
    for package in PACKAGES:
        try:
            distributions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            distributions[package] = None
    modules = {}
    for package in PACKAGES:
        module = importlib.import_module(package)
        module_file = getattr(module, "__file__", None)
        modules[package] = {
            "version_attribute": getattr(module, "__version__", None),
            "module_file": file_record(module_file) if module_file else None,
            "import_status": "pass",
        }
    if distributions["ffparam"] != "1.2.0":
        raise RuntimeError(f"Expected FFParam 1.2.0, found {distributions['ffparam']!r}")
    if modules["psi4"]["version_attribute"] != "1.9.1":
        raise RuntimeError(
            f"Expected Psi4 1.9.1, found {modules['psi4']['version_attribute']!r}"
        )
    report = {
        "schema_version": "1.0",
        "status": "pass",
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "security": "No environment variables, credentials, tokens, or command histories are recorded.",
        "host": {
            "platform": platform.platform(),
            "uname": list(platform.uname()),
            "logical_cpu_count": os.cpu_count(),
            "memory_bytes": memory_bytes(),
        },
        "python": {
            "version": sys.version,
            "executable": file_record(sys.executable),
            "distributions": distributions,
            "modules": modules,
        },
        "official_ffparam_archive": file_record(args.ffparam_archive),
        "ffparam_source_tree": tree_manifest(pathlib.Path(args.ffparam_source_dir)),
        "executables": {
            "crest": executable_record(args.crest, ["--version"]),
            "xtb": executable_record(args.xtb, ["--version"]),
            "gromacs": executable_record(args.gmx, ["--version"]),
        },
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)
    print(sha256_file(output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
