#!/usr/bin/env python3
"""Capture and validate a sealed isolated gmx_MMPBSA CPU toolchain."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    from run_gmx_mmpbsa_canary import (
        GMX_MMPBSA_COMMIT,
        REQUIRED_EXECUTABLES,
        validate_toolchain_record,
    )
except ModuleNotFoundError:
    contract_path = Path(__file__).with_name("run_gmx_mmpbsa_canary.py")
    contract_spec = importlib.util.spec_from_file_location("run_gmx_mmpbsa_canary", contract_path)
    if contract_spec is None or contract_spec.loader is None:
        raise
    contract_module = importlib.util.module_from_spec(contract_spec)
    contract_spec.loader.exec_module(contract_module)
    GMX_MMPBSA_COMMIT = contract_module.GMX_MMPBSA_COMMIT
    REQUIRED_EXECUTABLES = contract_module.REQUIRED_EXECUTABLES
    validate_toolchain_record = contract_module.validate_toolchain_record


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def index_conda_packages(records: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for record in records:
        name = str(record.get("name", "")).lower()
        if not name:
            raise ValueError("conda package record has no name")
        if name in indexed:
            raise ValueError(f"duplicate conda package name: {name}")
        indexed[name] = {
            "version": str(record.get("version", "")),
            "build_string": str(record.get("build_string", record.get("build", ""))),
            "channel": str(record.get("channel", record.get("schannel", ""))),
            "url": str(record.get("url", "")),
        }
    return indexed


def _read_conda_records(prefix: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    metadata = prefix / "conda-meta"
    if not metadata.is_dir():
        raise ValueError(f"missing conda metadata: {metadata}")
    records = []
    for path in sorted(metadata.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        record["metadata_file"] = path.name
        record["metadata_sha256"] = file_sha256(path)
        records.append(record)
    return index_conda_packages(records), records


def _run(command: list[str], env: dict[str, str], timeout: int = 120) -> str:
    completed = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        timeout=timeout,
    )
    return completed.stdout


def _pip_versions(python: Path, env: dict[str, str]) -> tuple[dict[str, str], str]:
    code = (
        "import importlib.metadata as m,json; "
        "print(json.dumps({d.metadata['Name'].lower():d.version for d in m.distributions()},sort_keys=True))"
    )
    versions = json.loads(_run([str(python), "-c", code], env))
    freeze = _run([str(python), "-m", "pip", "freeze", "--all"], env)
    return versions, freeze


def _package_version(name: str, conda: dict[str, dict[str, Any]], pip: dict[str, str]) -> str:
    key = name.lower()
    if key in pip:
        return pip[key]
    if key in conda:
        return conda[key]["version"]
    return "missing"


def capture(prefix: Path, output: Path) -> dict[str, Any]:
    prefix = prefix.resolve()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    bin_dir = prefix / "bin"
    python = bin_dir / "python"
    if not python.is_file():
        raise ValueError(f"missing isolated Python: {python}")
    marker = prefix / "GMX_MMPBSA_SOURCE_COMMIT.txt"
    if not marker.is_file() or marker.read_text(encoding="ascii").strip() != GMX_MMPBSA_COMMIT:
        raise ValueError("gmx_MMPBSA source commit marker is missing or invalid")
    env = {
        "PATH": str(bin_dir) + os.pathsep + os.environ.get("PATH", ""),
        "CONDA_PREFIX": str(prefix),
        "LC_ALL": "C",
        "LANG": "C",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    conda, raw_conda = _read_conda_records(prefix)
    pip, pip_freeze = _pip_versions(python, env)
    python_version = _run([str(python), "-c", "import platform; print(platform.python_version())"], env).strip()
    executable_records = {}
    command_outputs = {}
    version_arguments = {
        "gmx_MMPBSA": ["--version"],
        "gmx": ["--version"],
        "mpirun": ["--version"],
        "sander": ["-h"],
        "cpptraj": ["--version"],
        "tleap": ["-h"],
    }
    for name in REQUIRED_EXECUTABLES:
        executable = shutil.which(name, path=str(bin_dir))
        if not executable:
            raise ValueError(f"missing required executable: {name}")
        path = Path(executable).resolve()
        executable_records[name] = {"path": str(path), "sha256": file_sha256(path)}
        command_outputs[name] = _run([str(path), *version_arguments[name]], env)
    toolchain = {
        "gmx_mmpbsa": _package_version("gmx-mmpbsa", conda, pip),
        "gmx_mmpbsa_git_commit": marker.read_text(encoding="ascii").strip(),
        "python": python_version,
        "ambertools": _package_version("ambertools", conda, pip),
        "gromacs": _package_version("gromacs", conda, pip),
        "openmpi": _package_version("openmpi", conda, pip),
        "mpi4py": _package_version("mpi4py", conda, pip),
        "numpy": _package_version("numpy", conda, pip),
        "pandas": _package_version("pandas", conda, pip),
        "matplotlib": _package_version("matplotlib", conda, pip),
        "seaborn": _package_version("seaborn", conda, pip),
        "scipy": _package_version("scipy", conda, pip),
        "tqdm": _package_version("tqdm", conda, pip),
        "parmed": _package_version("parmed", conda, pip),
        "gpu_required": False,
        "executables": executable_records,
    }
    validation = validate_toolchain_record(toolchain)
    report = {
        "schema_version": "1.0",
        "report_type": "gmx_mmpbsa_cpu_toolchain",
        "status": validation["status"],
        "prefix": str(prefix),
        "toolchain": toolchain,
        "conda_packages": conda,
        "conda_metadata_records": raw_conda,
        "pip_freeze": pip_freeze.splitlines(),
        "version_command_outputs": command_outputs,
        "thread_environment": {key: env[key] for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(output)
    output.with_suffix(output.suffix + ".sha256").write_text(f"{digest}  {output.name}\n", encoding="ascii")
    return {"status": "pass", "report": str(output), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(capture(args.prefix, args.output), sort_keys=True))


if __name__ == "__main__":
    main()
