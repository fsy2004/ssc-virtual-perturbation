#!/usr/bin/env python3
"""Independently reconstruct the frozen O6U parameterization toolchain record."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_ARCHIVE_SHA256 = "d9508f3a1590ba9fbfb1d048e832d6726ef6131fd40a70572f5662c5bdc2cbdb"
EXPECTED_DISTRIBUTIONS = {
    "ffparam": "1.2.0",
    "openmm": "8.1.2",
    "rdkit": "2026.3.5",
    "numpy": "1.26.4",
    "scipy": "1.13.1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reconstruct_tree(root: Path) -> tuple[list[dict[str, object]], str]:
    entries: list[dict[str, object]] = []
    aggregate = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        digest = sha256(path)
        size = path.stat().st_size
        entry = {"path": relative, "size_bytes": size, "sha256": digest}
        entries.append(entry)
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(str(size).encode("ascii"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
    return entries, aggregate.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    record_path = args.record.resolve()
    report_path = args.report.resolve()
    if not record_path.is_file() or record_path.stat().st_size == 0:
        fail(f"Missing or empty toolchain record: {record_path}")
    if report_path.exists():
        fail(f"Refusing to overwrite validation report: {report_path}")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    if record.get("schema_version") != "1.0" or record.get("status") != "pass":
        fail("Toolchain record schema/status is invalid")
    if record.get("security") != "No environment variables, credentials, tokens, or command histories are recorded.":
        fail("Toolchain security declaration differs")

    archive_record = record.get("official_ffparam_archive", {})
    archive_path = Path(str(archive_record.get("resolved_path", ""))).resolve()
    if not archive_path.is_file():
        fail(f"FFParam archive is absent: {archive_path}")
    archive_hash = sha256(archive_path)
    if (
        archive_hash != EXPECTED_ARCHIVE_SHA256
        or archive_record.get("sha256") != archive_hash
        or archive_record.get("size_bytes") != archive_path.stat().st_size
    ):
        fail("FFParam 1.2.0 archive identity differs")

    source_record = record.get("ffparam_source_tree", {})
    source_root = Path(str(source_record.get("root", ""))).resolve()
    if not source_root.is_dir():
        fail(f"FFParam source root is absent: {source_root}")
    entries, tree_hash = reconstruct_tree(source_root)
    if (
        source_record.get("file_count") != len(entries)
        or source_record.get("files") != entries
        or source_record.get("tree_sha256") != tree_hash
    ):
        fail("FFParam source-tree reconstruction differs")

    python_record = record.get("python", {})
    distributions = python_record.get("distributions", {})
    observed_distributions: dict[str, str | None] = {}
    for package in ("ffparam", "psi4", "openmm", "rdkit", "numpy", "scipy"):
        try:
            observed_distributions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            observed_distributions[package] = None
    if distributions != observed_distributions:
        fail("Installed Python distribution versions differ from the frozen record")
    for package, expected in EXPECTED_DISTRIBUTIONS.items():
        if observed_distributions.get(package) != expected:
            fail(f"Unexpected {package} distribution: {observed_distributions.get(package)!r}")

    module_results: dict[str, dict[str, object]] = {}
    module_records = python_record.get("modules", {})
    for package in ("ffparam", "psi4", "openmm", "rdkit", "numpy", "scipy"):
        module = importlib.import_module(package)
        module_path = Path(str(getattr(module, "__file__", ""))).resolve()
        frozen = module_records.get(package, {})
        frozen_file = frozen.get("module_file", {})
        observed_hash = sha256(module_path)
        if (
            not module_path.is_file()
            or frozen.get("import_status") != "pass"
            or frozen_file.get("resolved_path") != str(module_path)
            or frozen_file.get("sha256") != observed_hash
            or frozen.get("version_attribute") != getattr(module, "__version__", None)
        ):
            fail(f"Imported module identity differs: {package}")
        module_results[package] = {
            "module_path": str(module_path),
            "sha256": observed_hash,
            "version_attribute": getattr(module, "__version__", None),
        }
    if module_results["psi4"]["version_attribute"] != "1.9.1":
        fail(f"Unexpected Psi4 module version: {module_results['psi4']['version_attribute']!r}")

    executable_results: dict[str, dict[str, object]] = {}
    for role in ("crest", "xtb", "gromacs"):
        frozen = record.get("executables", {}).get(role, {})
        executable = Path(str(frozen.get("resolved_path", ""))).resolve()
        if not executable.is_file() or sha256(executable) != frozen.get("sha256"):
            fail(f"Executable identity differs: {role}")
        command = frozen.get("version_command")
        if not isinstance(command, list) or not command or command[0] != str(executable):
            fail(f"Frozen version command is invalid: {role}")
        completed = subprocess.run(command, text=True, capture_output=True, timeout=60, check=False)
        if (
            completed.returncode != frozen.get("version_returncode")
            or completed.stdout.strip() != frozen.get("version_stdout")
            or completed.stderr.strip() != frozen.get("version_stderr")
        ):
            fail(f"Executable version output differs: {role}")
        executable_results[role] = {
            "path": str(executable),
            "sha256": sha256(executable),
            "version_returncode": completed.returncode,
        }

    report = {
        "schema_version": "1.0",
        "report_type": "independent_o6u_parameterization_toolchain_validation",
        "status": "pass_independently_reconstructed",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_record": {"path": str(record_path), "sha256": sha256(record_path)},
        "ffparam_archive": {"path": str(archive_path), "sha256": archive_hash},
        "ffparam_source_tree": {"root": str(source_root), "file_count": len(entries), "tree_sha256": tree_hash},
        "python_distributions": observed_distributions,
        "python_modules": module_results,
        "executables": executable_results,
        "release_boundary": "Toolchain identity only; QM targets, ligand parameters, CHARMM-GUI, and MD remain blocked.",
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report_sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
