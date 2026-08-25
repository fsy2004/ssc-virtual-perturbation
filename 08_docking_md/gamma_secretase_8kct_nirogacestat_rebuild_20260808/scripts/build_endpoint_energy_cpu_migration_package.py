#!/usr/bin/env python3
"""Build a hash-bound, immutable CPU migration package after all MD gates pass."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
from pathlib import Path, PurePosixPath
from typing import Any


REPLICAS = ("rep01", "rep02", "rep03")
CODE_FILES = (
    "validate_secondary_endpoint_energy_plan.py",
    "validate_secondary_endpoint_execution_defaults.py",
    "seal_secondary_endpoint_all_three_gate.py",
    "prepare_secondary_endpoint_energy_inputs.py",
    "freeze_endpoint_energy_membrane_geometry.py",
    "install_gmx_mmpbsa_1_6_5_cpu.sh",
    "capture_gmx_mmpbsa_toolchain.py",
    "run_gmx_mmpbsa_canary.py",
    "collect_endpoint_cpu_inventory.py",
    "plan_secondary_endpoint_resources.py",
    "run_secondary_endpoint_energy_cpu.py",
    "summarize_secondary_endpoint_energy.py",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_replica_mapping(values: list[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"invalid replica mapping: {value}")
        replica, raw_path = value.split("=", 1)
        if replica not in REPLICAS:
            raise ValueError(f"unknown replica: {replica}")
        if replica in parsed:
            raise ValueError(f"duplicate replica mapping: {replica}")
        parsed[replica] = Path(raw_path)
    missing = [replica for replica in REPLICAS if replica not in parsed]
    if missing:
        raise ValueError("missing replica mappings: " + ", ".join(missing))
    return parsed


def safe_relative_path(value: str) -> Path:
    normalized = value.replace("\\", "/")
    posix = PurePosixPath(normalized)
    if value.startswith(("/", "\\")) or Path(value).drive or not posix.parts or any(part in ("", ".", "..") for part in posix.parts):
        raise ValueError("topology path must be a safe relative path")
    return Path(*posix.parts)


def manifest_file_records(root: Path) -> list[dict[str, Any]]:
    records = []
    for path in sorted((path for path in root.rglob("*") if path.is_file()), key=lambda item: item.relative_to(root).as_posix()):
        relative = path.relative_to(root).as_posix()
        records.append({"path": relative, "bytes": path.stat().st_size, "sha256": file_sha256(path)})
    return records


def validate_completion_gate(gate: dict[str, Any]) -> dict[str, Any]:
    if gate.get("status") != "pass":
        raise ValueError("all-three completion gate status must pass")
    if gate.get("all_three_500ns_complete") is not True:
        raise ValueError("all three 500 ns realizations are not complete")
    if gate.get("all_required_gates_passed") is not True:
        raise ValueError("required integrity/PBC/membrane/energy gates have not all passed")
    if gate.get("production_runners_active") not in (False, 0, []):
        raise ValueError("production runner is still active")
    return {"status": "pass"}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _verify_preparation_directory(replica: str, source: Path) -> None:
    manifest_path = source / "PREPARATION_MANIFEST.json"
    manifest = _load_json(manifest_path)
    if manifest.get("status") != "pass" or manifest.get("replica") != replica or int(manifest.get("frame_count", -1)) != 300:
        raise ValueError(f"{replica}: preparation manifest is not passing")
    for name, record in manifest.get("outputs", {}).items():
        path = source / name
        if not path.is_file() or file_sha256(path) != record.get("sha256"):
            raise ValueError(f"{replica}: prepared output is missing or drifted: {name}")


def build_package(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = args.output_dir.resolve()
    archive = args.archive.resolve()
    if output_dir.exists() or archive.exists() or archive.with_suffix(archive.suffix + ".sha256").exists():
        raise FileExistsError("refusing to overwrite migration package or archive")
    prep_dirs = parse_replica_mapping(args.prep)
    topology_dirs = parse_replica_mapping(args.topology_bundle)
    # Parse topology-relative values by their explicit replica key, independent of argument order.
    topology_paths = {replica: safe_relative_path(str(path)) for replica, path in parse_replica_mapping(args.topology_relative).items()}
    gate = _load_json(args.completion_gate)
    validate_completion_gate(gate)
    geometry = _load_json(args.geometry)
    if geometry.get("status") != "pass" or float(geometry.get("mctrdz_angstrom", 1.0)) != 0.0:
        raise ValueError("common membrane geometry is not passing")

    output_dir.mkdir(parents=True)
    for replica in REPLICAS:
        prep_source = prep_dirs[replica].resolve()
        _verify_preparation_directory(replica, prep_source)
        shutil.copytree(prep_source, output_dir / "prepared" / replica)
        topology_source = topology_dirs[replica].resolve()
        if not (topology_source / topology_paths[replica]).is_file():
            raise ValueError(f"{replica}: topology bundle lacks {topology_paths[replica]}")
        shutil.copytree(topology_source, output_dir / "topology" / replica)

    metadata = output_dir / "metadata"
    metadata.mkdir()
    for source in (
        args.completion_gate,
        args.geometry,
        args.plan,
        args.amendment,
        args.execution_defaults,
        args.execution_supplement,
        args.native_contact_set,
    ):
        shutil.copy2(source, metadata / source.name)
    scripts_dir = Path(__file__).resolve().parent
    code_dir = output_dir / "code"
    code_dir.mkdir()
    for name in CODE_FILES:
        source = scripts_dir / name
        if not source.is_file():
            raise ValueError(f"migration code file is missing: {name}")
        shutil.copy2(source, code_dir / name)

    topology_map = {replica: topology_paths[replica].as_posix() for replica in REPLICAS}
    (metadata / "topology_relative_paths.json").write_text(
        json.dumps(topology_map, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme = (
        "O6U secondary endpoint-energy CPU migration package\n\n"
        "This package is immutable. Verify MIGRATION_MANIFEST.json before use.\n"
        "Install the isolated CPU environment, capture its toolchain record, run the rep01 three-frame canary, "
        "freeze the resource plan from canary RSS, then run the formal batch.\n"
        "Do not run concurrently with MD production and do not overwrite raw or prepared inputs.\n"
    )
    (output_dir / "README_CPU_MIGRATION.txt").write_text(readme, encoding="utf-8")
    records = manifest_file_records(output_dir)
    manifest = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_cpu_migration_package",
        "status": "sealed_after_all_md_gates_before_endpoint_energy_results",
        "replicas": list(REPLICAS),
        "topology_relative_paths": topology_map,
        "completion_gate_sha256": file_sha256(args.completion_gate),
        "geometry_sha256": file_sha256(args.geometry),
        "plan_sha256": file_sha256(args.plan),
        "amendment_sha256": file_sha256(args.amendment),
        "execution_defaults_sha256": file_sha256(args.execution_defaults),
        "execution_supplement_sha256": file_sha256(args.execution_supplement),
        "native_contact_set_sha256": file_sha256(args.native_contact_set),
        "files": records,
    }
    manifest_path = output_dir / "MIGRATION_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_digest = file_sha256(manifest_path)
    (output_dir / "MIGRATION_MANIFEST.json.sha256").write_text(
        f"{manifest_digest}  MIGRATION_MANIFEST.json\n", encoding="ascii"
    )
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, mode="w", format=tarfile.PAX_FORMAT) as handle:
        handle.add(output_dir, arcname=output_dir.name, recursive=True)
    archive_digest = file_sha256(archive)
    archive.with_suffix(archive.suffix + ".sha256").write_text(
        f"{archive_digest}  {archive.name}\n", encoding="ascii"
    )
    return {
        "status": "pass",
        "package": str(output_dir),
        "manifest_sha256": manifest_digest,
        "archive": str(archive),
        "archive_sha256": archive_digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prep", action="append", required=True, help="repXX=/absolute/preparation/dir")
    parser.add_argument("--topology-bundle", action="append", required=True, help="repXX=/absolute/topology/bundle")
    parser.add_argument("--topology-relative", action="append", required=True, help="repXX=relative/path/to/topol.top")
    parser.add_argument("--completion-gate", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--amendment", type=Path, required=True)
    parser.add_argument("--execution-defaults", type=Path, required=True)
    parser.add_argument("--execution-supplement", type=Path, required=True)
    parser.add_argument("--native-contact-set", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_package(args), sort_keys=True))


if __name__ == "__main__":
    main()
