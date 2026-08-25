#!/usr/bin/env python3
"""Run the frozen endpoint-energy model with three replicas in CPU partitions."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


def _load_contract_module() -> Any:
    path = Path(__file__).with_name("run_gmx_mmpbsa_canary.py")
    spec = importlib.util.spec_from_file_location("run_gmx_mmpbsa_canary", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load endpoint-energy contract: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONTRACT = _load_contract_module()
MODELS = CONTRACT.MODEL_NAMES
REPLICAS = ("rep01", "rep02", "rep03")


def _load_execution_defaults_validator() -> Any:
    path = Path(__file__).with_name("validate_secondary_endpoint_execution_defaults.py")
    spec = importlib.util.spec_from_file_location("validate_endpoint_defaults", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load execution-defaults validator: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_summary_contract() -> Any:
    path = Path(__file__).with_name("summarize_secondary_endpoint_energy.py")
    spec = importlib.util.spec_from_file_location("endpoint_summary_contract", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load endpoint summary contract: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_formal_command(
    prefix: Path,
    cpu_set: str,
    mpi_ranks: int,
    input_file: Path,
    structure: Path,
    trajectory: Path,
    index: Path,
    topology: Path,
    reference: Path,
    final_text: Path,
    final_csv: Path,
    decomp_text: Path,
    decomp_csv: Path,
) -> list[str]:
    if mpi_ranks <= 0:
        raise ValueError("mpi_ranks must be positive")
    return [
        "taskset", "-c", cpu_set,
        str(prefix / "bin" / "mpirun"), "--bind-to", "none", "-np", str(mpi_ranks),
        str(prefix / "bin" / "gmx_MMPBSA"), "MPI", "-O",
        "-i", str(input_file),
        "-cs", str(structure),
        "-ct", str(trajectory),
        "-ci", str(index),
        "-cg", "0", "1",
        "-cp", str(topology),
        "-cr", str(reference),
        "-o", str(final_text),
        "-eo", str(final_csv),
        "-do", str(decomp_text),
        "-deo", str(decomp_csv),
        "-prefix", "_GMXMMPBSA_",
        "-nogui",
    ]


def _expand_cpu_set(value: str) -> set[int]:
    cpus: set[int] = set()
    for item in value.split(","):
        item = item.strip()
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", item)
        if match is None:
            raise ValueError(f"invalid CPU set: {value}")
        first = int(match.group(1))
        last = int(match.group(2) or first)
        if last < first:
            raise ValueError(f"invalid CPU set: {value}")
        cpus.update(range(first, last + 1))
    return cpus


def validate_model_resource_plan(model: dict[str, Any]) -> dict[str, Any]:
    if int(model.get("concurrent_jobs", -1)) != 3 or model.get("replica_batch") != list(REPLICAS):
        raise ValueError("resource plan must run exactly three replicas concurrently")
    ranks = int(model.get("mpi_ranks_per_job", -1))
    if ranks <= 0 or int(model.get("total_mpi_ranks", -1)) != 3 * ranks:
        raise ValueError("resource plan MPI rank counts are inconsistent")
    cpu_sets = model.get("cpu_sets")
    if not isinstance(cpu_sets, dict) or set(cpu_sets) != set(REPLICAS):
        raise ValueError("resource plan CPU sets are incomplete")
    expanded = {replica: _expand_cpu_set(str(cpu_sets[replica])) for replica in REPLICAS}
    if any(len(values) != ranks for values in expanded.values()):
        raise ValueError("each CPU set must contain one CPU per MPI rank")
    if any(expanded[left] & expanded[right] for index, left in enumerate(REPLICAS) for right in REPLICAS[index + 1:]):
        raise ValueError("resource plan CPU sets overlap")
    return {"status": "pass", "mpi_ranks_per_job": ranks}


def validate_manifest_records(root: Path, records: list[dict[str, Any]]) -> dict[str, Any]:
    if not isinstance(records, list) or not records:
        raise ValueError("migration manifest has no hash records")
    seen = set()
    for record in records:
        relative = str(record.get("path", ""))
        if not relative or relative in seen or relative.startswith(("/", "\\")) or ".." in Path(relative).parts:
            raise ValueError("migration manifest path/hash record is invalid")
        seen.add(relative)
        path = root / relative
        expected = record.get("sha256")
        if not path.is_file() or not isinstance(expected, str) or file_sha256(path) != expected:
            raise ValueError(f"migration package hash mismatch: {relative}")
        if int(record.get("bytes", -1)) != path.stat().st_size:
            raise ValueError(f"migration package hash record size mismatch: {relative}")
    return {"status": "pass", "file_count": len(records)}


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _prepared_file(package_root: Path, replica: str, manifest: dict[str, Any], name: str) -> Path:
    path = package_root / "prepared" / replica / name
    expected = manifest.get("outputs", {}).get(name, {}).get("sha256")
    if not path.is_file() or file_sha256(path) != expected:
        raise ValueError(f"{replica}: prepared input hash mismatch: {name}")
    return path.resolve()


def _make_environment(prefix: Path) -> dict[str, str]:
    env = dict(os.environ)
    env.update({
        "PATH": str(prefix / "bin") + os.pathsep + env.get("PATH", ""),
        "CONDA_PREFIX": str(prefix),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OMPI_ALLOW_RUN_AS_ROOT": "1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
    })
    return env


def execute(args: argparse.Namespace) -> dict[str, Any]:
    package_root = args.package_root.resolve()
    output_dir = args.output_dir.resolve()
    prefix = args.prefix.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    migration = _load_json(package_root / "MIGRATION_MANIFEST.json")
    if migration.get("status") != "sealed_after_all_md_gates_before_endpoint_energy_results":
        raise ValueError("migration package is not sealed at the required gate")
    validate_manifest_records(package_root, migration.get("files", []))
    defaults_path = package_root / "metadata" / "secondary_endpoint_energy_execution_defaults_v1.json"
    supplement_path = package_root / "metadata" / "SECONDARY_ENDPOINT_ENERGY_EXECUTION_SUPPLEMENT_20260820.md"
    plip_path = package_root / "metadata" / "8KCT_O6U_native_contacts.interactions.normalized.json"
    _load_execution_defaults_validator().validate(defaults_path, supplement_path, plip_path)
    toolchain = _load_json(args.toolchain_report)
    CONTRACT.validate_toolchain_record(toolchain.get("toolchain", {}))
    canary = _load_json(args.canary_report)
    CONTRACT.validate_canary_report(canary)
    resource = _load_json(args.resource_plan)
    if resource.get("status") != "frozen_before_formal_endpoint_energy_results" or resource.get("gpu_required") is not False:
        raise ValueError("formal resource plan is not frozen or incorrectly requires a GPU")
    geometry = _load_json(args.geometry)
    mthick = float(geometry.get("mthick_angstrom", 0.0))
    model_inputs = CONTRACT.render_model_inputs(mthick, 300)
    topology_map = migration.get("topology_relative_paths", {})
    if set(topology_map) != set(REPLICAS):
        raise ValueError("migration topology map is incomplete")
    prep_manifests = {
        replica: _load_json(package_root / "prepared" / replica / "PREPARATION_MANIFEST.json")
        for replica in REPLICAS
    }
    for replica, manifest in prep_manifests.items():
        if manifest.get("status") != "pass" or manifest.get("replica") != replica or int(manifest.get("frame_count", -1)) != 300:
            raise ValueError(f"{replica}: preparation manifest is invalid")

    output_dir.mkdir(parents=True)
    jobs: dict[str, dict[str, dict[str, Any]]] = {}
    for model in MODELS:
        model_resource = resource.get("models", {}).get(model, {})
        validation = validate_model_resource_plan(model_resource)
        ranks = validation["mpi_ranks_per_job"]
        jobs[model] = {}
        for replica in REPLICAS:
            work = output_dir / model / replica
            work.mkdir(parents=True)
            topology_stage = work / "topology"
            shutil.copytree(package_root / "topology" / replica, topology_stage)
            topology = topology_stage / Path(topology_map[replica])
            if not topology.is_file():
                raise ValueError(f"{replica}: staged topology is missing")
            input_file = work / f"{model}.in"
            input_file.write_text(model_inputs[model], encoding="ascii")
            manifest = prep_manifests[replica]
            command = build_formal_command(
                prefix=prefix,
                cpu_set=model_resource["cpu_sets"][replica],
                mpi_ranks=ranks,
                input_file=input_file.resolve(),
                structure=_prepared_file(package_root, replica, manifest, f"{replica}_endpoint_structure.gro"),
                trajectory=_prepared_file(package_root, replica, manifest, f"{replica}_endpoint_300frames_midplane0.xtc"),
                index=_prepared_file(package_root, replica, manifest, "endpoint_groups.ndx"),
                topology=topology.resolve(),
                reference=_prepared_file(package_root, replica, manifest, f"{replica}_endpoint_complex_reference.pdb"),
                final_text=(work / "FINAL_RESULTS_MMPBSA.dat").resolve(),
                final_csv=(work / "FINAL_RESULTS_MMPBSA.csv").resolve(),
                decomp_text=(work / "FINAL_DECOMP_MMPBSA.dat").resolve(),
                decomp_csv=(work / "FINAL_DECOMP_MMPBSA.csv").resolve(),
            )
            jobs[model][replica] = {
                "work": work,
                "command": command,
                "input_sha256": file_sha256(input_file),
                "cpu_set": model_resource["cpu_sets"][replica],
                "mpi_ranks": ranks,
            }
    release = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_formal_run_release",
        "status": "released_before_formal_results",
        "model_order": list(MODELS),
        "replicas_concurrent_within_model": list(REPLICAS),
        "models_sequential": True,
        "migration_manifest_sha256": file_sha256(package_root / "MIGRATION_MANIFEST.json"),
        "toolchain_report_sha256": file_sha256(args.toolchain_report),
        "canary_report_sha256": file_sha256(args.canary_report),
        "resource_plan_sha256": file_sha256(args.resource_plan),
        "geometry_sha256": file_sha256(args.geometry),
        "execution_defaults_sha256": file_sha256(defaults_path),
        "execution_supplement_sha256": file_sha256(supplement_path),
        "native_plip_contact_set_sha256": file_sha256(plip_path),
        "jobs": {
            model: {
                replica: {key: value for key, value in record.items() if key != "work"}
                for replica, record in replicas.items()
            }
            for model, replicas in jobs.items()
        },
    }
    release_path = output_dir / "FORMAL_RUN_RELEASE.json"
    release_path.write_text(json.dumps(release, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    release_digest = file_sha256(release_path)
    release_path.with_suffix(".json.sha256").write_text(f"{release_digest}  {release_path.name}\n", encoding="ascii")

    if not Path("/usr/bin/time").is_file() or shutil.which("taskset") is None:
        raise ValueError("formal execution requires /usr/bin/time and taskset")
    env = _make_environment(prefix)
    completed_jobs: dict[str, dict[str, Any]] = {}
    for model in MODELS:
        running = {}
        handles = {}
        for replica in REPLICAS:
            record = jobs[model][replica]
            work = record["work"]
            log_path = work / "run.log"
            time_path = work / "time_verbose.txt"
            handle = log_path.open("w", encoding="utf-8")
            command = ["/usr/bin/time", "-v", "-o", str(time_path), *record["command"]]
            running[replica] = subprocess.Popen(
                command,
                cwd=work,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
            handles[replica] = handle
        returncodes = {}
        for replica in REPLICAS:
            returncodes[replica] = running[replica].wait()
            handles[replica].close()
        if any(code != 0 for code in returncodes.values()):
            failure_path = output_dir / f"{model}_BATCH_FAILURE.json"
            failure_path.write_text(
                json.dumps({"model": model, "returncodes": returncodes, "next_models_started": False}, indent=2) + "\n",
                encoding="utf-8",
            )
            raise RuntimeError(f"{model}: one or more frozen replica jobs failed; no later model was started")
        completed_jobs[model] = {}
        for replica in REPLICAS:
            work = jobs[model][replica]["work"]
            final_text = work / "FINAL_RESULTS_MMPBSA.dat"
            final_csv = work / "FINAL_RESULTS_MMPBSA.csv"
            decomp_text = work / "FINAL_DECOMP_MMPBSA.dat"
            decomp_csv = work / "FINAL_DECOMP_MMPBSA.csv"
            values = CONTRACT.parse_result_components(final_text.read_text(encoding="utf-8", errors="replace"))
            CONTRACT.parse_result_components(decomp_text.read_text(encoding="utf-8", errors="replace"))
            summary_contract = _load_summary_contract()
            decomp_rows = summary_contract.parse_decomposition_csv(
                decomp_csv.read_text(encoding="utf-8", errors="replace")
            )
            decomp_validation = summary_contract.validate_fixed_decomposition_rows(
                decomp_rows, frame_count=300
            )
            topologies = CONTRACT.classify_generated_topologies(
                [path for path in work.glob("_GMXMMPBSA_*.*top") if path.is_file()]
            )
            completed_jobs[model][replica] = {
                "returncode": 0,
                "finite_component_count": len(values),
                "peak_rss_bytes_per_rank": CONTRACT.parse_peak_rss_bytes(
                    (work / "time_verbose.txt").read_text(encoding="utf-8"),
                    jobs[model][replica]["mpi_ranks"],
                ),
                "outputs": {
                    final_text.name: {"sha256": file_sha256(final_text), "bytes": final_text.stat().st_size},
                    final_csv.name: {"sha256": file_sha256(final_csv), "bytes": final_csv.stat().st_size},
                    decomp_text.name: {"sha256": file_sha256(decomp_text), "bytes": decomp_text.stat().st_size},
                    decomp_csv.name: {"sha256": file_sha256(decomp_csv), "bytes": decomp_csv.stat().st_size},
                },
                "decomposition_status": "pass",
                "decomposition_residues_in_frozen_order": decomp_validation[
                    "residues_in_frozen_order"
                ],
                "generated_topology_sha256": {label: file_sha256(path) for label, path in topologies.items()},
            }
    completion = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_formal_run_completion",
        "status": "pass",
        "release_sha256": release_digest,
        "jobs": completed_jobs,
        "raw_outputs_immutable": True,
    }
    completion_path = output_dir / "FORMAL_RUN_COMPLETION.json"
    completion_path.write_text(json.dumps(completion, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(completion_path)
    completion_path.with_suffix(".json.sha256").write_text(f"{digest}  {completion_path.name}\n", encoding="ascii")
    return {"status": "pass", "completion": str(completion_path), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--toolchain-report", type=Path, required=True)
    parser.add_argument("--canary-report", type=Path, required=True)
    parser.add_argument("--resource-plan", type=Path, required=True)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(execute(args), sort_keys=True))


if __name__ == "__main__":
    main()
