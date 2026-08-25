#!/usr/bin/env python3
"""Freeze and validate the three-frame gmx_MMPBSA technical canary."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


GMX_MMPBSA_COMMIT = "64e994c71aaff315f3c82dd0852919aecb1ab62e"
MODEL_NAMES = ("PB_membrane_indi4",)
REQUIRED_EXECUTABLES = ("gmx_MMPBSA", "mpirun", "sander", "cpptraj", "tleap", "gmx")
FROZEN_DECOMP_PRINT_RES = "B/261,268,272,282,287,380-381,431-432,502"
FROZEN_DECOMP_RESIDUES = (
    "R:B:VAL:261",
    "R:B:LEU:268",
    "R:B:VAL:272",
    "R:B:LEU:282",
    "R:B:ILE:287",
    "R:B:LYS:380",
    "R:B:LEU:381",
    "R:B:ALA:431",
    "R:B:LEU:432",
    "L:B:O6U:502",
)
EXACT_TOOLCHAIN = {
    "gmx_mmpbsa": "1.6.5",
    "gmx_mmpbsa_git_commit": GMX_MMPBSA_COMMIT,
    "python": "3.11.8",
    "ambertools": "23.3",
    "gromacs": "2023.4",
    "openmpi": "4.1.6",
    "mpi4py": "4.0.1",
    "numpy": "1.26.4",
    "pandas": "1.5.3",
    "matplotlib": "3.7.3",
    "seaborn": "0.11.2",
    "scipy": "1.14.1",
    "tqdm": "4.67.1",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prespecified_canary_frames() -> list[dict[str, float | int]]:
    return [
        {"output_index_zero_based": 0, "target_time_ns": 200.5},
        {"output_index_zero_based": 150, "target_time_ns": 350.5},
        {"output_index_zero_based": 299, "target_time_ns": 499.5},
    ]


def build_gmx_mmpbsa_command(
    prefix: Path,
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
    return [
        str(prefix / "bin" / "mpirun"),
        "-np",
        "3",
        str(prefix / "bin" / "gmx_MMPBSA"),
        "MPI",
        "-O",
        "-i",
        str(input_file),
        "-cs",
        str(structure),
        "-ct",
        str(trajectory),
        "-ci",
        str(index),
        "-cg",
        "0",
        "1",
        "-cp",
        str(topology),
        "-cr",
        str(reference),
        "-o",
        str(final_text),
        "-eo",
        str(final_csv),
        "-do",
        str(decomp_text),
        "-deo",
        str(decomp_csv),
        "-prefix",
        "_GMXMMPBSA_",
        "-nogui",
    ]


def parse_result_components(text: str) -> list[float]:
    if re.search(r"(?i)(?:^|[\s=,])(?:nan|[+-]?inf(?:inity)?)(?=$|[\s,])", text):
        raise ValueError("result contains non-finite component values")
    values = [
        float(value)
        for value in re.findall(r"(?<![A-Za-z0-9_.])[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?", text)
    ]
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("result does not contain finite component values")
    return values


def classify_generated_topologies(paths: list[Path]) -> dict[str, Path]:
    markers = {"complex": "COM", "receptor": "REC", "ligand": "LIG"}
    result: dict[str, Path] = {}
    for label, marker in markers.items():
        candidates = [path for path in paths if re.search(rf"(?:^|_){marker}(?:_|\.|$)", path.name, re.IGNORECASE)]
        if len(candidates) != 1:
            raise ValueError(f"generated {label} topology count must be exactly one")
        result[label] = candidates[0]
    return result


def parse_peak_rss_bytes(text: str, mpi_ranks: int) -> int:
    if mpi_ranks <= 0:
        raise ValueError("mpi_ranks must be positive")
    match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if match is None or int(match.group(1)) <= 0:
        raise ValueError("GNU time maximum resident set size is missing or invalid")
    return int(match.group(1)) * 1024


def _general_block(frame_count: int, model_name: str) -> str:
    if frame_count <= 0:
        raise ValueError("frame_count must be positive")
    return (
        "&general\n"
        f"  sys_name=\"O6U_{model_name}\",\n"
        "  startframe=1,\n"
        f"  endframe={frame_count},\n"
        "  interval=1,\n"
        "  keep_files=2,\n"
        "  PBRadii=7,\n"
        "/\n"
    )


def _pb_block(indi: float, mthick_angstrom: float) -> str:
    if not math.isfinite(mthick_angstrom) or mthick_angstrom <= 0:
        raise ValueError("mthick_angstrom must be finite and positive")
    return (
        "&pb\n"
        "  memopt=1,\n"
        "  emem=7.0,\n"
        f"  indi={indi:.1f},\n"
        "  exdi=80.0,\n"
        "  istrng=0.150,\n"
        "  poretype=1,\n"
        "  mctrdz=0.0,\n"
        f"  mthick={mthick_angstrom:.6f},\n"
        "  radiopt=0,\n"
        "  fillratio=1.25,\n"
        "  inp=2,\n"
        "  sasopt=0,\n"
        "  solvopt=2,\n"
        "  ipb=1,\n"
        "  bcopt=10,\n"
        "  nfocus=1,\n"
        "  linit=1000,\n"
        "  eneopt=1,\n"
        "  cutfd=7.0,\n"
        "  cutnb=99.0,\n"
        "  maxarcdot=15000,\n"
        "  npbverb=1,\n"
        "/\n"
    )


def _gb_block(igb: int) -> str:
    if igb not in (5, 8):
        raise ValueError("frozen GB models require igb 5 or 8")
    return (
        "&gb\n"
        f"  igb={igb},\n"
        "  intdiel=4.0,\n"
        "  extdiel=80.0,\n"
        "  saltcon=0.150,\n"
        "/\n"
    )


def _decomp_block() -> str:
    return (
        "&decomp\n"
        "  idecomp=2,\n"
        "  dec_verbose=0,\n"
        "  csv_format=1,\n"
        f"  print_res=\"{FROZEN_DECOMP_PRINT_RES}\",\n"
        "/\n"
    )


def render_model_inputs(mthick_angstrom: float, frame_count: int) -> dict[str, str]:
    return {
        "PB_membrane_indi4": _general_block(frame_count, "PB_membrane_indi4") + _pb_block(4.0, mthick_angstrom) + _decomp_block(),
    }


def _version_tuple(value: str) -> tuple[int, ...]:
    try:
        return tuple(int(part) for part in value.split("."))
    except ValueError as exc:
        raise ValueError(f"invalid package version: {value}") from exc


def validate_toolchain_record(record: dict[str, Any]) -> dict[str, Any]:
    for field, expected in EXACT_TOOLCHAIN.items():
        if str(record.get(field)) != expected:
            raise ValueError(f"{field} must be exactly {expected}")
    if str(record.get("parmed")) != "4.3.0":
        raise ValueError("parmed must be exactly 4.3.0")
    if record.get("gpu_required") is not False:
        raise ValueError("endpoint toolchain must record gpu_required=false")
    executables = record.get("executables")
    if not isinstance(executables, dict) or set(executables) != set(REQUIRED_EXECUTABLES):
        raise ValueError("required executable set is incomplete or has drifted")
    for name in REQUIRED_EXECUTABLES:
        digest = executables[name].get("sha256") if isinstance(executables[name], dict) else None
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError(f"{name}: missing executable SHA256")
    return {"status": "pass", "toolchain": dict(EXACT_TOOLCHAIN), "gpu_required": False}


def example_passing_canary_report() -> dict[str, Any]:
    models = {}
    for offset, name in enumerate(MODEL_NAMES):
        models[name] = {
            "frame_count": 3,
            "finite_components": True,
            "component_values": [-10.0 - offset, 2.0, 4.0],
            "generated_topology_sha256": {
                "complex": "a" * 64,
                "receptor": "b" * 64,
                "ligand": "c" * 64,
            },
            "decomposition_status": "pass",
            "decomposition_output_sha256": "d" * 64,
            "decomposition_residues_in_frozen_order": list(FROZEN_DECOMP_RESIDUES),
        }
    return {
        "models": models,
        "complex_atom_count": 5000,
        "receptor_atom_count": 4924,
        "ligand_atom_count": 76,
        "ligand_total_charge": 0.0,
        "ligand_partial_charges_all_zero": False,
        "toolchain_status": "pass",
    }


def validate_canary_report(report: dict[str, Any]) -> dict[str, Any]:
    if report.get("toolchain_status") != "pass":
        raise ValueError("toolchain status must pass")
    receptor = int(report.get("receptor_atom_count", -1))
    ligand = int(report.get("ligand_atom_count", -1))
    complex_count = int(report.get("complex_atom_count", -1))
    if ligand != 76 or receptor <= 0 or complex_count != receptor + ligand:
        raise ValueError("receptor/ligand/complex atom counts are inconsistent")
    if abs(float(report.get("ligand_total_charge", math.inf))) > 1e-6:
        raise ValueError("ligand total charge is not neutral")
    if report.get("ligand_partial_charges_all_zero") is not False:
        raise ValueError("ligand partial charges must be nonzero and sum to zero")
    models = report.get("models")
    if not isinstance(models, dict) or set(models) != set(MODEL_NAMES):
        raise ValueError("canary must contain exactly the frozen executable endpoint-energy model set")
    for name in MODEL_NAMES:
        model = models[name]
        if int(model.get("frame_count", -1)) != 3:
            raise ValueError(f"{name}: canary must contain exactly three frames")
        values = model.get("component_values")
        if model.get("finite_components") is not True or not isinstance(values, list) or not values:
            raise ValueError(f"{name}: finite component validation is missing")
        if any(not math.isfinite(float(value)) for value in values):
            raise ValueError(f"{name}: component values must be finite")
        hashes = model.get("generated_topology_sha256")
        if not isinstance(hashes, dict) or set(hashes) != {"complex", "receptor", "ligand"}:
            raise ValueError(f"{name}: generated topology hashes are incomplete")
        if any(not isinstance(value, str) or len(value) != 64 for value in hashes.values()):
            raise ValueError(f"{name}: generated topology SHA256 is invalid")
        if model.get("decomposition_status") != "pass":
            raise ValueError(f"{name}: decomposition canary did not pass")
        decomp_hash = model.get("decomposition_output_sha256")
        if not isinstance(decomp_hash, str) or len(decomp_hash) != 64:
            raise ValueError(f"{name}: decomposition output SHA256 is invalid")
        if model.get("decomposition_residues_in_frozen_order") != list(FROZEN_DECOMP_RESIDUES):
            raise ValueError(f"{name}: decomposition residue mapping drifted")
    return {"status": "pass", "models": list(MODEL_NAMES), "frame_count_per_model": 3}


def write_model_inputs(output_dir: Path, mthick_angstrom: float, frame_count: int) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True)
    outputs = {}
    for name, text in render_model_inputs(mthick_angstrom, frame_count).items():
        path = output_dir / f"{name}.in"
        path.write_text(text, encoding="ascii")
        outputs[name] = {"path": str(path), "sha256": file_sha256(path)}
    manifest = {
        "schema_version": "1.0",
        "status": "frozen_before_canary_results",
        "frame_count": frame_count,
        "mthick_angstrom": mthick_angstrom,
        "models": outputs,
    }
    path = output_dir / "MODEL_INPUT_MANIFEST.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"status": "pass", "manifest": str(path), "sha256": file_sha256(path)}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _validate_execution_defaults(defaults: Path, supplement: Path, plip: Path) -> dict[str, Any]:
    validator_path = Path(__file__).with_name("validate_secondary_endpoint_execution_defaults.py")
    spec = importlib.util.spec_from_file_location("validate_endpoint_defaults", validator_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load execution-defaults validator: {validator_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.validate(defaults, supplement, plip)


def _load_decomposition_contract() -> Any:
    path = Path(__file__).with_name("summarize_secondary_endpoint_energy.py")
    spec = importlib.util.spec_from_file_location("endpoint_decomposition_contract", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load decomposition contract: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _verified_preparation_file(prep_dir: Path, manifest: dict[str, Any], name: str) -> Path:
    path = prep_dir / name
    expected = manifest.get("outputs", {}).get(name, {}).get("sha256")
    if not path.is_file() or not isinstance(expected, str) or file_sha256(path) != expected:
        raise ValueError(f"prepared input is missing or hash-drifted: {name}")
    return path.resolve()


def _topology_properties(path: Path) -> dict[str, Any]:
    import parmed

    topology = parmed.load_file(str(path))
    charges = [float(atom.charge) for atom in topology.atoms]
    return {
        "atom_count": len(charges),
        "total_charge": sum(charges),
        "partial_charges_all_zero": all(abs(value) <= 1e-12 for value in charges),
    }


def execute_canary(
    prefix: Path,
    prep_dir: Path,
    topology_bundle: Path,
    topology_relative_path: Path,
    geometry_path: Path,
    toolchain_report_path: Path,
    execution_defaults_path: Path,
    execution_supplement_path: Path,
    native_contact_set_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    prefix = prefix.resolve()
    prep_dir = prep_dir.resolve()
    topology_bundle = topology_bundle.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    prep_manifest_path = prep_dir / "PREPARATION_MANIFEST.json"
    prep_manifest = _load_json(prep_manifest_path)
    if (
        prep_manifest.get("status") != "pass"
        or prep_manifest.get("replica") != "rep01"
        or int(prep_manifest.get("frame_count", -1)) != 300
        or prep_manifest.get("canary_output_indices_zero_based") != [0, 150, 299]
        or prep_manifest.get("pbc_invariance", {}).get("status") != "pass"
    ):
        raise ValueError("rep01 preparation manifest does not satisfy the canary contract")
    names = {
        "structure": "rep01_endpoint_structure.gro",
        "reference": "rep01_endpoint_complex_reference.pdb",
        "trajectory": "rep01_endpoint_canary_3frames_midplane0.xtc",
        "index": "endpoint_groups.ndx",
    }
    prepared = {key: _verified_preparation_file(prep_dir, prep_manifest, name) for key, name in names.items()}
    geometry = _load_json(geometry_path)
    if geometry.get("status") != "pass" or float(geometry.get("mctrdz_angstrom", math.inf)) != 0.0:
        raise ValueError("frozen membrane geometry is not passing or centered at z=0")
    mthick = float(geometry.get("mthick_angstrom", math.nan))
    if not math.isfinite(mthick) or mthick <= 0:
        raise ValueError("frozen mthick is invalid")
    toolchain_report = _load_json(toolchain_report_path)
    if toolchain_report.get("status") != "pass":
        raise ValueError("toolchain report has not passed")
    validate_toolchain_record(toolchain_report.get("toolchain", {}))
    _validate_execution_defaults(
        execution_defaults_path, execution_supplement_path, native_contact_set_path
    )
    source_topology = topology_bundle / topology_relative_path
    if not source_topology.is_file():
        raise ValueError("topology bundle does not contain the requested topology")

    output_dir.mkdir(parents=True)
    model_inputs = render_model_inputs(mthick, 3)
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
    model_reports: dict[str, dict[str, Any]] = {}
    reference_counts: dict[str, Any] | None = None
    for model_name in MODEL_NAMES:
        model_dir = output_dir / model_name
        model_dir.mkdir()
        staged_topology = model_dir / "topology"
        shutil.copytree(topology_bundle, staged_topology)
        topology = staged_topology / topology_relative_path
        input_file = model_dir / f"{model_name}.in"
        input_file.write_text(model_inputs[model_name], encoding="ascii")
        final_text = model_dir / "FINAL_RESULTS_MMPBSA.dat"
        final_csv = model_dir / "FINAL_RESULTS_MMPBSA.csv"
        decomp_text = model_dir / "FINAL_DECOMP_MMPBSA.dat"
        decomp_csv = model_dir / "FINAL_DECOMP_MMPBSA.csv"
        command = build_gmx_mmpbsa_command(
            prefix,
            input_file.resolve(),
            prepared["structure"],
            prepared["trajectory"],
            prepared["index"],
            topology.resolve(),
            prepared["reference"],
            final_text.resolve(),
            final_csv.resolve(),
            decomp_text.resolve(),
            decomp_csv.resolve(),
        )
        time_record = model_dir / "time_verbose.txt"
        if not Path("/usr/bin/time").is_file():
            raise ValueError("/usr/bin/time is required to calibrate formal resource use")
        executed_command = ["/usr/bin/time", "-v", "-o", str(time_record), *command]
        log_path = model_dir / "run.log"
        with log_path.open("w", encoding="utf-8") as log_handle:
            completed = subprocess.run(
                executed_command,
                cwd=model_dir,
                env=env,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=4 * 60 * 60,
                check=False,
            )
        if completed.returncode != 0:
            raise RuntimeError(f"{model_name}: canary exited with code {completed.returncode}; inspect retained log")
        if not final_text.is_file() or not final_csv.is_file() or not decomp_text.is_file() or not decomp_csv.is_file():
            raise ValueError(f"{model_name}: expected final outputs are missing")
        component_values = parse_result_components(final_text.read_text(encoding="utf-8", errors="replace"))
        parse_result_components(decomp_text.read_text(encoding="utf-8", errors="replace"))
        decomp_contract = _load_decomposition_contract()
        decomp_rows = decomp_contract.parse_decomposition_csv(
            decomp_csv.read_text(encoding="utf-8", errors="replace")
        )
        decomp_validation = decomp_contract.validate_fixed_decomposition_rows(
            decomp_rows, frame_count=3
        )
        topology_paths = classify_generated_topologies(
            [path for path in model_dir.glob("_GMXMMPBSA_*.*top") if path.is_file()]
        )
        properties = {label: _topology_properties(path) for label, path in topology_paths.items()}
        if properties["ligand"]["atom_count"] != 76:
            raise ValueError(f"{model_name}: generated ligand does not contain 76 atoms")
        if properties["complex"]["atom_count"] != properties["receptor"]["atom_count"] + 76:
            raise ValueError(f"{model_name}: generated topology atom counts are inconsistent")
        if reference_counts is None:
            reference_counts = properties
        elif any(properties[label]["atom_count"] != reference_counts[label]["atom_count"] for label in properties):
            raise ValueError(f"{model_name}: generated topology counts drifted between models")
        model_reports[model_name] = {
            "frame_count": 3,
            "finite_components": True,
            "component_values": component_values,
            "peak_rss_bytes_per_rank": parse_peak_rss_bytes(time_record.read_text(encoding="utf-8"), mpi_ranks=3),
            "decomposition_status": "pass",
            "decomposition_output_sha256": file_sha256(decomp_csv),
            "decomposition_residues_in_frozen_order": decomp_validation[
                "residues_in_frozen_order"
            ],
            "generated_topology_sha256": {label: file_sha256(path) for label, path in topology_paths.items()},
            "generated_topology_properties": properties,
            "command_argv": command,
            "returncode": completed.returncode,
            "outputs": {
                final_text.name: {"sha256": file_sha256(final_text), "bytes": final_text.stat().st_size},
                final_csv.name: {"sha256": file_sha256(final_csv), "bytes": final_csv.stat().st_size},
                log_path.name: {"sha256": file_sha256(log_path), "bytes": log_path.stat().st_size},
                decomp_text.name: {"sha256": file_sha256(decomp_text), "bytes": decomp_text.stat().st_size},
                decomp_csv.name: {"sha256": file_sha256(decomp_csv), "bytes": decomp_csv.stat().st_size},
            },
        }
        if time_record.is_file():
            model_reports[model_name]["outputs"][time_record.name] = {
                "sha256": file_sha256(time_record),
                "bytes": time_record.stat().st_size,
            }
    assert reference_counts is not None
    report = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_three_frame_canary",
        "status": "pass",
        "models": model_reports,
        "complex_atom_count": reference_counts["complex"]["atom_count"],
        "receptor_atom_count": reference_counts["receptor"]["atom_count"],
        "ligand_atom_count": reference_counts["ligand"]["atom_count"],
        "ligand_total_charge": reference_counts["ligand"]["total_charge"],
        "ligand_partial_charges_all_zero": reference_counts["ligand"]["partial_charges_all_zero"],
        "toolchain_status": "pass",
        "sources": {
            "preparation_manifest_sha256": file_sha256(prep_manifest_path),
            "geometry_sha256": file_sha256(geometry_path),
            "toolchain_report_sha256": file_sha256(toolchain_report_path),
            "execution_defaults_sha256": file_sha256(execution_defaults_path),
            "execution_supplement_sha256": file_sha256(execution_supplement_path),
            "native_contact_set_sha256": file_sha256(native_contact_set_path),
        },
    }
    validate_canary_report(report)
    report_path = output_dir / "CANARY_REPORT.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(report_path)
    report_path.with_suffix(".json.sha256").write_text(f"{digest}  {report_path.name}\n", encoding="ascii")
    return {"status": "pass", "report": str(report_path), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-inputs", action="store_true")
    group.add_argument("--validate-report", type=Path)
    group.add_argument("--execute", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--mthick-angstrom", type=float)
    parser.add_argument("--frame-count", type=int, default=3)
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--prep-dir", type=Path)
    parser.add_argument("--topology-bundle", type=Path)
    parser.add_argument("--topology-relative-path", type=Path, default=Path("topol.top"))
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--toolchain-report", type=Path)
    parser.add_argument("--execution-defaults", type=Path)
    parser.add_argument("--execution-supplement", type=Path)
    parser.add_argument("--native-contact-set", type=Path)
    args = parser.parse_args()
    if args.write_inputs:
        if args.output_dir is None or args.mthick_angstrom is None:
            parser.error("--write-inputs requires --output-dir and --mthick-angstrom")
        result = write_model_inputs(args.output_dir, args.mthick_angstrom, args.frame_count)
    elif args.validate_report is not None:
        report = json.loads(args.validate_report.read_text(encoding="utf-8"))
        result = validate_canary_report(report)
    else:
        required = {
            "--prefix": args.prefix,
            "--prep-dir": args.prep_dir,
            "--topology-bundle": args.topology_bundle,
            "--geometry": args.geometry,
            "--toolchain-report": args.toolchain_report,
            "--execution-defaults": args.execution_defaults,
            "--execution-supplement": args.execution_supplement,
            "--native-contact-set": args.native_contact_set,
            "--output-dir": args.output_dir,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            parser.error("--execute requires " + ", ".join(missing))
        result = execute_canary(
            args.prefix,
            args.prep_dir,
            args.topology_bundle,
            args.topology_relative_path,
            args.geometry,
            args.toolchain_report,
            args.execution_defaults,
            args.execution_supplement,
            args.native_contact_set,
            args.output_dir,
        )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
