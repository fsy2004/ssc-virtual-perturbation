#!/usr/bin/env python3
"""Run one fail-closed O6U MP2/6-31G(d) geometry-optimization canary.

Only a previously selected, identity-audited conformer may enter this job.  A
successful run releases later representative calculations for review; it does
not approve charges, bonded terms, the force field, or production MD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

from prepare_o6u_crest_input import EXPECTED_ATOMS, load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble


OPTIMIZER_STRATEGY_VERSION = "cartesian_rfo_trust020_v1"
OPTIMIZER_OPTIONS = {
    # OptKing 0.5.0 documents Cartesian coordinates as the mature fallback
    # when redundant-internal-coordinate back-transformation is unstable.
    "opt_coordinates": "cartesian",
    "step_type": "rfo",
    "intrafrag_step_limit": 0.20,
    "intrafrag_step_limit_min": 0.01,
    "intrafrag_step_limit_max": 0.20,
    "dynamic_level": 0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_record(path: Path, record: dict[str, object]) -> None:
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--start-xyz", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--memory-gib", type=int, default=96)
    parser.add_argument("--scratch-root", type=Path, default=Path("/root/autodl-tmp/psi4_scratch"))
    parser.add_argument(
        "--role",
        choices=("canary", "representative_target"),
        default="canary",
        help="Scientific role of this optimization; the model chemistry is identical for both roles.",
    )
    args = parser.parse_args()

    if not 1 <= args.threads <= 24:
        raise SystemExit("This project canary requires 1-24 threads")
    if not 16 <= args.memory_gib <= 128:
        raise SystemExit("This project canary requires 16-128 GiB of declared Psi4 memory")
    source_path = args.source_sdf.resolve()
    start_path = args.start_xyz.resolve()
    output_dir = args.output_dir.resolve()
    for path in (source_path, start_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty input: {path}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    source = load_single_sdf(source_path)
    identity = validate_identity(source)
    frames = read_xyz_ensemble(start_path)
    if len(frames) != 1:
        raise SystemExit("MP2 canary start must contain exactly one XYZ frame")
    frame = frames[0]
    expected_elements = [atom.GetSymbol() for atom in source.GetAtoms()]
    if frame["elements"] != expected_elements or len(expected_elements) != EXPECTED_ATOMS:
        raise SystemExit("MP2 canary XYZ element sequence differs from immutable O6U SDF")
    coordinates = frame["coordinates"]
    if not all(math.isfinite(float(value)) for value in coordinates.reshape(-1)):
        raise SystemExit("MP2 canary XYZ contains non-finite coordinates")

    scratch = (args.scratch_root.resolve() / f"o6u_mp2_canary_{os.getpid()}")
    scratch.mkdir(parents=True, exist_ok=False)
    os.environ["PSI_SCRATCH"] = str(scratch)

    import psi4

    artifact_stem = (
        "o6u_mp2_631gd_optimization_canary"
        if args.role == "canary"
        else "o6u_mp2_631gd_optimization_representative_target"
    )
    raw_output = output_dir / f"{artifact_stem}.psi4.out"
    optimized_xyz = output_dir / f"{artifact_stem}.optimized.xyz"
    record_path = output_dir / f"{artifact_stem}.json"
    geometry_lines = ["0 1"]
    geometry_lines.extend(
        f"{element:<2s} {float(x): .10f} {float(y): .10f} {float(z): .10f}"
        for element, (x, y, z) in zip(expected_elements, coordinates, strict=True)
    )
    geometry_lines.extend(["units angstrom", "no_com", "no_reorient", "symmetry c1"])
    molecule = psi4.geometry("\n".join(geometry_lines))

    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gib} GiB")
    psi4.core.set_output_file(str(raw_output), False)
    psi4.set_options(
        {
            "basis": "6-31G(d)",
            "reference": "rhf",
            "scf_type": "df",
            "mp2_type": "df",
            "freeze_core": True,
            "guess": "sad",
            "e_convergence": 1.0e-8,
            "d_convergence": 1.0e-8,
            "maxiter": 200,
            "geom_maxiter": 100,
            "g_convergence": "gau_tight",
            **OPTIMIZER_OPTIONS,
        }
    )
    record: dict[str, object] = {
        "schema_version": "1.1",
        "report_type": f"o6u_mp2_631gd_optimization_{args.role}",
        "role": args.role,
        "status": "running",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "psi4_version": psi4.__version__,
        "identity": identity,
        "charge_e": 0,
        "multiplicity": 1,
        "method": "DF-MP2/6-31G(d), frozen core, RHF reference",
        "optimizer_convergence": "gau_tight",
        "optimizer_strategy_version": OPTIMIZER_STRATEGY_VERSION,
        "optimizer_options": OPTIMIZER_OPTIONS,
        "threads": args.threads,
        "memory_gib": args.memory_gib,
        "source_sdf": {"path": str(source_path), "sha256": sha256(source_path)},
        "start_xyz": {"path": str(start_path), "sha256": sha256(start_path), "comment": str(frame["comment"])},
        "raw_output": str(raw_output),
        "scratch": str(scratch),
        "pid": os.getpid(),
        "release_boundary": (
            "Canary only; no force-field parameter or production MD is approved."
            if args.role == "canary"
            else "QM representative target only; no fitted force-field parameter or production MD is approved."
        ),
    }
    write_record(record_path, record)
    started = time.time()
    exit_code = 1
    try:
        energy, wavefunction = psi4.optimize("mp2", molecule=molecule, return_wfn=True)
        molecule.save_xyz_file(str(optimized_xyz), True)
        record.update(
            {
                "status": (
                    "pass_optimization_canary"
                    if args.role == "canary"
                    else "pass_optimization_representative_target"
                ),
                "final_energy_hartree": float(energy),
                "wavefunction_energy_hartree": float(wavefunction.energy()),
                "optimized_xyz": {"path": str(optimized_xyz), "sha256": sha256(optimized_xyz)},
            }
        )
        exit_code = 0
    except BaseException as exc:
        record.update({"status": "fail", "error_type": type(exc).__name__, "error": str(exc)})
    finally:
        record["elapsed_seconds"] = time.time() - started
        psi4.core.flush_outfile()
        if raw_output.is_file():
            record["raw_output_sha256"] = sha256(raw_output)
            record["raw_output_bytes"] = raw_output.stat().st_size
        write_record(record_path, record)
        print(json.dumps({"status": record["status"], "elapsed_seconds": record["elapsed_seconds"], "record_sha256": sha256(record_path)}, sort_keys=True))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
