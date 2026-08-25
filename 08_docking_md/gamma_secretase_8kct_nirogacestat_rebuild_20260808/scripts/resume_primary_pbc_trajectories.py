#!/usr/bin/env python3
"""Resume PBC validation from hash-bound retained trajectory intermediates.

This recovery path is intentionally narrower than a fresh preparation.  It
validates the production completion gate and retained trajectories, archives
only interrupted mindist outputs, and then completes the invariance/check/
provenance stages without replacing the valid 03/05/06/07 intermediates.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path
from typing import Any

from make_analysis_trajectories import compare_minimum_image_distances, execute, parse_ndx
from prepare_primary_pbc_trajectories import (
    GROUP_NAMES,
    PRODUCTION_RELEASE_SHA256,
    REPLICAS,
    TPR_SHA256,
    _run_check,
    pipeline_contract,
    sha256,
    utc_now,
    validate_completion_gate,
)


RETAINED_INPUTS = (
    "03_processed_first_frame.gro",
    "05_centered_reboxed.xtc",
    "06_fitted_analysis.xtc",
    "07_fixed_200_500ns.xtc",
)
PARTIAL_MINDIST_OUTPUTS = (
    "09_raw_minimum_image_protein_O6U_heavy.xvg",
    "10_processed_minimum_image_protein_O6U_heavy.xvg",
    "mindist_raw.command.json",
    "mindist_processed.command.json",
)
FINAL_OUTPUTS = (
    "11_pbc_distance_invariance.json",
    "trajectory_provenance.pre_qc.json",
    "trajectory_provenance.pre_qc.json.sha256",
)


def validate_retained_inputs(output_dir: Path) -> dict[str, dict[str, Any]]:
    final = [name for name in FINAL_OUTPUTS if (output_dir / name).exists()]
    if any(name.startswith("trajectory_provenance") for name in final):
        raise FileExistsError(f"refusing to replace final provenance: {final}")
    if "11_pbc_distance_invariance.json" in final:
        raise FileExistsError(f"refusing to replace final PBC output: {final}")
    records: dict[str, dict[str, Any]] = {}
    for name in RETAINED_INPUTS:
        path = output_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing retained PBC intermediate: {path}")
        records[name] = {"path": str(path), "bytes": path.stat().st_size}
    return records


def archive_partial_mindist_outputs(
    output_dir: Path, archive_dir: Path
) -> list[dict[str, Any]]:
    final = [name for name in FINAL_OUTPUTS if (output_dir / name).exists()]
    if final:
        raise FileExistsError(f"refusing to archive beside final PBC output: {final}")
    candidates = [output_dir / name for name in PARTIAL_MINDIST_OUTPUTS]
    present = [path for path in candidates if path.exists()]
    if not present:
        return []
    archive_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for source in present:
        if not source.is_file():
            raise ValueError(f"partial output is not a regular file: {source}")
        record = {
            "name": source.name,
            "source": str(source),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
        }
        target = archive_dir / source.name
        shutil.move(str(source), str(target))
        record["archived_to"] = str(target)
        records.append(record)
    return records


def resume(root: Path, replica: str, gmx: str, recovery_id: str) -> dict[str, Any]:
    root = root.resolve()
    output_dir = (
        root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / replica
    )
    retained_inputs = validate_retained_inputs(output_dir)
    completion_path = root / replica / "PRODUCTION_COMPLETION_500NS.json"
    gate = validate_completion_gate(root, replica, completion_path)
    ndx = root / "builds" / "analysis.ndx"
    groups = parse_ndx(ndx)
    missing = [name for name in GROUP_NAMES.values() if name not in groups]
    if missing:
        raise ValueError(f"analysis index lacks frozen groups: {missing}")

    paths = {
        "processed_first_frame_reference": output_dir / RETAINED_INPUTS[0],
        "center_and_rebox": output_dir / RETAINED_INPUTS[1],
        "fit_analysis_selection": output_dir / RETAINED_INPUTS[2],
        "fixed_window_extract": output_dir / RETAINED_INPUTS[3],
    }
    for label, key in (
        ("centered_resume_preflight", "center_and_rebox"),
        ("fitted_resume_preflight", "fit_analysis_selection"),
        ("fixed_resume_preflight", "fixed_window_extract"),
    ):
        _run_check(gmx, paths[key], output_dir, label)

    recovery_dir = root / "audit" / "pbc_resume" / recovery_id / replica
    archived = archive_partial_mindist_outputs(output_dir, recovery_dir)
    recovery_manifest = {
        "schema_version": "1.0",
        "report_type": "pbc_resume_partial_archive",
        "status": "archived_before_resume",
        "replica": replica,
        "recovery_id": recovery_id,
        "created_at_utc": utc_now(),
        "retained_inputs": retained_inputs,
        "archived_partial_outputs": archived,
    }
    manifest_path = recovery_dir / "RECOVERY_ARCHIVE.json"
    manifest_path.write_text(json.dumps(recovery_manifest, indent=2) + "\n", encoding="utf-8")
    manifest_path.with_suffix(".json.sha256").write_text(
        f"{sha256(manifest_path)}  {manifest_path.name}\n", encoding="ascii"
    )

    work = root / replica / "work"
    tpr = work / "production.tpr"
    raw = work / "production.xtc"
    raw_xvg = output_dir / "09_raw_minimum_image_protein_O6U_heavy.xvg"
    processed_xvg = output_dir / "10_processed_minimum_image_protein_O6U_heavy.xvg"
    selections = [
        groups[GROUP_NAMES["pbc_invariance_protein_heavy"]],
        groups[GROUP_NAMES["pbc_invariance_ligand_heavy"]],
    ]
    for label, trajectory, target in (
        ("raw", raw, raw_xvg),
        ("processed", paths["center_and_rebox"], processed_xvg),
    ):
        execute(
            [
                gmx,
                "mindist",
                "-s",
                str(tpr),
                "-f",
                str(trajectory),
                "-n",
                str(ndx),
                "-od",
                str(target),
                "-b",
                "0.000",
                "-e",
                "500000.000",
            ],
            selections,
            output_dir,
            output_dir / f"mindist_{label}.command.json",
        )

    pbc = compare_minimum_image_distances(raw_xvg, processed_xvg, 0.01)
    pbc.update(
        {
            "schema_version": "1.0",
            "report_type": "pbc_minimum_image_invariance",
            "realization_id": replica,
            "analysis_window_ns": [200.0, 500.0],
            "evaluated_window_ns": [0.0, 500.0],
            "tolerance_nm": 0.01,
            "completion_gate": gate,
            "resume_recovery_manifest": {
                "path": str(manifest_path),
                "sha256": sha256(manifest_path),
            },
            "created_at_utc": utc_now(),
        }
    )
    pbc_path = output_dir / "11_pbc_distance_invariance.json"
    pbc_path.write_text(json.dumps(pbc, indent=2) + "\n", encoding="utf-8")
    if pbc.get("status") != "pass":
        raise RuntimeError(f"PBC distance invariance failed for {replica}")

    for label, key in (
        ("centered", "center_and_rebox"),
        ("fitted", "fit_analysis_selection"),
        ("fixed", "fixed_window_extract"),
    ):
        _run_check(gmx, paths[key], output_dir, label)

    retained = {
        name: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        for name, path in paths.items()
    }
    retained.update(
        {
            "raw_mindist": {
                "path": str(raw_xvg),
                "bytes": raw_xvg.stat().st_size,
                "sha256": sha256(raw_xvg),
            },
            "processed_mindist": {
                "path": str(processed_xvg),
                "bytes": processed_xvg.stat().st_size,
                "sha256": sha256(processed_xvg),
            },
            "pbc_report": {
                "path": str(pbc_path),
                "bytes": pbc_path.stat().st_size,
                "sha256": sha256(pbc_path),
            },
        }
    )
    provenance = {
        "schema_version": "2.1",
        "report_type": "primary_pbc_trajectory_preparation",
        "status": "pass_pending_scientific_qc_seal",
        "realization_id": replica,
        "raw_files_immutable": True,
        "production_release_sha256": PRODUCTION_RELEASE_SHA256,
        "production_tpr_sha256": TPR_SHA256[replica],
        "source_completion_gate": gate,
        "source_xtc": {
            "path": str(raw),
            "bytes": raw.stat().st_size,
            "sha256": gate["source_xtc"]["sha256"],
            "hash_source": "production_500ns_completion_report",
        },
        "analysis_index": {"path": str(ndx), "sha256": sha256(ndx)},
        "pipeline": pipeline_contract(),
        "resume_recovery_manifest": {
            "path": str(manifest_path),
            "sha256": sha256(manifest_path),
        },
        "retained_outputs": retained,
        "downstream_endpoint_eligible": False,
        "created_at_utc": utc_now(),
    }
    provenance_path = output_dir / "trajectory_provenance.pre_qc.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    digest = sha256(provenance_path)
    provenance_path.with_suffix(".json.sha256").write_text(
        f"{digest}  {provenance_path.name}\n", encoding="ascii"
    )
    return {
        "status": provenance["status"],
        "replica": replica,
        "provenance": str(provenance_path),
        "sha256": digest,
        "recovery_manifest": str(manifest_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", required=True, type=Path)
    parser.add_argument("--replica", required=True, choices=REPLICAS)
    parser.add_argument("--recovery-id", required=True)
    parser.add_argument("--gmx", default="/root/GROMACS-2025.2/bin/gmx")
    args = parser.parse_args()
    if not Path(args.gmx).is_file() or not os.access(args.gmx, os.X_OK):
        raise SystemExit(f"GROMACS executable is missing or not executable: {args.gmx}")
    print(
        json.dumps(
            resume(args.release_root, args.replica, args.gmx, args.recovery_id),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
