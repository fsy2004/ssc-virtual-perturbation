#!/usr/bin/env python3
"""Prepare formal deterministic PBC-safe trajectories before scientific QC.

This resolves the legacy circular gate in which structural/membrane
stationarity was required before the trajectories needed to evaluate it had
been created.  Raw production files are opened read-only.  The exact frozen
PBC order and 200--500 ns window are unchanged; outputs remain ineligible for
endpoint analysis until the later all-three QC seal passes.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any

from make_analysis_trajectories import (
    command,
    compare_minimum_image_distances,
    execute,
    parse_ndx,
)


REPLICAS = ("rep01", "rep02", "rep03")
PRODUCTION_RELEASE_SHA256 = "a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa"
TPR_SHA256 = {
    "rep01": "fd11c7287d5670c81ccb44fcb5b4215344726989f66f4f55db33643ba618678f",
    "rep02": "887f273a06cf4414692589584479d892cd1e5c0054b49e375a392ee852307d1a",
    "rep03": "abe512a4971c6cc26a61c3d9fbea39df8b40074265e1d30b11a488bd3ffac9ad",
}
REQUIRED_ARTIFACTS = (
    "production.tpr",
    "production.xtc",
    "production.edr",
    "production.log",
    "production.gro",
    "production.cpt",
)
GROUP_NAMES = {
    "system": "System",
    "complex": "Protein_O6U",
    "fit": "PSEN1_Core",
    "analysis_output": "Protein_O6U",
    "pbc_invariance_protein_heavy": "Protein_Heavy",
    "pbc_invariance_ligand_heavy": "O6U_Heavy",
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def pipeline_contract() -> list[dict[str, Any]]:
    return [
        {"mode": "whole", "retain": False},
        {"mode": "cluster_complex_if_required", "retain": False},
        {"mode": "processed_first_frame_reference", "retain": True},
        {"mode": "nojump", "retain": False},
        {"mode": "center_and_rebox", "retain": True},
        {"mode": "fit_analysis_selection", "retain": True},
        {"mode": "fixed_window_extract", "retain": True},
    ]


def validate_completion_gate(root: Path, replica: str, report_path: Path) -> dict[str, Any]:
    if replica not in REPLICAS:
        raise ValueError(f"unknown replica: {replica}")
    report = load_json(report_path)
    expected = {
        "schema_version": "1.0",
        "report_type": "production_500ns_completion",
        "status": "pass",
        "replica": replica,
        "final_step": 125000000,
        "final_time_ps": 500000.0,
        "production_release_sha256": PRODUCTION_RELEASE_SHA256,
        "production_tpr_sha256": TPR_SHA256.get(replica, report.get("production_tpr_sha256")),
    }
    for key, value in expected.items():
        if report.get(key) != value:
            raise ValueError(f"{replica}: completion report mismatch at {key}")
    readability = report.get("checks", {}).get("gmx_readability", {})
    if any(readability.get(key) != "pass" for key in ("cpt", "edr", "gro", "xtc")):
        raise ValueError(f"{replica}: completion report lacks passing GROMACS readability")
    artifacts = report.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError(f"{replica}: completion artifacts are invalid")
    work = root / replica / "work"
    for name in REQUIRED_ARTIFACTS:
        record = artifacts.get(name)
        path = work / name
        if not isinstance(record, dict) or not path.is_file():
            raise ValueError(f"{replica}: missing completion artifact {name}")
        if path.stat().st_size != record.get("bytes"):
            raise ValueError(f"{replica}: live byte count differs for {name}")
        if name != "production.xtc" and sha256(path) != record.get("sha256"):
            raise ValueError(f"{replica}: live SHA-256 differs for {name}")
    if artifacts["production.tpr"]["sha256"] != expected["production_tpr_sha256"]:
        raise ValueError(f"{replica}: completion TPR artifact hash differs")
    return {
        "status": "pass",
        "report": {"path": str(report_path), "sha256": sha256(report_path)},
        "source_xtc": artifacts["production.xtc"],
    }


def _run_check(gmx: str, path: Path, output_dir: Path, label: str) -> dict[str, Any]:
    result = subprocess.run(
        [gmx, "check", "-f", str(path)],
        cwd=output_dir,
        text=True,
        capture_output=True,
        errors="replace",
        check=False,
    )
    record = {
        "argv": [gmx, "check", "-f", str(path)],
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "checked_at_utc": utc_now(),
    }
    (output_dir / f"check_{label}.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    if result.returncode != 0:
        raise RuntimeError(f"gmx check failed for {path}")
    return record


def prepare(root: Path, replica: str, gmx: str) -> dict[str, Any]:
    root = root.resolve()
    completion_path = root / replica / "PRODUCTION_COMPLETION_500NS.json"
    gate = validate_completion_gate(root, replica, completion_path)
    ndx = root / "builds" / "analysis.ndx"
    groups = parse_ndx(ndx)
    missing = [name for name in GROUP_NAMES.values() if name not in groups]
    if missing:
        raise ValueError(f"analysis index lacks frozen groups: {missing}")
    work = root / replica / "work"
    tpr = work / "production.tpr"
    raw = work / "production.xtc"
    output_dir = root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / replica
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite derived trajectory directory: {output_dir}")
    output_dir.mkdir(parents=True)
    paths = {
        "whole": output_dir / "01_whole.xtc",
        "cluster_complex_if_required": output_dir / "02_clustered_complex.xtc",
        "processed_first_frame_reference": output_dir / "03_processed_first_frame.gro",
        "nojump": output_dir / "04_nojump.xtc",
        "center_and_rebox": output_dir / "05_centered_reboxed.xtc",
        "fit_analysis_selection": output_dir / "06_fitted_analysis.xtc",
        "fixed_window_extract": output_dir / "07_fixed_200_500ns.xtc",
    }
    steps = [
        ("whole", raw, paths["whole"], [groups[GROUP_NAMES["system"]]]),
        ("cluster_complex_if_required", paths["whole"], paths["cluster_complex_if_required"], [groups[GROUP_NAMES["complex"]], groups[GROUP_NAMES["system"]]]),
        ("processed_first_frame_reference", paths["cluster_complex_if_required"], paths["processed_first_frame_reference"], [groups[GROUP_NAMES["system"]]]),
        ("nojump", paths["cluster_complex_if_required"], paths["nojump"], [groups[GROUP_NAMES["system"]]]),
        ("center_and_rebox", paths["nojump"], paths["center_and_rebox"], [groups[GROUP_NAMES["complex"]], groups[GROUP_NAMES["system"]]]),
        ("fit_analysis_selection", paths["center_and_rebox"], paths["fit_analysis_selection"], [groups[GROUP_NAMES["fit"]], groups[GROUP_NAMES["analysis_output"]]]),
        ("fixed_window_extract", paths["fit_analysis_selection"], paths["fixed_window_extract"], [0]),
    ]
    delete_after = {
        "cluster_complex_if_required": paths["whole"],
        "nojump": paths["cluster_complex_if_required"],
        "center_and_rebox": paths["nojump"],
    }
    for index, (mode, source, output, selections) in enumerate(steps, start=1):
        argv = command(gmx, tpr, source, output, ndx, mode, (200.0, 500.0), paths["processed_first_frame_reference"])
        execute(argv, selections, output_dir, output_dir / f"{index:02d}_{mode}.command.json")
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"derived output missing after {mode}: {output}")
        obsolete = delete_after.get(mode)
        if obsolete is not None:
            obsolete.unlink()

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
            [gmx, "mindist", "-s", str(tpr), "-f", str(trajectory), "-n", str(ndx), "-od", str(target), "-b", "0.000", "-e", "500000.000"],
            selections,
            output_dir,
            output_dir / f"mindist_{label}.command.json",
        )
    pbc = compare_minimum_image_distances(raw_xvg, processed_xvg, 0.01)
    pbc.update({
        "schema_version": "1.0",
        "report_type": "pbc_minimum_image_invariance",
        "realization_id": replica,
        "analysis_window_ns": [200.0, 500.0],
        "evaluated_window_ns": [0.0, 500.0],
        "tolerance_nm": 0.01,
        "completion_gate": gate,
        "created_at_utc": utc_now(),
    })
    pbc_path = output_dir / "11_pbc_distance_invariance.json"
    pbc_path.write_text(json.dumps(pbc, indent=2) + "\n", encoding="utf-8")
    if pbc.get("status") != "pass":
        raise RuntimeError(f"PBC distance invariance failed for {replica}")

    for label in ("centered", "fitted", "fixed"):
        _run_check(
            gmx,
            paths[{"centered": "center_and_rebox", "fitted": "fit_analysis_selection", "fixed": "fixed_window_extract"}[label]],
            output_dir,
            label,
        )
    retained = {
        name: {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for name, path in paths.items()
        if path.is_file()
    }
    retained.update({
        "raw_mindist": {"path": str(raw_xvg), "bytes": raw_xvg.stat().st_size, "sha256": sha256(raw_xvg)},
        "processed_mindist": {"path": str(processed_xvg), "bytes": processed_xvg.stat().st_size, "sha256": sha256(processed_xvg)},
        "pbc_report": {"path": str(pbc_path), "bytes": pbc_path.stat().st_size, "sha256": sha256(pbc_path)},
    })
    provenance = {
        "schema_version": "2.0",
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
    return {"status": provenance["status"], "replica": replica, "provenance": str(provenance_path), "sha256": digest}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", required=True, type=Path)
    parser.add_argument("--replica", required=True, choices=REPLICAS)
    parser.add_argument("--gmx", default="/root/GROMACS-2025.2/bin/gmx")
    args = parser.parse_args()
    if not Path(args.gmx).is_file() or not os.access(args.gmx, os.X_OK):
        raise SystemExit(f"GROMACS executable is missing or not executable: {args.gmx}")
    print(json.dumps(prepare(args.release_root, args.replica, args.gmx), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
