#!/usr/bin/env python3
"""Evaluate the bound C17 +15 geometric-recovery point as a fitting input.

The validated numerical engine accepts the ordinary scan schema. This wrapper
creates a temporary schema adapter after verifying the closed fitting manifest,
records both hashes, and restores the original report binding in the final
output. No energy, geometry, method, or parameter is changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    return {"path": str(resolved), "sha256": sha256(resolved), "size_bytes": resolved.stat().st_size}


def adapt_recovery_report(source: dict) -> dict:
    observed = (source.get("status"), source.get("rotor_id"), source.get("signed_step_index"))
    if observed[0] != "pass_relaxed_mp2_torsion_geometric_recovery" or observed[1] != "ROT_C17_C15":
        raise ValueError(f"Unexpected geometric-recovery scope: {observed}")
    if observed[2] not in (None, 1):
        raise ValueError(f"Unexpected geometric-recovery step: {observed[2]}")
    adapted = dict(source)
    adapted["status"] = "pass_relaxed_mp2_torsion_scan_point"
    adapted["signed_step_index"] = 1
    return adapted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psf", required=True)
    parser.add_argument("--base-rtf", required=True)
    parser.add_argument("--base-prm", required=True)
    parser.add_argument("--ligand-rtf", required=True)
    parser.add_argument("--ligand-prm", required=True)
    parser.add_argument("--reference-xyz", required=True)
    parser.add_argument("--displaced-xyz", required=True)
    parser.add_argument("--reference-qm-energy-hartree", type=float, required=True)
    parser.add_argument("--recovery-report", required=True)
    parser.add_argument("--closed-fitting-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    recovery_path = Path(args.recovery_report).resolve()
    xyz_path = Path(args.displaced_xyz).resolve()
    manifest_path = Path(args.closed_fitting_manifest).resolve()
    source = json.loads(recovery_path.read_text())
    adapted = adapt_recovery_report(source)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "pass_hash_bound_adaptive_scan_dataset_closed_for_preliminary_fitting":
        raise ValueError("Closed fitting manifest did not pass")
    matches = [
        point for point in manifest["points"]
        if point["rotor_id"] == "ROT_C17_C15" and int(point["step"]) == 1
    ]
    if len(matches) != 1:
        raise ValueError("Closed fitting manifest does not contain exactly one C17 +15 point")
    bound = matches[0]
    if bound["status"] != "pass_independent_relaxed_mp2_torsion_geometric_recovery":
        raise ValueError("Manifest C17 +15 point is not independently passed geometric recovery")
    if sha256(recovery_path) != bound["scan_sha256"]:
        raise ValueError("Recovery report hash does not match closed fitting manifest")
    if sha256(xyz_path) != bound["optimized_xyz"]["sha256"]:
        raise ValueError("Recovery XYZ hash does not match closed fitting manifest")

    engine = Path(__file__).with_name("run_o6u_tier_b_initial_cgenff_mismatch_prescreen_v2.py")
    if not engine.is_file():
        raise FileNotFoundError(engine)
    with tempfile.TemporaryDirectory(prefix="o6u_c17_plus15_recovery_fitting_") as temporary:
        temporary_root = Path(temporary)
        adapted_path = temporary_root / "status_only_adapted_report.json"
        adapted_path.write_text(json.dumps(adapted, indent=2, sort_keys=True) + "\n")
        engine_output = temporary_root / "engine"
        command = [
            sys.executable, str(engine),
            "--psf", args.psf,
            "--base-rtf", args.base_rtf,
            "--base-prm", args.base_prm,
            "--ligand-rtf", args.ligand_rtf,
            "--ligand-prm", args.ligand_prm,
            "--reference-xyz", args.reference_xyz,
            "--displaced-xyz", args.displaced_xyz,
            "--reference-qm-energy-hartree", str(args.reference_qm_energy_hartree),
            "--displaced-qm-report", str(adapted_path),
            "--rotor-id", "ROT_C17_C15",
            "--signed-step-index", "1",
            "--output-dir", str(engine_output),
        ]
        subprocess.run(command, check=True, text=True, capture_output=True)
        engine_report_path = engine_output / "O6U_TIER_B_INITIAL_CGENFF_MISMATCH_PRESCREEN.json"
        engine_report = json.loads(engine_report_path.read_text())
        engine_report_sha = sha256(engine_report_path)
        adapted_sha = sha256(adapted_path)
    if engine_report.get("status") != "pass_tier_b_initial_cgenff_mismatch_prescreen":
        raise ValueError("Numerical engine did not pass")
    engine_report["schema_version"] = "1.0"
    engine_report["report_type"] = "o6u_geometric_recovery_initial_cgenff_fitting_input"
    engine_report["status"] = "pass_geometric_recovery_initial_cgenff_fitting_input"
    engine_report["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    engine_report["scope"] = {
        "compound": "nirogacestat/PF-03084014/BRD-K61691541",
        "residue": "O6U",
        "dataset": "closed compound-specific adaptive scan",
        "rotor_id": "ROT_C17_C15",
        "signed_step_index": 1,
        "purpose": "Fill the provenance-aware C17 +15 preliminary fitting-input residual",
    }
    engine_report["inputs"]["displaced_qm_report"] = record(recovery_path)
    engine_report["closed_fitting_manifest"] = record(manifest_path)
    changed_fields = ["status"]
    if source.get("signed_step_index") is None:
        changed_fields.append("signed_step_index")
    engine_report["schema_adapter"] = {
        "original_status": source["status"],
        "temporary_status": adapted["status"],
        "original_report_sha256": sha256(recovery_path),
        "temporary_adapted_report_sha256": adapted_sha,
        "signed_step_index_source": "closed independently validated fitting manifest" if source.get("signed_step_index") is None else "original recovery report",
        "changed_fields": changed_fields,
        "energy_geometry_method_changed": False,
    }
    engine_report["numerical_engine"] = {
        "path": str(engine.resolve()),
        "sha256": sha256(engine),
        "temporary_engine_report_sha256": engine_report_sha,
    }
    engine_report["method"]["parameter_mutation"] = False
    engine_report["interpretation_boundary"] = (
        "Same-geometry initial-CGenFF residual for the independently validated geometric-recovery "
        "point only. It is not a fitted parameter set or production-MD approval."
    )
    output_dir.mkdir(parents=True)
    report_path = output_dir / "O6U_GEOMETRIC_RECOVERY_INITIAL_CGENFF_FITTING_INPUT.json"
    temporary_path = report_path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(engine_report, indent=2, sort_keys=True) + "\n")
    temporary_path.replace(report_path)
    print(json.dumps({
        "status": engine_report["status"],
        "comparison": engine_report["comparison"],
        "sha256": sha256(report_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
