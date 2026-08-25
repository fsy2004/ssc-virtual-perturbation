#!/usr/bin/env python3
"""Generate one hash-bound initial-CGenFF fitting input from a passed scan point.

The numerical calculation is delegated unchanged to the validated Tier-B
same-geometry evaluator.  This wrapper only constrains the closed adaptive
dataset scope and labels the result as a preliminary fitting input.  It does
not mutate, fit, or validate force-field parameters.
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


ALLOWED_STEPS = {
    "ROT_C09_N04": {-2, -1, 1, 2},
    "ROT_C15_N05": {-3, -2, -1, 1, 2, 3},
    "ROT_C17_C15": {-3, -2, -1, 1, 2, 3},
    "ROT_C24_C14": {-3, -2, -1, 1, 2, 3},
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_scope(rotor_id: str, signed_step_index: int) -> None:
    if rotor_id not in ALLOWED_STEPS:
        raise ValueError(f"Rotor is outside the closed adaptive fitting dataset: {rotor_id}")
    if signed_step_index not in ALLOWED_STEPS[rotor_id]:
        raise ValueError(f"Step is outside the authorized adaptive scope: {rotor_id} {signed_step_index}")


def reframe_report(
    report: dict,
    *,
    rotor_id: str,
    signed_step_index: int,
    engine_path: Path,
    engine_sha: str,
    engine_report_sha: str,
) -> dict:
    validate_scope(rotor_id, signed_step_index)
    if report.get("status") != "pass_tier_b_initial_cgenff_mismatch_prescreen":
        raise ValueError("Numerical evaluator did not pass")
    scope = report.get("scope", {})
    if scope.get("rotor_id") != rotor_id or scope.get("signed_step_index") != signed_step_index:
        raise ValueError("Numerical evaluator scope mismatch")

    reframed = dict(report)
    reframed["schema_version"] = "1.0"
    reframed["report_type"] = "o6u_adaptive_scan_initial_cgenff_fitting_input"
    reframed["status"] = "pass_adaptive_scan_initial_cgenff_fitting_input"
    reframed["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    reframed["scope"] = {
        "compound": "nirogacestat/PF-03084014/BRD-K61691541",
        "residue": "O6U",
        "dataset": "closed compound-specific adaptive scan",
        "rotor_id": rotor_id,
        "signed_step_index": signed_step_index,
        "purpose": "Preliminary quantitative fitting input and mismatch triage",
    }
    reframed["method"] = dict(reframed.get("method", {}))
    reframed["method"]["parameter_mutation"] = False
    reframed["numerical_engine"] = {
        "path": str(engine_path.resolve()),
        "sha256": engine_sha,
        "temporary_engine_report_sha256": engine_report_sha,
        "metadata_override_only": True,
    }
    reframed["interpretation_boundary"] = (
        "This is a same-geometry initial-CGenFF residual used as a preliminary fitting input. "
        "It is not a fitted parameter set, an independent parameter validation, or approval "
        "for production MD, affinity, or efficacy claims."
    )
    return reframed


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
    parser.add_argument("--displaced-qm-report", required=True)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--signed-step-index", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    validate_scope(args.rotor_id, args.signed_step_index)

    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    engine = Path(__file__).with_name("run_o6u_tier_b_initial_cgenff_mismatch_prescreen_v2.py")
    if not engine.is_file():
        raise FileNotFoundError(engine)

    forwarded = [
        "--psf", args.psf,
        "--base-rtf", args.base_rtf,
        "--base-prm", args.base_prm,
        "--ligand-rtf", args.ligand_rtf,
        "--ligand-prm", args.ligand_prm,
        "--reference-xyz", args.reference_xyz,
        "--displaced-xyz", args.displaced_xyz,
        "--reference-qm-energy-hartree", str(args.reference_qm_energy_hartree),
        "--displaced-qm-report", args.displaced_qm_report,
        "--rotor-id", args.rotor_id,
        "--signed-step-index", str(args.signed_step_index),
    ]
    with tempfile.TemporaryDirectory(prefix="o6u_adaptive_fitting_input_") as temporary:
        temporary_output = Path(temporary) / "engine"
        subprocess.run(
            [sys.executable, str(engine), *forwarded, "--output-dir", str(temporary_output)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        engine_report = temporary_output / "O6U_TIER_B_INITIAL_CGENFF_MISMATCH_PRESCREEN.json"
        report = json.loads(engine_report.read_text())
        engine_report_sha = sha256(engine_report)

    report = reframe_report(
        report,
        rotor_id=args.rotor_id,
        signed_step_index=args.signed_step_index,
        engine_path=engine,
        engine_sha=sha256(engine),
        engine_report_sha=engine_report_sha,
    )
    output_dir.mkdir()
    report_path = output_dir / "O6U_ADAPTIVE_SCAN_INITIAL_CGENFF_FITTING_INPUT.json"
    temporary_path = report_path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary_path.replace(report_path)
    print(json.dumps({
        "status": report["status"],
        "report": str(report_path),
        "sha256": sha256(report_path),
        "comparison": report["comparison"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
