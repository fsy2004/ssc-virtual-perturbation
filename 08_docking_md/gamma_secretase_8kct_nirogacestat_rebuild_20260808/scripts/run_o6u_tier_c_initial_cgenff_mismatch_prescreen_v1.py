#!/usr/bin/env python3
"""Run a scope-correct Tier-C same-geometry initial-CGenFF mismatch prescreen.

The numerical evaluator is the previously validated Tier-B prescreen engine.
This wrapper changes only scope metadata, preserves all bound inputs and
numeric results, and removes the temporary engine report after hashing it.
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

    allowed_rotors = {
        "ROT_C20_C18",
        "ROT_C17_N06",
        "ROT_C26_C19",
        "ROT_C14_C19",
        "ROT_N06_C20",
    }
    if args.rotor_id not in allowed_rotors:
        raise ValueError(f"Not a Tier-C rotor: {args.rotor_id}")

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
    with tempfile.TemporaryDirectory(prefix="o6u_tier_c_prescreen_") as temporary:
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

    if report.get("status") != "pass_tier_b_initial_cgenff_mismatch_prescreen":
        raise ValueError("Numerical prescreen engine did not pass")
    if report.get("scope", {}).get("rotor_id") != args.rotor_id:
        raise ValueError("Rotor identity mismatch in numerical engine output")
    if report.get("scope", {}).get("signed_step_index") != args.signed_step_index:
        raise ValueError("Signed-step mismatch in numerical engine output")

    report["schema_version"] = "1.0"
    report["report_type"] = "o6u_tier_c_initial_cgenff_mismatch_prescreen"
    report["status"] = "pass_tier_c_initial_cgenff_mismatch_prescreen"
    report["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    report["scope"] = {
        "compound": "nirogacestat/PF-03084014/BRD-K61691541",
        "residue": "O6U",
        "tier": "C conditional transferability check",
        "rotor_id": args.rotor_id,
        "signed_step_index": args.signed_step_index,
        "purpose": (
            "Decide whether the completed displaced point demonstrates an "
            "initial-CGenFF energy mismatch that justifies minimum additional QM."
        ),
    }
    report["numerical_engine"] = {
        "path": str(engine.resolve()),
        "sha256": sha256(engine),
        "temporary_engine_report_sha256": engine_report_sha,
        "metadata_override_only": True,
    }
    report["interpretation_boundary"] = (
        "This same-geometry initial-force-field comparison is a conditional scan-scope "
        "prescreen. It does not fit or validate parameters, affinity, efficacy, or "
        "production MD readiness. Additional QM requires a documented mismatch trigger."
    )

    output_dir.mkdir()
    report_path = output_dir / "O6U_TIER_C_INITIAL_CGENFF_MISMATCH_PRESCREEN.json"
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
