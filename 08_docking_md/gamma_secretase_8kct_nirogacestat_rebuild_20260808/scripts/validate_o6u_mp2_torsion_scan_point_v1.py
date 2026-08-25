#!/usr/bin/env python3
"""Independently validate one completed hash-bound O6U MP2 torsion point."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from run_o6u_mp2_torsion_canary_v1 import circular_delta, dihedral, read_xyz, sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--canary-validation", required=True, type=Path)
    parser.add_argument("--scan-report", required=True, type=Path)
    parser.add_argument("--expected-rotor-id", required=True)
    parser.add_argument("--expected-step-index", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    auth_path = args.authorization.resolve()
    canary_path = args.canary_validation.resolve()
    scan_path = args.scan_report.resolve()
    auth = json.loads(auth_path.read_text())
    canary = json.loads(canary_path.read_text())
    scan = json.loads(scan_path.read_text())
    checks: dict[str, bool] = {}

    checks["authorization_status"] = auth.get("status") == "pass_bonded_torsion_scope_authorized_canary_only"
    checks["canary_validation_status"] = canary.get("status") == "pass_independent_relaxed_mp2_torsion_geometric_recovery"
    checks["scan_stage_released"] = canary.get("full_torsion_scan_stage_released") is True
    checks["authorization_hash"] = scan.get("authorization_sha256") == sha256(auth_path)
    checks["canary_validation_hash"] = scan.get("canary_validation_sha256") == sha256(canary_path)
    checks["scan_status"] = scan.get("status") == "pass_relaxed_mp2_torsion_scan_point"
    checks["rotor_identity"] = scan.get("rotor_id") == args.expected_rotor_id
    checks["step_identity"] = scan.get("signed_step_index") == args.expected_step_index
    checks["method"] = scan.get("method") == "DF-MP2/6-31G(d), frozen-core RHF, GAU_TIGHT"
    checks["psi4_version"] = scan.get("psi4_version") == "1.9.1"
    checks["geometric_version"] = scan.get("geometric_version") == "1.1.1"
    checks["optimizer"] = scan.get("optimizer") == "geomeTRIC 1.1.1 TRIC exact-set constraint conmethod=1"
    checks["optimizer_keywords"] = scan.get("optimizer_keywords") == {
        "coordsys": "tric",
        "constraints": {
            "set": [{
                "type": "dihedral",
                "indices": [value - 1 for value in scan.get("one_based_ordinals", [])],
                "value": scan.get("target_angle_deg"),
            }]
        },
        "conmethod": 1,
        "maxiter": 200,
        "convergence_set": "GAU_TIGHT",
    }
    checks["resource_route"] = scan.get("threads") == 24 and scan.get("memory_gib") == 44 and scan.get("scan_parallelism") == 1
    checks["finite_energy"] = math.isfinite(float(scan.get("final_energy_hartree", math.nan)))
    checks["reported_constraint_error"] = float(scan.get("target_error_deg", math.inf)) <= 0.5
    checks["scratch_success_cleanup"] = scan.get("scratch_cleanup") == "removed_after_success"
    checks["not_production_approved"] = scan.get("production_approved") is False

    raw_path = Path(scan["raw_output"])
    xyz_path = Path(scan["optimized_xyz"]["path"])
    start_path = Path(scan["start_xyz"]["path"])
    checks["raw_output_hash"] = raw_path.is_file() and sha256(raw_path) == scan.get("raw_output_sha256")
    checks["optimized_xyz_hash"] = xyz_path.is_file() and sha256(xyz_path) == scan["optimized_xyz"].get("sha256")
    checks["start_xyz_hash"] = start_path.is_file() and sha256(start_path) == scan["start_xyz"].get("sha256")
    raw_text = raw_path.read_text(errors="replace") if raw_path.is_file() else ""
    checks["raw_method_evidence"] = all(
        marker in raw_text
        for marker in ("DF-MP2", "6-31G(D)", "RHF Reference", "charge = 0, multiplicity = 1")
    )
    checks["raw_convergence_evidence"] = "Optimization converged!" in raw_text

    rotor = next(item for item in auth["rotatable_heavy_atom_torsions"] if item["rotor_id"] == args.expected_rotor_id)
    ids = [int(value) for value in rotor["one_based_ordinals"]]
    elements, coordinates, _ = read_xyz(xyz_path)
    frame_elements, _, _ = read_xyz(Path(auth["inputs"]["frame342"]["path"]))
    recomputed_angle = dihedral(coordinates, ids)
    recomputed_error = abs(circular_delta(recomputed_angle, float(scan["target_angle_deg"])))
    checks["atom_identity_order"] = elements == frame_elements
    checks["recomputed_dihedral"] = recomputed_error <= 0.5 and abs(circular_delta(recomputed_angle, float(scan["final_dihedral_deg"]))) <= 1e-5

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_torsion_scan_point_independent_validation",
        "status": "pass_independent_mp2_torsion_scan_point" if all(checks.values()) else "fail",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "validator_sha256": sha256(Path(__file__).resolve()),
        "authorization": {"path": str(auth_path), "sha256": sha256(auth_path)},
        "canary_validation": {"path": str(canary_path), "sha256": sha256(canary_path)},
        "scan_report": {"path": str(scan_path), "sha256": sha256(scan_path)},
        "rotor_id": args.expected_rotor_id,
        "signed_step_index": args.expected_step_index,
        "final_energy_hartree": scan.get("final_energy_hartree"),
        "final_dihedral_deg": scan.get("final_dihedral_deg"),
        "recomputed_dihedral_deg": recomputed_angle,
        "recomputed_target_error_deg": recomputed_error,
        "checks": checks,
        "next_adjacent_point_released": all(checks.values()),
        "production_approved": False,
    }
    report_path = output_dir / "O6U_MP2_TORSION_SCAN_POINT_INDEPENDENT_VALIDATION.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
