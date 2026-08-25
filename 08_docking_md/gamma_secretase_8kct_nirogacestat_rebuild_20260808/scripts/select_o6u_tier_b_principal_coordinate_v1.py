#!/usr/bin/env python3
"""Select the compound-specific Tier-B principal coordinate from frozen evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rec(path: Path) -> dict[str, object]:
    path = path.resolve()
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha(path)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--authorization", required=True)
    parser.add_argument("--adaptive-scope", required=True)
    parser.add_argument("--prescreen", required=True)
    parser.add_argument("--c09-plus15-validation", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    paths = {name: Path(value).resolve() for name, value in vars(args).items() if name != "output_dir"}
    expected = {
        "authorization": "b87cf8117aa4abf37c64ca2248fa85e3d4a1e7100338535b8b3611d79c60d2d7",
        "adaptive_scope": "519d587648c6458e23dca0e7cc9041a3a711feb6d793b55bd3efcdc532aaa49e",
        "prescreen": "db016b95bb2ddffea96d5ffdc6e83ce90ebdf0ad46861aad342deb1c130a41be",
        "c09_plus15_validation": "85b63aac47de3c9aa3f74dd2965492070fc14f5930192364f590b83a6768da21",
    }
    for name, path in paths.items():
        if sha(path) != expected[name]:
            raise ValueError(f"Frozen input hash mismatch: {name}")
    authorization = json.loads(paths["authorization"].read_text())
    scope = json.loads(paths["adaptive_scope"].read_text())
    prescreen = json.loads(paths["prescreen"].read_text())
    validation = json.loads(paths["c09_plus15_validation"].read_text())
    if authorization["status"] != "pass_bonded_torsion_scope_authorized_canary_only":
        raise ValueError("Authorization gate is not passed")
    if scope["status"] != "pass_compound_specific_adaptive_scan_scope_authorized":
        raise ValueError("Adaptive-scope gate is not passed")
    if prescreen["status"] != "pass_tier_b_initial_cgenff_mismatch_prescreen":
        raise ValueError("Initial-CGenFF prescreen is not passed")
    if validation["status"] != "pass_independent_mp2_torsion_scan_point":
        raise ValueError("C09 +15 independent validation is not passed")
    rows = {
        row["rotor_id"]: row
        for row in authorization["rotatable_heavy_atom_torsions"]
        if row["rotor_id"] in {"ROT_C09_N04", "ROT_C14_N04"}
    }
    if set(rows) != {"ROT_C09_N04", "ROT_C14_N04"}:
        raise ValueError("Tier-B rotor definitions are incomplete")
    c09_spread = float(rows["ROT_C09_N04"]["observed_circular_spread_deg_across_frames1_342_641_768_and_native8KCT"])
    c14_spread = float(rows["ROT_C14_N04"]["observed_circular_spread_deg_across_frames1_342_641_768_and_native8KCT"])
    mismatch = float(prescreen["comparison"]["absolute_mm_minus_qm_delta_kcal_mol"])
    if not all(math.isfinite(x) for x in (c09_spread, c14_spread, mismatch)):
        raise ValueError("Non-finite Tier-B selection evidence")
    if c09_spread <= c14_spread:
        raise ValueError("Frozen conformer/native coverage does not favor C09")
    if prescreen["scope"]["rotor_id"] != "ROT_C09_N04":
        raise ValueError("Prescreen rotor identity mismatch")

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_tier_b_principal_coordinate_selection",
        "status": "pass_tier_b_principal_coordinate_selected",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "compound": "nirogacestat/PF-03084014/BRD-K61691541",
        "residue": "O6U",
        "inputs": {name: rec(path) for name, path in paths.items()},
        "evidence": {
            "ROT_C09_N04_observed_circular_spread_deg": c09_spread,
            "ROT_C14_N04_observed_circular_spread_deg": c14_spread,
            "ROT_C09_N04_plus15_initial_cgenff_absolute_delta_error_kcal_mol": mismatch,
            "ROT_C09_N04_plus15_qm_delta_kcal_mol": prescreen["comparison"]["qm_delta_kcal_mol"],
            "ROT_C09_N04_plus15_initial_cgenff_mm_delta_kcal_mol": prescreen["comparison"]["initial_cgenff_mm_delta_kcal_mol"],
        },
        "decision": {
            "principal_coordinate": "ROT_C09_N04",
            "scan_points": [0, 15, -15, 30, -30],
            "already_completed_and_independently_validated": [15],
            "holdout_coordinate": "ROT_C14_N04",
            "holdout_policy": "geometry/energy holdout unless a compact coupled +/-15 cross-check is triggered by persistent fitted-model mismatch",
            "reason": "C09 spans more of the retained/native conformational coverage and its validated +15 point shows a measurable same-geometry initial-CGenFF energy mismatch.",
        },
        "methodology_anchors": [
            "10.1002/jcc.21367",
            "10.1021/acs.jctc.5c00046",
            "10.1002/jcc.23422",
            "10.1063/5.0196848",
        ],
        "disease_and_compound_anchors": [
            "10.1136/ard.2010.134742",
            "10.1002/art.30254",
            "10.1158/1535-7163.MCT-10-0034",
            "10.1038/s41594-024-01439-8",
        ],
        "release_boundary": "Scope selection only. No parameter release, affinity/efficacy inference, or production-MD approval.",
    }
    output = output_dir / "O6U_TIER_B_PRINCIPAL_COORDINATE_SELECTION.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "report": str(output), "sha256": sha(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
