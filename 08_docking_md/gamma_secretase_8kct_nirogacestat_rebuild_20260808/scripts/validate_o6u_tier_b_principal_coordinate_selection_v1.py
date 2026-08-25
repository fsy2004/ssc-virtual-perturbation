#!/usr/bin/env python3
"""Independent structural validation of the frozen Tier-B selection report."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    selection_path = Path(args.selection).resolve()
    selection = json.loads(selection_path.read_text())
    if selection["status"] != "pass_tier_b_principal_coordinate_selected":
        raise ValueError("Selection status is not passed")
    decision = selection["decision"]
    evidence = selection["evidence"]
    if decision["principal_coordinate"] != "ROT_C09_N04":
        raise ValueError("Unexpected principal coordinate")
    if decision["holdout_coordinate"] != "ROT_C14_N04":
        raise ValueError("Unexpected holdout coordinate")
    if decision["scan_points"] != [0, 15, -15, 30, -30]:
        raise ValueError("Tier-B scan scope changed")
    if evidence["ROT_C09_N04_observed_circular_spread_deg"] <= evidence["ROT_C14_N04_observed_circular_spread_deg"]:
        raise ValueError("Reported conformational evidence does not support selection")
    if evidence["ROT_C09_N04_plus15_initial_cgenff_absolute_delta_error_kcal_mol"] <= 0:
        raise ValueError("Prescreen mismatch evidence is absent")
    expected_inputs = {
        "authorization": "b87cf8117aa4abf37c64ca2248fa85e3d4a1e7100338535b8b3611d79c60d2d7",
        "adaptive_scope": "519d587648c6458e23dca0e7cc9041a3a711feb6d793b55bd3efcdc532aaa49e",
        "prescreen": "db016b95bb2ddffea96d5ffdc6e83ce90ebdf0ad46861aad342deb1c130a41be",
        "c09_plus15_validation": "85b63aac47de3c9aa3f74dd2965492070fc14f5930192364f590b83a6768da21",
    }
    for name, expected in expected_inputs.items():
        item = selection["inputs"][name]
        path = Path(item["path"])
        if item["sha256"] != expected or sha(path) != expected:
            raise ValueError(f"Independent input verification failed: {name}")
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_tier_b_principal_coordinate_selection_independent_validation",
        "status": "pass_independent_tier_b_principal_coordinate_selection",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection": {"path": str(selection_path), "sha256": sha(selection_path)},
        "validated_principal_coordinate": "ROT_C09_N04",
        "validated_holdout_coordinate": "ROT_C14_N04",
        "validated_scan_points": [0, 15, -15, 30, -30],
        "production_approved": False,
    }
    output = output_dir / "O6U_TIER_B_PRINCIPAL_COORDINATE_SELECTION_INDEPENDENT_VALIDATION.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "report": str(output), "sha256": sha(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
