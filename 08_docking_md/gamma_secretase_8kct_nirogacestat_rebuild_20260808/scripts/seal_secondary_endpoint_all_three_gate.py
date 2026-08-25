#!/usr/bin/env python3
"""Seal the all-three post-MD eligibility gate for endpoint preprocessing."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPLICAS = ("rep01", "rep02", "rep03")
PRODUCTION_RELEASE_SHA256 = "a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_completion(replica: str, payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("status") != "pass" or payload.get("report_type") != "production_500ns_completion":
        raise ValueError(f"{replica}: production 500 ns completion status is not passing")
    if payload.get("replica") != replica:
        raise ValueError(f"{replica}: completion report replica mismatch")
    if int(payload.get("final_step", -1)) != 125000000 or float(payload.get("final_time_ps", math.nan)) != 500000.0:
        raise ValueError(f"{replica}: exact 500 ns endpoint is not present")
    if payload.get("production_release_sha256") != PRODUCTION_RELEASE_SHA256:
        raise ValueError(f"{replica}: production release hash mismatch")
    return {"status": "pass"}


def validate_analysis_evidence(category: str, payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("status") == "pass":
        if category == "pbc" and float(payload.get("maximum_absolute_difference_nm", math.inf)) > 0.01:
            raise ValueError("pbc: distance invariance exceeds 0.01 nm")
        return {"status": "pass"}
    if category == "membrane" and all(
        payload.get(field) == "pass"
        for field in ("technical_status", "sampling_status", "preproduction_status")
    ):
        return {"status": "pass"}
    raise ValueError(f"{category}: required analysis evidence is not passing")


def _parse_mapping(values: list[str], label: str) -> dict[str, Path]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{label}: invalid mapping")
        replica, raw = value.split("=", 1)
        if replica not in REPLICAS or replica in result:
            raise ValueError(f"{label}: unknown or duplicate replica {replica}")
        result[replica] = Path(raw)
    missing = [replica for replica in REPLICAS if replica not in result]
    if missing:
        raise ValueError(f"{label}: missing {', '.join(missing)}")
    return result


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def seal(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    mappings = {
        "completion": _parse_mapping(args.completion, "completion"),
        "trajectory_integrity": _parse_mapping(args.trajectory_integrity, "trajectory_integrity"),
        "pbc": _parse_mapping(args.pbc, "pbc"),
        "membrane": _parse_mapping(args.membrane, "membrane"),
        "energy": _parse_mapping(args.energy, "energy"),
    }
    runner = _load(args.runner_status)
    if runner.get("production_runners_active") not in (False, 0, []):
        raise ValueError("one or more production runners remain active")
    evidence = {}
    for replica in REPLICAS:
        evidence[replica] = {}
        completion_path = mappings["completion"][replica]
        validate_completion(replica, _load(completion_path))
        evidence[replica]["completion"] = {"path": str(completion_path), "sha256": file_sha256(completion_path)}
        for category in ("trajectory_integrity", "pbc", "membrane", "energy"):
            path = mappings[category][replica]
            validate_analysis_evidence(category, _load(path))
            evidence[replica][category] = {"path": str(path), "sha256": file_sha256(path)}
    report = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_all_three_gate",
        "status": "pass",
        "eligible_replicas": list(REPLICAS),
        "all_three_500ns_complete": True,
        "all_required_gates_passed": True,
        "production_runners_active": False,
        "production_release_sha256": PRODUCTION_RELEASE_SHA256,
        "runner_status": {"path": str(args.runner_status), "sha256": file_sha256(args.runner_status)},
        "evidence": evidence,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(args.output)
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(f"{digest}  {args.output.name}\n", encoding="ascii")
    return {"status": "pass", "output": str(args.output), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    for option in ("completion", "trajectory-integrity", "pbc", "membrane", "energy"):
        parser.add_argument(f"--{option}", action="append", required=True, metavar="REP=REPORT.json")
    parser.add_argument("--runner-status", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(seal(args), sort_keys=True))


if __name__ == "__main__":
    main()
