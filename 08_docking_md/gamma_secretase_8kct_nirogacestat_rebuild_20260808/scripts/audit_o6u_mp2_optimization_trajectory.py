#!/usr/bin/env python3
"""Create a read-only technical snapshot of an active O6U MP2 optimization."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from validate_o6u_mp2_optimization_canary import GAU_TIGHT_LIMITS, parse_convergence_rows


FAILURE_MARKERS = (
    "Optimization failed",
    "Back transformation failed",
    "PsiException",
    "Traceback",
)
COMPLETION_MARKER = "Optimization is complete"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def best(rows: list[dict[str, float | int]], key: str) -> dict[str, float | int]:
    row = min(rows, key=lambda item: abs(float(item[key])))
    return {"step": int(row["step"]), "value": float(row[key])}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--raw-output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    record_path = args.record.resolve()
    raw_path = args.raw_output.resolve()
    report_path = args.report.resolve()
    for path in (record_path, raw_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required input: {path}")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    raw = raw_path.read_text(encoding="utf-8", errors="replace")
    rows = parse_convergence_rows(raw)
    if not rows:
        raise SystemExit("No OptKing convergence rows were reconstructed")
    steps = [int(row["step"]) for row in rows]
    if steps != list(range(1, steps[-1] + 1)):
        raise SystemExit(f"OptKing step sequence is incomplete or duplicated: {steps}")

    failure_hits = {marker: raw.count(marker) for marker in FAILURE_MARKERS if marker in raw}
    completion_count = raw.count(COMPLETION_MARKER)
    last = rows[-1]
    ratios = {
        "max_force": float(last["max_force_au"]) / GAU_TIGHT_LIMITS["max_force"],
        "rms_force": float(last["rms_force_au"]) / GAU_TIGHT_LIMITS["rms_force"],
        "max_displacement": float(last["max_displacement_au"]) / GAU_TIGHT_LIMITS["max_displacement"],
        "rms_displacement": float(last["rms_displacement_au"]) / GAU_TIGHT_LIMITS["rms_displacement"],
    }
    uphill_steps = [int(row["step"]) for row in rows if float(row["delta_energy_hartree"]) > 0]
    lowest_energy_row = min(rows, key=lambda row: float(row["energy_hartree"]))
    maximum_consecutive_uphill = 0
    current_uphill = 0
    for row in rows:
        if float(row["delta_energy_hartree"]) > 0:
            current_uphill += 1
            maximum_consecutive_uphill = max(maximum_consecutive_uphill, current_uphill)
        else:
            current_uphill = 0

    if failure_hits:
        status = "technical_failure_marker_present"
    elif completion_count:
        status = "completion_marker_present_requires_independent_final_validation"
    elif record.get("status") == "running":
        status = "healthy_running_unconverged"
    else:
        status = "unexpected_record_state"

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_optimization_trajectory_snapshot",
        "status": status,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_record": {"path": str(record_path), "sha256": sha256(record_path), "status": record.get("status")},
        "raw_output": {
            "path": str(raw_path),
            "sha256_at_snapshot": sha256(raw_path),
            "size_bytes_at_snapshot": raw_path.stat().st_size,
            "mtime_ns_at_snapshot": raw_path.stat().st_mtime_ns,
        },
        "step_count": len(rows),
        "last_step": last,
        "gau_tight_limits": GAU_TIGHT_LIMITS,
        "last_to_limit_ratios": ratios,
        "best_observed": {
            "energy": {
                "step": int(lowest_energy_row["step"]),
                "value": float(lowest_energy_row["energy_hartree"]),
            },
            "max_force": best(rows, "max_force_au"),
            "rms_force": best(rows, "rms_force_au"),
            "max_displacement": best(rows, "max_displacement_au"),
            "rms_displacement": best(rows, "rms_displacement_au"),
        },
        "uphill_energy_steps": uphill_steps,
        "maximum_consecutive_uphill_steps": maximum_consecutive_uphill,
        "failure_markers": failure_hits,
        "completion_marker_count": completion_count,
        "interpretation": (
            "This is a read-only technical trajectory snapshot. It does not release any QM target, "
            "force-field fitting, CHARMM-GUI construction, or MD stage."
        ),
        "rows": rows,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = report_path.with_suffix(report_path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, report_path)
    print(json.dumps({"status": status, "step_count": len(rows), "report_sha256": sha256(report_path)}, sort_keys=True))
    return 1 if failure_hits or status == "unexpected_record_state" else 0


if __name__ == "__main__":
    raise SystemExit(main())
