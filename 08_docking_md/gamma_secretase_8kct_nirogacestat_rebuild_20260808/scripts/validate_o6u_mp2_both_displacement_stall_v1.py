#!/usr/bin/env python3
"""Authorize stopping the live OptKing BOTH canary after a displacement stall.

This validator is restricted to the immutable byte prefix captured by the
independent restart extractor.  It never labels the geometry converged and it
cannot release parameter fitting, system construction, or MD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


FLOAT_RE = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")
MAX_FORCE = 1.5e-5
RMS_FORCE = 1.0e-5
MAX_DISP = 6.0e-5
RMS_DISP = 4.0e-5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def parse_steps(text: str) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in text.splitlines():
        if "~" not in line:
            continue
        values = FLOAT_RE.findall(line)
        if len(values) != 7:
            continue
        step = int(values[0])
        energy, delta_e, max_force, rms_force, max_disp, rms_disp = map(float, values[1:])
        rows.append({
            "step": step,
            "energy": energy,
            "delta_e": delta_e,
            "max_force": max_force,
            "rms_force": rms_force,
            "max_disp": max_disp,
            "rms_disp": rms_disp,
        })
    return rows


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    return ordered[n // 2] if n % 2 else 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--snapshot-report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tail-steps", type=int, default=24)
    args = parser.parse_args()
    record_path = args.record.resolve()
    snapshot_path = args.snapshot_report.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    if args.tail_steps < 20:
        raise SystemExit("Displacement-stall assessment requires at least 20 steps")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if record.get("status") != "running" or record.get("optimizer_strategy_version") != "optking_both_dynamic_level1_recovery_v1":
        raise SystemExit("Input is not the live OptKing BOTH recovery canary")
    if record.get("method") != "DF-MP2/6-31G(d), frozen core, RHF reference" or record.get("optimizer_convergence") != "gau_tight":
        raise SystemExit("Frozen model chemistry or convergence differs")
    pid = record.get("pid")
    if not isinstance(pid, int) or pid <= 0:
        raise SystemExit("Run record has no valid PID")
    try:
        os.kill(pid, 0)
    except OSError as exc:
        raise SystemExit("Recorded live PID is not running") from exc
    if snapshot.get("status") != "pass_technical_restart_candidate_not_converged" or snapshot.get("role") != "snapshot_only":
        raise SystemExit("Snapshot is not an independent running-run snapshot")
    if snapshot.get("record", {}).get("sha256") != sha256(record_path):
        raise SystemExit("Snapshot is not bound to the current live record")

    raw = snapshot.get("raw_output", {})
    raw_path = Path(str(raw.get("path", ""))).resolve()
    byte_count = raw.get("bytes_at_extraction")
    if not isinstance(byte_count, int) or byte_count <= 0 or not raw_path.is_file():
        raise SystemExit("Snapshot raw-output metadata is invalid")
    with raw_path.open("rb") as handle:
        prefix = handle.read(byte_count)
    if len(prefix) != byte_count or sha256_bytes(prefix) != raw.get("sha256_at_extraction"):
        raise SystemExit("Immutable raw-output prefix differs from the snapshot")

    rows = parse_steps(prefix.decode("utf-8", errors="replace"))
    if len(rows) < args.tail_steps:
        raise SystemExit("Too few complete OptKing step rows")
    tail = rows[-args.tail_steps:]
    first = tail[: len(tail) // 2]
    second = tail[len(tail) // 2 :]
    energies = [float(row["energy"]) for row in tail]
    max_disp = [float(row["max_disp"]) for row in tail]
    metrics = {
        "completed_step_rows": len(rows),
        "tail_steps": len(tail),
        "tail_first_step": int(tail[0]["step"]),
        "tail_last_step": int(tail[-1]["step"]),
        "energy_range_hartree": max(energies) - min(energies),
        "all_max_force_pass": all(float(row["max_force"]) <= MAX_FORCE for row in tail),
        "all_rms_force_pass": all(float(row["rms_force"]) <= RMS_FORCE for row in tail),
        "all_rms_displacement_pass": all(float(row["rms_disp"]) <= RMS_DISP for row in tail),
        "maximum_displacement_pass_count": sum(value <= MAX_DISP for value in max_disp),
        "maximum_displacement_min": min(max_disp),
        "maximum_displacement_median_first_half": median([float(row["max_disp"]) for row in first]),
        "maximum_displacement_median_second_half": median([float(row["max_disp"]) for row in second]),
    }
    pass_stall = (
        metrics["energy_range_hartree"] <= 5.0e-7
        and metrics["all_max_force_pass"]
        and metrics["all_rms_force_pass"]
        and metrics["all_rms_displacement_pass"]
        and metrics["maximum_displacement_pass_count"] == 0
        and metrics["maximum_displacement_median_second_half"] >= 0.90 * metrics["maximum_displacement_median_first_half"]
    )
    if not pass_stall:
        raise SystemExit(f"Frozen displacement-stall criteria not satisfied: {metrics}")

    output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_both_displacement_stall_termination_authorization",
        "status": "pass_displacement_stall_stop_authorized_no_convergence",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "pid": pid,
        "metrics": metrics,
        "thresholds": {"max_force": MAX_FORCE, "rms_force": RMS_FORCE, "max_disp": MAX_DISP, "rms_disp": RMS_DISP},
        "run_record": {"path": str(record_path), "sha256": sha256(record_path)},
        "snapshot_report": {"path": str(snapshot_path), "sha256": sha256(snapshot_path)},
        "authorized_action": "Stop only the recorded stalled PID, close and hash its raw output, then restart from the identity-validated snapshot using Psi4's geomeTRIC engine with TRIC coordinates and unchanged model chemistry and gau_tight convergence.",
        "references": [
            "https://psi4.github.io/psi4docs/master/optking.html",
            "https://psi4.github.io/psi4docs/master/opt.html",
            "https://doi.org/10.1063/1.4952956",
        ],
        "release_boundary": "Technical optimizer change only; no convergence, parameter fitting, CHARMM-GUI, or MD is approved.",
    }
    out = output_dir / "O6U_MP2_BOTH_DISPLACEMENT_STALL_TERMINATION_AUTHORIZATION.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report": str(out), "sha256": sha256(out), "metrics": metrics}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
