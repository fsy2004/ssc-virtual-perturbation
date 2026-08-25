#!/usr/bin/env python3
"""Authorize termination of a live O6U MP2 optimization only after a hard plateau.

This validator reads only the immutable byte prefix recorded by the independent
snapshot extractor.  It does not claim convergence and cannot release force-field
fitting or MD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path


ENERGY_RE = re.compile(r"^\s*Total Energy\s+=\s+(-?\d+\.\d+)\s+\[Eh\]", re.M)
GRADIENT_HEADER = "-Total Gradient:"
EXPECTED_ATOMS = 76


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def gradient_norms(text: str) -> list[tuple[float, float]]:
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if GRADIENT_HEADER in line]
    norms: list[tuple[float, float]] = []
    for start in starts:
        values: list[float] = []
        for line in lines[start + 3 : start + 3 + EXPECTED_ATOMS]:
            fields = line.split()
            if len(fields) < 4:
                break
            try:
                values.extend(float(value) for value in fields[-3:])
            except ValueError:
                break
        if len(values) == EXPECTED_ATOMS * 3:
            norms.append((max(abs(value) for value in values), math.sqrt(sum(value * value for value in values) / len(values))))
    return norms


def relative_span(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return (max(values) - min(values)) / abs(mean) if mean else math.inf


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--snapshot-report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--plateau-steps", type=int, default=12)
    args = parser.parse_args()
    record_path = args.record.resolve()
    snapshot_path = args.snapshot_report.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    if args.plateau_steps < 8:
        raise SystemExit("Plateau assessment requires at least 8 completed gradients")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if record.get("status") != "running" or record.get("optimizer_strategy_version") != "cartesian_dynamic_level4_recovery_v2":
        raise SystemExit("Input is not the live Cartesian dynamic-level-4 recovery")
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
        raise SystemExit("Snapshot report is not an independent running-run snapshot")
    if snapshot.get("record", {}).get("sha256") != sha256(record_path):
        raise SystemExit("Snapshot report is not bound to the current live record")
    raw = snapshot.get("raw_output", {})
    raw_path = Path(str(raw.get("path", ""))).resolve()
    byte_count = raw.get("bytes_at_extraction")
    if not isinstance(byte_count, int) or byte_count <= 0 or not raw_path.is_file():
        raise SystemExit("Snapshot byte-prefix metadata is invalid")
    with raw_path.open("rb") as handle:
        prefix = handle.read(byte_count)
    if len(prefix) != byte_count or sha256_bytes(prefix) != raw.get("sha256_at_extraction"):
        raise SystemExit("Immutable raw-output prefix differs from the snapshot report")

    text = prefix.decode("utf-8", errors="replace")
    energies = [float(value) for value in ENERGY_RE.findall(text)]
    gradients = gradient_norms(text)
    count = args.plateau_steps
    if len(energies) < count or len(gradients) < count:
        raise SystemExit("Too few complete energy/gradient points for plateau assessment")
    e_tail = energies[-count:]
    max_tail = [value[0] for value in gradients[-count:]]
    rms_tail = [value[1] for value in gradients[-count:]]
    metrics = {
        "completed_energy_points": len(energies),
        "completed_gradient_points": len(gradients),
        "plateau_steps": count,
        "energy_range_hartree": max(e_tail) - min(e_tail),
        "max_gradient_last_hartree_per_bohr": max_tail[-1],
        "rms_gradient_last_hartree_per_bohr": rms_tail[-1],
        "max_gradient_relative_span": relative_span(max_tail),
        "rms_gradient_relative_span": relative_span(rms_tail),
    }
    pass_plateau = (
        metrics["energy_range_hartree"] <= 1.0e-9
        and metrics["max_gradient_relative_span"] <= 5.0e-3
        and metrics["rms_gradient_relative_span"] <= 5.0e-3
        and metrics["max_gradient_last_hartree_per_bohr"] > 1.5e-5
        and metrics["rms_gradient_last_hartree_per_bohr"] > 1.0e-5
    )
    if not pass_plateau:
        raise SystemExit(f"Prospective plateau criteria not satisfied: {metrics}")

    output_dir.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_cartesian_plateau_termination_authorization",
        "status": "pass_plateau_stop_authorized_no_convergence",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "pid": pid,
        "metrics": metrics,
        "run_record": {"path": str(record_path), "sha256": sha256(record_path)},
        "snapshot_report": {"path": str(snapshot_path), "sha256": sha256(snapshot_path)},
        "frozen_invariants": {
            "method": record["method"],
            "convergence": record["optimizer_convergence"],
            "charge_e": record["charge_e"],
            "multiplicity": record["multiplicity"],
        },
        "authorized_action": "Interrupt only the recorded plateaued PID and restart from the independently extracted geometry with OptKing BOTH coordinates and dynamic_level 1.",
        "references": [
            "https://psi4.github.io/psi4docs/master/optking.html",
            "https://doi.org/10.1063/1.4952956",
        ],
        "release_boundary": "Technical plateau termination only; no convergence, parameter fitting, CHARMM-GUI, or MD is approved.",
    }
    out = output_dir / "O6U_MP2_CARTESIAN_PLATEAU_TERMINATION_AUTHORIZATION.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report": str(out), "sha256": sha256(out), "metrics": metrics}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

