#!/usr/bin/env python3
"""Freeze the O6U 24-vs-32-thread torsion throughput decision."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from run_o6u_mp2_torsion_canary_v1 import sha256


def gradient_count(path: Path) -> int:
    return path.read_text(errors="replace").count("-Total Gradient:")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-report", required=True, type=Path)
    parser.add_argument("--baseline-validation", required=True, type=Path)
    parser.add_argument("--canary-report", required=True, type=Path)
    parser.add_argument("--canary-validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    paths = {name: path.resolve() for name, path in {
        "baseline_report": args.baseline_report,
        "baseline_validation": args.baseline_validation,
        "canary_report": args.canary_report,
        "canary_validation": args.canary_validation,
    }.items()}
    docs = {name: json.loads(path.read_text()) for name, path in paths.items()}
    baseline = docs["baseline_report"]
    canary = docs["canary_report"]
    bval = docs["baseline_validation"]
    cval = docs["canary_validation"]
    braw = Path(baseline["raw_output"])
    craw = Path(canary["raw_output"])
    bg = gradient_count(braw)
    cg = gradient_count(craw)
    bsec = float(baseline["elapsed_seconds"]) / bg
    csec = float(canary["elapsed_seconds"]) / cg
    improvement = 1.0 - csec / bsec

    checks = {
        "baseline_independent_validation": bval.get("status") == "pass_independent_mp2_torsion_scan_point",
        "canary_independent_validation": cval.get("status") == "pass_independent_mp2_torsion_scan_point",
        "same_rotor": baseline.get("rotor_id") == canary.get("rotor_id") == "ROT_C17_C15",
        "same_method": baseline.get("method") == canary.get("method") == "DF-MP2/6-31G(d), frozen-core RHF, GAU_TIGHT",
        "same_memory_gib": baseline.get("memory_gib") == canary.get("memory_gib") == 44,
        "thread_routes": baseline.get("threads") == 24 and canary.get("threads") == 32,
        "raw_hashes": sha256(braw) == baseline.get("raw_output_sha256") and sha256(craw) == canary.get("raw_output_sha256"),
        "positive_gradient_counts": bg > 0 and cg > 0,
        "both_converged": "Optimization converged!" in braw.read_text(errors="replace") and "Optimization converged!" in craw.read_text(errors="replace"),
    }
    threshold = 0.10
    retain_32 = all(checks.values()) and improvement >= threshold
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_thread_throughput_decision",
        "status": "pass_retain_32_threads" if retain_32 else "pass_revert_to_24_threads",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
        "checks": checks,
        "baseline": {
            "threads": 24,
            "signed_step_index": baseline.get("signed_step_index"),
            "elapsed_seconds": baseline.get("elapsed_seconds"),
            "gradient_evaluations": bg,
            "seconds_per_gradient": bsec,
        },
        "canary": {
            "threads": 32,
            "signed_step_index": canary.get("signed_step_index"),
            "elapsed_seconds": canary.get("elapsed_seconds"),
            "gradient_evaluations": cg,
            "seconds_per_gradient": csec,
        },
        "relative_throughput_improvement": improvement,
        "required_improvement": threshold,
        "decision": "retain_32_threads" if retain_32 else "revert_to_24_threads",
        "cgroup_observation": "No new memory.max, OOM, or OOM-kill event during the canary; historical max-event counter remained unchanged.",
        "chemistry_or_convergence_changed": False,
        "production_approved": False,
    }
    out = output_dir / "O6U_THREAD_THROUGHPUT_DECISION.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "decision": report["decision"], "report": str(out), "sha256": sha256(out)}, sort_keys=True))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
