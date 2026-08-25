#!/usr/bin/env python3
"""Freeze a CPU/RAM/disk-aware MPI plan after the technical canary."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


GIB = 1024 ** 3
MODELS = ["PB_membrane_indi4"]


def build_resource_plan(
    logical_cpus: int,
    mem_total_bytes: int,
    disk_free_bytes: int,
    canary_peak_rss_bytes: dict[str, int],
    projected_disk_bytes: int,
) -> dict[str, Any]:
    if logical_cpus < 5:
        raise ValueError("logical_cpus: at least five are required")
    if mem_total_bytes <= 0:
        raise ValueError("mem_total_bytes must be positive")
    if disk_free_bytes <= 0 or projected_disk_bytes < 0:
        raise ValueError("disk byte counts are invalid")

    cpu_reserved = 2
    cpu_usable = logical_cpus - cpu_reserved
    memory_reserved = max(16 * GIB, int(mem_total_bytes * 0.05))
    memory_usable = mem_total_bytes - memory_reserved
    disk_reserved = max(50 * GIB, int(disk_free_bytes * 0.10))
    if disk_free_bytes < projected_disk_bytes + disk_reserved:
        raise ValueError(
            "disk: projected derived inputs and outputs would violate the safety reserve"
        )

    missing = [model for model in MODELS if model not in canary_peak_rss_bytes]
    if missing:
        raise ValueError(f"missing canary peak RSS for {', '.join(missing)}")

    model_plans: dict[str, Any] = {}
    for model in MODELS:
        peak = canary_peak_rss_bytes[model]
        if not isinstance(peak, int) or peak <= 0:
            raise ValueError(f"{model}: canary peak RSS must be a positive integer")
        memory_rank_cap = memory_usable // peak
        total_rank_cap = min(cpu_usable, memory_rank_cap, 300 * 3)
        if total_rank_cap < 3:
            raise ValueError(
                f"{model}: resources cannot provide at least one MPI rank per replica"
            )
        concurrent_jobs = 3
        ranks_per_job = total_rank_cap // concurrent_jobs
        total_ranks = ranks_per_job * concurrent_jobs
        cpu_sets = {}
        for job_index, replica in enumerate(("rep01", "rep02", "rep03")):
            first = cpu_reserved + job_index * ranks_per_job
            last = first + ranks_per_job - 1
            cpu_sets[replica] = str(first) if first == last else f"{first}-{last}"
        model_plans[model] = {
            "concurrent_jobs": concurrent_jobs,
            "mpi_ranks_per_job": ranks_per_job,
            "total_mpi_ranks": total_ranks,
            "canary_peak_rss_bytes_per_rank": peak,
            "estimated_peak_rss_bytes": total_ranks * peak,
            "limiting_resource": "memory" if memory_rank_cap < cpu_usable else "cpu",
            "replica_batch": ["rep01", "rep02", "rep03"],
            "cpu_sets": cpu_sets,
        }

    return {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_resource_plan",
        "status": "frozen_before_formal_endpoint_energy_results",
        "cpu": {
            "logical": logical_cpus,
            "reserved": cpu_reserved,
            "usable": cpu_usable,
        },
        "memory": {
            "total_bytes": mem_total_bytes,
            "reserved_bytes": memory_reserved,
            "usable_bytes": memory_usable,
        },
        "disk": {
            "free_bytes": disk_free_bytes,
            "reserved_bytes": disk_reserved,
            "projected_bytes": projected_disk_bytes,
        },
        "environment": {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
        "gpu_required": False,
        "models": model_plans,
    }


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--canary-report", type=Path, required=True)
    parser.add_argument("--projected-disk-gib", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix(args.output.suffix + ".sha256").exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    canary = json.loads(args.canary_report.read_text(encoding="utf-8"))
    plan = build_resource_plan(
        logical_cpus=int(inventory["logical_cpus"]),
        mem_total_bytes=int(inventory["mem_total_bytes"]),
        disk_free_bytes=int(inventory["disk_free_bytes"]),
        canary_peak_rss_bytes={
            model: int(canary["models"][model]["peak_rss_bytes_per_rank"])
            for model in MODELS
        },
        projected_disk_bytes=int(args.projected_disk_gib * GIB),
    )
    rendered = json.dumps(plan, indent=2, sort_keys=True) + "\n"
    args.output.write_text(rendered, encoding="utf-8")
    digest = file_sha256(args.output)
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest}  {args.output.name}\n", encoding="ascii"
    )
    print(json.dumps({"status": "pass", "output": str(args.output), "sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
