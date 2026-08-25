#!/usr/bin/env python3
"""Capture the CPU/RAM/disk inventory used to freeze formal MPI resources."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path


def parse_mem_total_bytes(text: str) -> int:
    match = re.search(r"^MemTotal:\s*(\d+)\s+kB\s*$", text, re.MULTILINE)
    if match is None or int(match.group(1)) <= 0:
        raise ValueError("MemTotal is missing or invalid")
    return int(match.group(1)) * 1024


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect(data_path: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    logical = os.cpu_count()
    if logical is None or logical < 5:
        raise ValueError("at least five logical CPUs are required")
    meminfo = Path("/proc/meminfo").read_text(encoding="ascii")
    disk = shutil.disk_usage(data_path)
    lscpu = subprocess.run(
        ["lscpu"], check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=30,
    ).stdout
    payload = {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_cpu_inventory",
        "status": "pass",
        "logical_cpus": logical,
        "mem_total_bytes": parse_mem_total_bytes(meminfo),
        "disk_path": str(data_path.resolve()),
        "disk_total_bytes": disk.total,
        "disk_used_bytes": disk.used,
        "disk_free_bytes": disk.free,
        "platform": platform.platform(),
        "lscpu": lscpu.splitlines(),
        "gpu_required": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(output)
    output.with_suffix(output.suffix + ".sha256").write_text(f"{digest}  {output.name}\n", encoding="ascii")
    return {"status": "pass", "output": str(output), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-path", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(collect(args.data_path, args.output), sort_keys=True))


if __name__ == "__main__":
    main()
