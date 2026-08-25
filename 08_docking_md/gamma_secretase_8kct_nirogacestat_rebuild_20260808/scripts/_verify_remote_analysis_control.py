#!/usr/bin/env python3
"""Verify every synchronized analysis-control file by remote SHA-256."""
from __future__ import annotations

import hashlib
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts import _sync_analysis_to_server as sync  # noqa: E402
from scripts.new_md_server import connect, run  # noqa: E402


def local_records() -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for name in sync.UPLOAD_SCRIPTS:
        local = sync.ROOT / "scripts" / name
        if local.is_file():
            records.append((f"{sync.REMOTE_BASE}/scripts/{name}", hashlib.sha256(local.read_bytes()).hexdigest()))
    for name in sync.UPLOAD_CONFIG:
        local = sync.ROOT / "config" / name
        if local.is_file():
            records.append((f"{sync.REMOTE_BASE}/config/{name}", hashlib.sha256(local.read_bytes()).hexdigest()))
    ndx = sync.ROOT / "builds" / "analysis.ndx"
    records.append((f"{sync.REMOTE_BASE}/builds/analysis.ndx", hashlib.sha256(ndx.read_bytes()).hexdigest()))
    return records


def main() -> int:
    expected = dict(local_records())
    command = "sha256sum -- " + " ".join(shlex.quote(path) for path in expected)
    connection = connect()
    try:
        code, out, err = run(connection, command, timeout=120)
    finally:
        connection.close()
    if code != 0:
        print(err, file=sys.stderr, end="")
        return code or 1
    observed: dict[str, str] = {}
    for line in out.splitlines():
        digest, path = line.split(maxsplit=1)
        observed[path] = digest
    mismatches = []
    for path, digest in expected.items():
        if observed.get(path) != digest:
            mismatches.append((path, digest, observed.get(path, "MISSING")))
    if mismatches:
        for path, wanted, got in mismatches:
            print(f"MISMATCH {path}: expected={wanted} observed={got}")
        return 1
    print(f"REMOTE CONTROL VERIFY PASS: {len(expected)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
