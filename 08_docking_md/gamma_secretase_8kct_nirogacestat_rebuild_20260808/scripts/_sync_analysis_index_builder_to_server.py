#!/usr/bin/env python3
"""Atomically upload and verify the formal analysis-index builder only."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
LOCAL = ROOT / "scripts" / "build_analysis_ndx.py"
REMOTE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/scripts/build_analysis_ndx.py"


def hash_handle(handle) -> str:
    digest = hashlib.sha256()
    while chunk := handle.read(1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    expected = hashlib.sha256(LOCAL.read_bytes()).hexdigest()
    temporary = REMOTE + ".uploading"
    client = connect()
    sftp = client.open_sftp()
    try:
        sftp.put(str(LOCAL), temporary)
        with sftp.open(temporary, "rb") as handle:
            if hash_handle(handle) != expected:
                raise RuntimeError("staged analysis-index builder hash mismatch")
        sftp.posix_rename(temporary, REMOTE)
        with sftp.open(REMOTE, "rb") as handle:
            observed = hash_handle(handle)
        if observed != expected:
            raise RuntimeError("final analysis-index builder hash mismatch")
    finally:
        try:
            sftp.remove(temporary)
        except OSError:
            pass
        sftp.close()
        client.close()
    print(f"REMOTE ANALYSIS INDEX BUILDER VERIFIED sha256={expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
