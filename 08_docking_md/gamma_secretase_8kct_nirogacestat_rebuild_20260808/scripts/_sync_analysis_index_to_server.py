#!/usr/bin/env python3
"""Atomically replace only the remote analysis index and verify its SHA-256."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
LOCAL = ROOT / "builds" / "analysis.ndx"
REMOTE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/builds/analysis.ndx"


def digest_stream(handle) -> str:
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
            staged = digest_stream(handle)
        if staged != expected:
            raise RuntimeError(f"staged hash mismatch: expected={expected} observed={staged}")
        sftp.posix_rename(temporary, REMOTE)
        with sftp.open(REMOTE, "rb") as handle:
            observed = digest_stream(handle)
        if observed != expected:
            raise RuntimeError(f"final hash mismatch: expected={expected} observed={observed}")
    finally:
        try:
            sftp.remove(temporary)
        except OSError:
            pass
        sftp.close()
        client.close()
    print(f"REMOTE ANALYSIS INDEX VERIFIED sha256={expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
