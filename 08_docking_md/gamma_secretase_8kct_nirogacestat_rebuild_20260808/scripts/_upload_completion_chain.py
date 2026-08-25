#!/usr/bin/env python3
"""Upload the rep01-completion chain script and verify SHA-256."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
REMOTE_BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
FILES = [("scripts/_remote_launch_rep01_completion_chain.sh", "scripts/_remote_launch_rep01_completion_chain.sh")]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_remote(sftp, path: str) -> str:
    digest = hashlib.sha256()
    with sftp.open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    c = connect()
    sftp = c.open_sftp()
    try:
        mismatches = []
        for local_rel, remote_rel in FILES:
            local = ROOT / local_rel
            if not local.is_file():
                print(f"SKIP missing {local_rel}")
                continue
            remote = f"{REMOTE_BASE}/{remote_rel}"
            sftp.put(str(local), remote)
            observed = sha256_remote(sftp, remote)
            expected = sha256_file(local)
            if observed != expected:
                mismatches.append((remote, expected, observed))
            print(f"UP {remote_rel} sha256={expected[:16]}...")
        if mismatches:
            for remote, expected, observed in mismatches:
                print(f"HASH MISMATCH {remote}: expected={expected} observed={observed}")
            sftp.close(); c.close(); return 1
        sftp.close(); c.close()
        print("Uploaded and SHA-256 verified.")
        return 0
    finally:
        try:
            sftp.close()
        except Exception:
            pass
        c.close()


if __name__ == "__main__":
    sys.exit(main())
