#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
REMOTE_BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
NAME = "MMGBSA_PBSA_SECONDARY_ANALYSIS_SENSITIVITY_WITHDRAWAL_20260822.md"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    local = ROOT / NAME
    if not local.is_file():
        raise SystemExit(f"missing local secondary amendment: {local}")
    expected = sha256_bytes(local.read_bytes())
    client = connect()
    sftp = client.open_sftp()
    try:
        remote = f"{REMOTE_BASE}/{NAME}"
        sftp.put(str(local), remote)
        with sftp.open(remote, "rb") as handle:
            observed = sha256_bytes(handle.read())
        if observed != expected:
            raise SystemExit(f"hash mismatch: expected {expected}, observed {observed}")
        print(f"UP {NAME} {observed}")
    finally:
        sftp.close()
        client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
