#!/usr/bin/env python3
"""Download and verify the fragment-closed formal analysis index."""
from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
REMOTE_ROOT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/builds"
FILES = [
    ("analysis.ndx", ROOT / "builds" / "analysis.ndx", "0c5b58650756bfa18d6a2b00b3688b697f42f6e64ea267f6cdeceb44e132189e"),
    ("analysis.ndx.provenance.json", ROOT / "builds" / "analysis.ndx.provenance.json", "4df5ada7d317aec6b0283fa1fed4b69b45b4e967c61166e4a7c6979bcb3d201f"),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    client = connect()
    sftp = client.open_sftp()
    try:
        for name, local, expected in FILES:
            temporary = local.with_suffix(local.suffix + ".downloading")
            sftp.get(f"{REMOTE_ROOT}/{name}", str(temporary))
            observed = sha256(temporary)
            if observed != expected:
                temporary.unlink(missing_ok=True)
                raise RuntimeError(f"hash mismatch for {name}: {observed}")
            os.replace(temporary, local)
            print(f"DOWNLOADED {name} sha256={observed}")
    finally:
        sftp.close()
        client.close()
    (ROOT / "builds" / "analysis.ndx.sha256").write_text(
        f"{FILES[0][2]}  analysis.ndx\n", encoding="ascii"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
