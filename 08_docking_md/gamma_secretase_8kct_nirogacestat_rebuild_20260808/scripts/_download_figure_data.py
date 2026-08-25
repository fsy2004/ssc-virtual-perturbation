#!/usr/bin/env python3
"""Download figure7_data.json to the local project for plotting."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect  # noqa: E402

REMOTE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/analysis/primary_postprocessing/20260822T170054Z_primary_qc_rep01_completion/figure7_data.json"
DEST = Path(__file__).resolve().parent.parent / "analysis" / "postprocessing" / "figure7_data.json"
DEST.parent.mkdir(parents=True, exist_ok=True)


def main() -> int:
    c = connect()
    try:
        sftp = c.open_sftp()
        with sftp.open(REMOTE, "rb") as f:
            data = f.read()
        DEST.write_bytes(data)
        sftp.close()
        print(f"downloaded {len(data)} bytes -> {DEST}")
        return 0
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
