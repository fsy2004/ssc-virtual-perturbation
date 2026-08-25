#!/usr/bin/env python3
"""Upload prepare_figure7_data.py and run it on the server."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect, run  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
REMOTE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/scripts/prepare_figure7_data.py"


def main() -> int:
    c = connect()
    try:
        sftp = c.open_sftp()
        local = ROOT / "scripts" / "prepare_figure7_data.py"
        sftp.put(str(local), REMOTE)
        sftp.close()
        code, out, err = run(
            c,
            "/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python "
            + REMOTE,
            timeout=300,
        )
        print(out, end="")
        if err:
            print("STDERR:", err, end="")
        return code
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
