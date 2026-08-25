#!/usr/bin/env python3
"""Read-only server inspection helper (never prints credentials)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect, run

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "pwd"
    c = connect()
    try:
        code, out, err = run(c, cmd, timeout=600)
        print(out)
        if err:
            print("STDERR:", err)
    finally:
        c.close()
