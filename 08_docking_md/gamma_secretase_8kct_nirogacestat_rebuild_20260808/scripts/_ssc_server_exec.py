#!/usr/bin/env python3
"""Upload-free server command runner: python scripts/_ssc_server_exec.py "<bash command>"
If the argument starts with '@', the command is read from the file after '@' (local path)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect, run  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print("usage: _ssc_server_exec.py '<command>' | '@command-file'", file=sys.stderr)
        return 2
    arg = args[0]
    if arg.startswith("@"):
        cmd = Path(arg[1:]).read_text(encoding="utf-8")
    else:
        cmd = arg
    c = connect()
    try:
        code, out, err = run(c, cmd, timeout=900)
        print(out, end="")
        if err:
            print("STDERR:", err, end="")
    finally:
        c.close()
    return code


if __name__ == "__main__":
    raise SystemExit(main())
