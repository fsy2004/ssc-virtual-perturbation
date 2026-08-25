#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

from new_md_server import connect, run


ROOT = Path(__file__).resolve().parent.parent
REMOTE_ROOT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
FILES = (
    "scripts/install_gorder_1_5_0.sh",
    "scripts/install_fatslim_0_2_2.sh",
    "scripts/_remote_launch_gorder_install.sh",
    "scripts/_remote_launch_fatslim_install.sh",
    "scripts/_remote_pause_postproduction_for_cpu_switch.sh",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    client = connect()
    try:
        code, _, error = run(client, f"mkdir -p {REMOTE_ROOT}/scripts", timeout=60)
        if code:
            raise RuntimeError(error)
        sftp = client.open_sftp()
        try:
            for relative in FILES:
                local = ROOT / relative
                remote = f"{REMOTE_ROOT}/{relative}"
                sftp.put(str(local), remote)
                sftp.chmod(remote, 0o755)
                expected = sha256(local)
                code, output, error = run(client, f"sha256sum {remote}", timeout=60)
                if code:
                    raise RuntimeError(error)
                actual = output.split()[0]
                if actual != expected:
                    raise RuntimeError(f"SHA256 mismatch for {relative}: {actual} != {expected}")
                print(f"{expected}  {remote}")
        finally:
            sftp.close()
    finally:
        client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
