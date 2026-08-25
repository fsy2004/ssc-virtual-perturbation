#!/usr/bin/env python3
"""Atomically upload the tested energy-term builder, validator, and template."""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
REMOTE_ROOT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
FILES = [
    (ROOT / "scripts" / "build_energy_terms_record.py", "scripts/build_energy_terms_record.py"),
    (ROOT / "scripts" / "gmx_energy_qc.py", "scripts/gmx_energy_qc.py"),
    (ROOT / "templates" / "gromacs_energy_terms.template.json", "templates/gromacs_energy_terms.template.json"),
    (ROOT / "config" / "production_protocol_hmr4fs_303K_v1.json", "config/production_protocol_hmr4fs_303K_v1.json"),
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def remote_digest(sftp, path: str) -> str:
    checksum = hashlib.sha256()
    with sftp.open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            checksum.update(chunk)
    return checksum.hexdigest()


def mkdirs(sftp, path: str) -> None:
    current = ""
    for part in path.strip("/").split("/"):
        current += "/" + part
        try:
            sftp.mkdir(current)
        except OSError:
            pass


def main() -> int:
    client = connect()
    sftp = client.open_sftp()
    try:
        for local, relative in FILES:
            remote = f"{REMOTE_ROOT}/{relative}"
            mkdirs(sftp, remote.rsplit("/", 1)[0])
            temporary = remote + ".uploading"
            sftp.put(str(local), temporary)
            expected = digest(local)
            if remote_digest(sftp, temporary) != expected:
                raise RuntimeError(f"staged hash mismatch: {relative}")
            sftp.posix_rename(temporary, remote)
            if remote_digest(sftp, remote) != expected:
                raise RuntimeError(f"final hash mismatch: {relative}")
            print(f"VERIFIED {relative} {expected}")
    finally:
        sftp.close()
        client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
