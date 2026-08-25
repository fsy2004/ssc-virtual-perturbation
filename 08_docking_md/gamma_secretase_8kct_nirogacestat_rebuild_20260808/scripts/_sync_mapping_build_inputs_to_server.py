#!/usr/bin/env python3
"""Upload only tested mapping builders and immutable mapping source inputs."""
from __future__ import annotations

import hashlib
import posixpath
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
REMOTE_ROOT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
FILES = [
    ("scripts/build_analysis_ndx.py", "scripts/build_analysis_ndx.py"),
    ("scripts/build_primary_mapping_records.py", "scripts/build_primary_mapping_records.py"),
    ("scripts/build_membrane_mapping.py", "scripts/build_membrane_mapping.py"),
    ("scripts/analyze_primary_structure_mdanalysis.py", "scripts/analyze_primary_structure_mdanalysis.py"),
    ("scripts/analyze_membrane_qc_mdanalysis.py", "scripts/analyze_membrane_qc_mdanalysis.py"),
    ("scripts/primary_postprocessing_common.py", "scripts/primary_postprocessing_common.py"),
    ("docking_native_redock/plip_native/8KCT_protonated.pdb", "docking_native_redock/plip_native/8KCT_protonated.pdb"),
    ("inputs/ligand_parameterization/O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv", "inputs/ligand_parameterization/O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv"),
    ("docking_native_redock/plip_native/run1/8KCT_O6U.xml", "docking_native_redock/plip_native/run1/8KCT_O6U.xml"),
    (
        "docking_native_redock/figures/native_8kct_o6u/8KCT_O6U_native_contacts.interactions.normalized.json",
        "docking_native_redock/figures/native_8kct_o6u/8KCT_O6U_native_contacts.interactions.normalized.json",
    ),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def remote_sha256(sftp, path: str) -> str:
    digest = hashlib.sha256()
    with sftp.open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def mkdirs(sftp, directory: str) -> None:
    parts = directory.strip("/").split("/")
    current = ""
    for part in parts:
        current += "/" + part
        try:
            sftp.mkdir(current)
        except OSError:
            pass


def main() -> int:
    client = connect()
    sftp = client.open_sftp()
    results: list[tuple[str, str]] = []
    try:
        for local_rel, remote_rel in FILES:
            local = ROOT / local_rel
            if not local.is_file():
                raise FileNotFoundError(local)
            remote = posixpath.join(REMOTE_ROOT, remote_rel)
            mkdirs(sftp, posixpath.dirname(remote))
            temporary = remote + ".uploading"
            sftp.put(str(local), temporary)
            expected = sha256(local)
            if remote_sha256(sftp, temporary) != expected:
                raise RuntimeError(f"staged hash mismatch: {remote}")
            sftp.posix_rename(temporary, remote)
            observed = remote_sha256(sftp, remote)
            if observed != expected:
                raise RuntimeError(f"final hash mismatch: {remote}")
            results.append((remote_rel, observed))
            print(f"VERIFIED {remote_rel} {observed}")
    finally:
        sftp.close()
        client.close()
    print(f"MAPPING INPUT SYNC PASS: {len(results)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
