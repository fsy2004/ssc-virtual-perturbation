#!/usr/bin/env python3
"""Download hash-verified formal mapping artifacts from the active release."""
from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
REMOTE_ROOT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
ARTIFACTS = [
    (
        "config/primary_atom_mapping_contacts.json",
        ROOT / "config" / "primary_atom_mapping_contacts.json",
        "16815da6dd38e3fb6a1e08ff42ff29a011755e1afe182931785d2bdba2abb831",
    ),
    (
        "config/membrane_qc_mapping.json",
        ROOT / "config" / "membrane_qc_mapping.json",
        "5b90b6b7faf8387edc72507204de19fc77abcdba49b3aad87386459bf0c0f630",
    ),
    (
        "audit/mapping_revisions/20260822T152808+0800_pre_formal_tpr_mapping/FORMAL_MAPPING_VALIDATION.json",
        ROOT / "reports" / "postproduction" / "FORMAL_MAPPING_VALIDATION_20260822.json",
        "67c8ec3226e86566f047e0f351e73eea17c949b84e89f188a4a5bbb72903d45b",
    ),
    (
        "config/gromacs_energy_terms.json",
        ROOT / "config" / "gromacs_energy_terms.json",
        "d390b7f9fe4d2bf6feb9521a898b250a7242a5699d926755ec6b74127b05aad8",
    ),
    (
        "audit/energy_term_inventory/ENERGY_TERM_RECORD_VALIDATION.json",
        ROOT / "reports" / "postproduction" / "ENERGY_TERM_RECORD_VALIDATION_20260822.json",
        "005a959db5228ed9be4bc0e814e6df0a3b2a863481bbae0552b6453679f57cba",
    ),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    client = connect()
    sftp = client.open_sftp()
    try:
        for remote_rel, local, expected in ARTIFACTS:
            local.parent.mkdir(parents=True, exist_ok=True)
            temporary = local.with_suffix(local.suffix + ".downloading")
            sftp.get(f"{REMOTE_ROOT}/{remote_rel}", str(temporary))
            observed = sha256(temporary)
            if observed != expected:
                temporary.unlink(missing_ok=True)
                raise RuntimeError(f"download hash mismatch for {remote_rel}: {observed}")
            os.replace(temporary, local)
            sidecar = local.with_suffix(local.suffix + ".sha256")
            sidecar.write_text(f"{observed}  {local.name}\n", encoding="ascii")
            print(f"DOWNLOADED {remote_rel} -> {local} sha256={observed}")
    finally:
        sftp.close()
        client.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
