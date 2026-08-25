#!/usr/bin/env python3
"""Seal a completed GPU release after adding dynamic validation reports."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("release", type=Path)
    args = parser.parse_args()
    release = args.release.resolve()
    protocol = json.loads((release / "production_protocol.json").read_text(encoding="utf-8"))
    canary_path = release / "CANARY_VALIDATION.json"
    canary = json.loads(canary_path.read_text(encoding="utf-8"))
    if canary.get("pass") is not True:
        raise ValueError("CANARY_VALIDATION.json is not a passing report")

    artifacts = []
    for path in sorted(release.rglob("*")):
        if path.is_file() and path.name not in {"RELEASE_MANIFEST.json", "RELEASE_MANIFEST.sha256"}:
            artifacts.append({
                "path": path.relative_to(release).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    manifest = {
        "schema_version": "1.1",
        "release_status": "gpu-ready-after-passing-canary",
        "protocol_id": protocol["protocol_id"],
        "archive_sha256": protocol["system"]["archive_sha256"],
        "canary_validation_sha256": sha256(canary_path),
        "realizations": protocol["realizations"],
        "artifacts": artifacts,
    }
    payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    part = release / "RELEASE_MANIFEST.json.part"
    part.write_bytes(payload)
    os.replace(part, release / "RELEASE_MANIFEST.json")
    manifest_sha = sha256(release / "RELEASE_MANIFEST.json")
    (release / "RELEASE_MANIFEST.sha256").write_text(
        f"{manifest_sha}  RELEASE_MANIFEST.json\n", encoding="utf-8", newline="\n"
    )
    print(manifest_sha)


if __name__ == "__main__":
    main()
