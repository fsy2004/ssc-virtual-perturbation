#!/usr/bin/env python3
"""Bind the frozen primary-analysis manifest to completed, hash-verified files."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPLICAS = ("rep01", "rep02", "rep03")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return payload


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical(root: Path, path: Path | str) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    candidate = candidate.resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"artifact escapes release root: {candidate}") from exc
    return candidate


def relative_record(root: Path, path: Path) -> dict[str, Any]:
    path = canonical(root, path)
    require(path.is_file() and path.stat().st_size > 0, f"missing or empty artifact: {path}")
    return {
        "path": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def verify_report_artifact(root: Path, report: dict[str, Any], rid: str, name: str) -> dict[str, Any]:
    expected_path = root / rid / "work" / name
    live = relative_record(root, expected_path)
    recorded = report.get("artifacts", {}).get(name)
    require(isinstance(recorded, dict), f"{rid} completion report lacks {name}")
    require(recorded.get("sha256") == live["sha256"], f"{rid} completion-report hash differs for {name}")
    require(int(recorded.get("bytes", -1)) == live["bytes"], f"{rid} completion-report byte count differs for {name}")
    return live


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def bind_manifest(draft_path: Path, release_root: Path, protocol_path: Path, output_path: Path) -> dict[str, Any]:
    root = release_root.resolve()
    draft_path = draft_path.resolve()
    protocol_path = protocol_path.resolve()
    output_path = output_path.resolve()
    require(root.is_dir(), f"release root is missing: {root}")
    require(not output_path.exists(), f"refusing to overwrite sealed manifest: {output_path}")

    payload = load_json(draft_path)
    require(payload.get("approval_status") == "draft_not_for_execution", "input manifest is not the frozen draft")
    realizations = payload.get("realizations")
    require(isinstance(realizations, list), "draft manifest lacks realizations")
    require([item.get("realization_id") for item in realizations] == list(REPLICAS), "draft realization order differs")

    protocol = load_json(protocol_path)
    protocol_replicas = protocol.get("realizations")
    require(isinstance(protocol_replicas, list), "production protocol lacks realizations")
    seeds = {item.get("id"): item.get("velocity_seed") for item in protocol_replicas}
    require(set(seeds) == set(REPLICAS), "production protocol realization identities differ")
    require(all(isinstance(seeds[rid], int) and seeds[rid] > 0 for rid in REPLICAS), "velocity seeds are invalid")
    require(len({seeds[rid] for rid in REPLICAS}) == 3, "velocity seeds are not unique")

    binding_evidence: list[dict[str, Any]] = []
    for item, rid in zip(realizations, REPLICAS, strict=True):
        completion_path = root / rid / "PRODUCTION_COMPLETION_500NS.json"
        provenance_path = root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / rid / "trajectory_provenance.pre_qc.json"
        completion = load_json(completion_path)
        provenance = load_json(provenance_path)
        require(completion.get("status") == "pass", f"{rid} 500 ns completion gate did not pass")
        require(provenance.get("status") == "pass_pending_scientific_qc_seal", f"{rid} PBC provenance status differs")
        require(provenance.get("realization_id") == rid, f"{rid} PBC provenance identity differs")

        topology = verify_report_artifact(root, completion, rid, "production.tpr")
        energy = verify_report_artifact(root, completion, rid, "production.edr")
        log = verify_report_artifact(root, completion, rid, "production.log")
        require(provenance.get("production_tpr_sha256") == topology["sha256"], f"{rid} TPR hash differs between completion and PBC provenance")

        centered_record = provenance.get("retained_outputs", {}).get("center_and_rebox")
        require(isinstance(centered_record, dict), f"{rid} PBC provenance lacks centered trajectory")
        centered_path = canonical(root, centered_record.get("path", ""))
        expected_centered = (root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / rid / "05_centered_reboxed.xtc").resolve()
        require(centered_path == expected_centered, f"{rid} centered trajectory path differs")
        centered = relative_record(root, centered_path)
        require(centered_record.get("sha256") == centered["sha256"], f"{rid} centered trajectory hash differs")
        require(int(centered_record.get("bytes", -1)) == centered["bytes"], f"{rid} centered trajectory byte count differs")

        item.update({
            "velocity_seed": seeds[rid],
            "topology": topology,
            "centered_system_trajectory": centered,
            "energy_edr": energy,
            "production_log": log,
        })
        binding_evidence.append({
            "realization_id": rid,
            "completion_report": relative_record(root, completion_path),
            "pbc_provenance": relative_record(root, provenance_path),
        })

    payload["approval_status"] = "approved_for_server_execution"
    acceptance = payload.get("acceptance_gates")
    require(isinstance(acceptance, dict), "draft manifest lacks acceptance gates")
    acceptance["approval_status"] = "approved_and_frozen_before_production"
    payload["binding_evidence"] = binding_evidence
    payload["sealed_at_utc"] = datetime.now(timezone.utc).isoformat()
    payload["seal_policy"] = "hash_bound_from_all_three_500ns_completion_and_pbc_provenance_records"

    text = json.dumps(payload, indent=2) + "\n"
    atomic_write(output_path, text)
    digest = sha256(output_path)
    atomic_write(output_path.with_suffix(".json.sha256"), f"{digest}  {output_path.name}\n")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draft", required=True, type=Path)
    parser.add_argument("--release-root", required=True, type=Path)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = bind_manifest(args.draft, args.release_root, args.protocol, args.output)
    print(json.dumps({
        "status": payload["approval_status"],
        "output": str(args.output.resolve()),
        "sha256": sha256(args.output.resolve()),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
