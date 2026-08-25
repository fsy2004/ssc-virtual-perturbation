#!/usr/bin/env python3
"""Independently freeze the completed five-member O6U MP2 geometry ensemble."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from prepare_o6u_crest_input import load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble


HARTREE_TO_KCAL_MOL = 627.5094740631
DUPLICATE_RMSD_ANGSTROM = 0.1
DUPLICATE_ENERGY_KCAL_MOL = 0.1
EXPECTED_SOURCE_SDF_SHA256 = "2cb9d769cde4157181a6199b83294cad56cade14ab34a5e86a6deb6790fc28d5"
EXPECTED_SELECTION_SHA256 = "c2afc4b19cf5159f864067b4a2342d57b8f74e59eaddab39ffd948ffe85d80442"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_object(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SystemExit(f"Expected JSON object: {path}")
    return value


def require_file(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SystemExit(f"Missing or empty required artifact: {resolved}")
    return resolved


def kabsch_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 3:
        raise SystemExit("RMSD coordinate shapes differ")
    left_centered = left - left.mean(axis=0)
    right_centered = right - right.mean(axis=0)
    covariance = left_centered.T @ right_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(u @ vt) < 0:
        correction[-1, -1] = -1
    rotation = u @ correction @ vt
    aligned = left_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - right_centered) ** 2, axis=1))))


def validate_member(
    frame: int,
    record_path: Path,
    validation_path: Path,
    role: str,
    source_sdf_sha256: str,
    expected_elements: list[str],
) -> dict[str, object]:
    record_path = require_file(record_path)
    validation_path = require_file(validation_path)
    record = load_object(record_path)
    validation = load_object(validation_path)
    expected_record_status = {
        "canary": "pass_optimization_canary",
        "representative_target": "pass_optimization_representative_target",
    }[role]
    expected_validation_status = {
        "canary": "pass_canary_independently_reconstructed",
        "representative_target": "pass_representative_independently_reconstructed",
    }[role]
    if record.get("schema_version") != "1.1" or record.get("role") != role:
        raise SystemExit(f"Frame {frame} record schema/role differs")
    if record.get("status") != expected_record_status or record.get("production_approved") is not False:
        raise SystemExit(f"Frame {frame} optimization is not a completed non-production pass")
    if record.get("method") != "DF-MP2/6-31G(d), frozen core, RHF reference":
        raise SystemExit(f"Frame {frame} model chemistry differs")
    if record.get("optimizer_strategy_version") != "cartesian_rfo_trust020_v1":
        raise SystemExit(f"Frame {frame} optimizer strategy differs")
    if record.get("source_sdf", {}).get("sha256") != source_sdf_sha256:
        raise SystemExit(f"Frame {frame} source SDF binding differs")
    if validation.get("schema_version") != "1.1" or validation.get("role") != role:
        raise SystemExit(f"Frame {frame} validation schema/role differs")
    if validation.get("status") != expected_validation_status or validation.get("production_approved") is not False:
        raise SystemExit(f"Frame {frame} independent validation is not PASS")
    if validation.get("record_sha256") != sha256(record_path):
        raise SystemExit(f"Frame {frame} validation does not bind its run record")

    optimized = require_file(Path(str(record.get("optimized_xyz", {}).get("path", ""))))
    raw_output = require_file(Path(str(record.get("raw_output", ""))))
    if record.get("optimized_xyz", {}).get("sha256") != sha256(optimized):
        raise SystemExit(f"Frame {frame} optimized XYZ hash differs")
    if record.get("raw_output_sha256") != sha256(raw_output):
        raise SystemExit(f"Frame {frame} raw output hash differs")
    if validation.get("optimized_xyz_sha256") != sha256(optimized):
        raise SystemExit(f"Frame {frame} validation optimized XYZ hash differs")
    if validation.get("raw_output_sha256") != sha256(raw_output):
        raise SystemExit(f"Frame {frame} validation raw output hash differs")

    frames = read_xyz_ensemble(optimized)
    if len(frames) != 1 or frames[0]["elements"] != expected_elements:
        raise SystemExit(f"Frame {frame} optimized XYZ identity/order differs")
    coordinates = np.asarray(frames[0]["coordinates"], dtype=float)
    if coordinates.shape != (len(expected_elements), 3) or not np.isfinite(coordinates).all():
        raise SystemExit(f"Frame {frame} optimized coordinates are invalid")
    energy = float(record.get("final_energy_hartree", math.nan))
    if not math.isfinite(energy) or abs(float(validation.get("final_energy_hartree", math.nan)) - energy) > 1.0e-10:
        raise SystemExit(f"Frame {frame} final energy differs between record and validation")
    return {
        "crest_frame_1based": frame,
        "role": role,
        "final_energy_hartree": energy,
        "record": {"path": str(record_path), "sha256": sha256(record_path)},
        "validation": {"path": str(validation_path), "sha256": sha256(validation_path)},
        "optimized_xyz": {"path": str(optimized), "sha256": sha256(optimized)},
        "raw_output": {"path": str(raw_output), "sha256": sha256(raw_output)},
        "coordinates": coordinates,
    }


def components(frames: list[int], duplicate_pairs: set[tuple[int, int]]) -> list[list[int]]:
    graph = {frame: set() for frame in frames}
    for left, right in duplicate_pairs:
        graph[left].add(right)
        graph[right].add(left)
    seen: set[int] = set()
    result: list[list[int]] = []
    for start in sorted(frames):
        if start in seen:
            continue
        stack = [start]
        group: list[int] = []
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            group.append(current)
            stack.extend(sorted(graph[current] - seen, reverse=True))
        result.append(sorted(group))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-report", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--canary-record", required=True, type=Path)
    parser.add_argument("--canary-validation", required=True, type=Path)
    parser.add_argument("--batch-state", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    selection_path = require_file(args.selection_report)
    source_path = require_file(args.source_sdf)
    canary_record_path = require_file(args.canary_record)
    canary_validation_path = require_file(args.canary_validation)
    batch_path = require_file(args.batch_state)
    report_path = args.report.resolve()
    if report_path.exists():
        raise SystemExit(f"Refusing to overwrite ensemble freeze report: {report_path}")
    if sha256(selection_path) != EXPECTED_SELECTION_SHA256 or sha256(source_path) != EXPECTED_SOURCE_SDF_SHA256:
        raise SystemExit("Frozen selection or source SDF hash differs")

    selection = load_object(selection_path)
    batch = load_object(batch_path)
    if selection.get("status") != "pass" or selection.get("selected_count") != 5:
        raise SystemExit("Representative selection is not the frozen five-member PASS")
    if batch.get("status") != "pass_all_representatives" or batch.get("production_approved") is not False:
        raise SystemExit("Four-member non-canary batch is not complete")
    if batch.get("selection_report", {}).get("sha256") != EXPECTED_SELECTION_SHA256:
        raise SystemExit("Batch state does not bind the frozen selection")
    if batch.get("source_sdf", {}).get("sha256") != EXPECTED_SOURCE_SDF_SHA256:
        raise SystemExit("Batch state does not bind the frozen O6U SDF")

    molecule = load_single_sdf(source_path)
    identity = validate_identity(molecule)
    expected_elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    heavy_indices = [index for index, element in enumerate(expected_elements) if element != "H"]
    global_frame = int(selection["crest_global_minimum_frame_1based"])
    members = [
        validate_member(
            global_frame,
            canary_record_path,
            canary_validation_path,
            "canary",
            EXPECTED_SOURCE_SDF_SHA256,
            expected_elements,
        )
    ]
    targets = batch.get("targets")
    if not isinstance(targets, list) or len(targets) != 4:
        raise SystemExit("Batch state does not contain four targets")
    for target in targets:
        if not isinstance(target, dict) or target.get("status") != "pass":
            raise SystemExit("Batch target is incomplete")
        frame = int(target["crest_frame_1based"])
        record_path = Path(str(target.get("record", "")))
        validation_path = Path(str(target.get("validation", "")))
        if target.get("record_sha256") != sha256(require_file(record_path)):
            raise SystemExit(f"Frame {frame} batch record hash differs")
        if target.get("validation_sha256") != sha256(require_file(validation_path)):
            raise SystemExit(f"Frame {frame} batch validation hash differs")
        members.append(
            validate_member(
                frame,
                record_path,
                validation_path,
                "representative_target",
                EXPECTED_SOURCE_SDF_SHA256,
                expected_elements,
            )
        )
    if len({int(item["crest_frame_1based"]) for item in members}) != 5:
        raise SystemExit("Five optimized representative frame identities are not unique")

    member_by_frame = {int(item["crest_frame_1based"]): item for item in members}
    pairwise: list[dict[str, object]] = []
    duplicate_pairs: set[tuple[int, int]] = set()
    frame_ids = sorted(member_by_frame)
    for left_index, left in enumerate(frame_ids):
        for right in frame_ids[left_index + 1 :]:
            left_member = member_by_frame[left]
            right_member = member_by_frame[right]
            rmsd = kabsch_rmsd(
                np.asarray(left_member["coordinates"])[heavy_indices],
                np.asarray(right_member["coordinates"])[heavy_indices],
            )
            delta = abs(float(left_member["final_energy_hartree"]) - float(right_member["final_energy_hartree"])) * HARTREE_TO_KCAL_MOL
            duplicate = rmsd < DUPLICATE_RMSD_ANGSTROM and delta < DUPLICATE_ENERGY_KCAL_MOL
            if duplicate:
                duplicate_pairs.add((left, right))
            pairwise.append(
                {
                    "left_frame_1based": left,
                    "right_frame_1based": right,
                    "heavy_atom_kabsch_rmsd_angstrom": rmsd,
                    "absolute_energy_difference_kcal_mol": delta,
                    "duplicate_requires_both_strict_criteria": duplicate,
                }
            )

    groups = components(frame_ids, duplicate_pairs)
    collapsed_groups: list[dict[str, object]] = []
    representatives: list[int] = []
    for group in groups:
        retained = min(group, key=lambda frame: (float(member_by_frame[frame]["final_energy_hartree"]), frame))
        representatives.append(retained)
        collapsed_groups.append(
            {
                "member_frames_1based": group,
                "retained_frame_1based": retained,
                "rule": "lowest MP2 energy, then lowest frame number",
            }
        )
    charge_target = min(
        representatives,
        key=lambda frame: (float(member_by_frame[frame]["final_energy_hartree"]), frame),
    )

    serializable_members = []
    for frame in frame_ids:
        item = dict(member_by_frame[frame])
        item.pop("coordinates")
        serializable_members.append(item)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_mp2_representative_ensemble_freeze",
        "status": "pass_five_member_ensemble_independently_reconstructed",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "selection_report": {"path": str(selection_path), "sha256": sha256(selection_path)},
        "source_sdf": {"path": str(source_path), "sha256": sha256(source_path)},
        "canary_record": {"path": str(canary_record_path), "sha256": sha256(canary_record_path)},
        "canary_validation": {"path": str(canary_validation_path), "sha256": sha256(canary_validation_path)},
        "batch_state": {"path": str(batch_path), "sha256": sha256(batch_path)},
        "duplicate_rule": {
            "heavy_atom_kabsch_rmsd_angstrom_strictly_below": DUPLICATE_RMSD_ANGSTROM,
            "absolute_energy_difference_kcal_mol_strictly_below": DUPLICATE_ENERGY_KCAL_MOL,
            "both_required": True,
        },
        "members": serializable_members,
        "pairwise": pairwise,
        "duplicate_groups": collapsed_groups,
        "unique_representative_frames_1based": sorted(representatives),
        "charge_water_target_frame_1based": charge_target,
        "charge_water_target_optimized_xyz": member_by_frame[charge_target]["optimized_xyz"],
        "charge_water_target_rule": "lowest-energy unique MP2/6-31G(d) optimized representative; ties resolve by frame number",
        "release_boundary": (
            "This freezes the optimized QM geometry ensemble and one geometry for charge/water target generation. "
            "It does not approve any fitted parameter, CHARMM-GUI build, or MD stage."
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "report_sha256": sha256(report_path), "unique_count": len(representatives), "charge_target_frame": charge_target}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
