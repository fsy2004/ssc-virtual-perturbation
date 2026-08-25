#!/usr/bin/env python3
"""Independently reconstruct geometry metrics behind O6U review panels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON root is not an object: {path}")
    return value


def verify_record(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    if not path.is_file() or path.stat().st_size != record.get("size_bytes") or sha256(path) != record.get("sha256"):
        raise RuntimeError(f"Artifact failed hash/size verification: {label}")
    return path


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            }
        )
    return atoms


def distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rendering-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-role", required=True, choices=("canary_geometry_review", "formal_geometry_review"))
    args = parser.parse_args()

    rendering_path = args.rendering_report.resolve()
    rendering = load_json(rendering_path)
    if (
        rendering.get("status") != "pass_rendering_only_direct_review_required"
        or rendering.get("role") != args.expected_role
        or rendering.get("production_approved") is not False
        or rendering.get("automatic_decision_applied") is not False
    ):
        raise SystemExit("Rendering report does not pass its exact role/status gate")
    table_path = verify_record(rendering.get("adjudication_table"), "rendering.adjudication_table")
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    row_by_id = {row["orientation_id"]: row for row in rows}
    records = rendering.get("orientations")
    if not isinstance(records, list) or len(records) != 20 or len(row_by_id) != 20:
        raise SystemExit("Rendering report or source table does not contain exactly 20 orientations")

    reconstructed: list[dict[str, object]] = []
    seen: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            raise SystemExit("Rendering orientation record is not an object")
        orientation_id = str(record.get("orientation_id", ""))
        if orientation_id in seen or orientation_id not in row_by_id:
            raise SystemExit(f"Duplicate or unexpected rendering orientation: {orientation_id}")
        seen.add(orientation_id)
        row = row_by_id[orientation_id]
        pdb_path = verify_record(record.get("representative_pdb"), f"{orientation_id}.representative_pdb")
        panel_path = verify_record(record.get("panel"), f"{orientation_id}.panel")
        with Image.open(panel_path) as image:
            image.verify()
        with Image.open(panel_path) as image:
            if image.width < 1000 or image.height < 500:
                raise SystemExit(f"Panel dimensions are unexpectedly small: {orientation_id}")

        atoms = parse_pdb(pdb_path)
        if len(atoms) != 79 or any(atom["resname"] != "O6U" for atom in atoms[:76]) or any(atom["resname"] != "TIP" for atom in atoms[76:]):
            raise SystemExit(f"PDB component identity differs: {orientation_id}")
        ligand_names = [str(atom["name"]) for atom in atoms[:76]]
        target_index = ligand_names.index(row["target_atom"])
        competitor_ligand_index = ligand_names.index(row["nearest_non_target_ligand_atom_at_2p0A"])
        water_indices = list(range(76, 79))
        intended_index = min(water_indices, key=lambda index: distance(atoms[target_index]["xyz"], atoms[index]["xyz"]))
        label = row["nearest_non_target_water_atom_at_2p0A"]
        competitor_candidates = [index for index in water_indices if atoms[index]["name"] == label or atoms[index]["element"] == label]
        if not competitor_candidates:
            raise SystemExit(f"Cannot map competitor water label: {orientation_id}")
        competitor_index = min(competitor_candidates, key=lambda index: distance(atoms[competitor_ligand_index]["xyz"], atoms[index]["xyz"]))
        intended_distance = distance(atoms[target_index]["xyz"], atoms[intended_index]["xyz"])
        competitor_distance = distance(atoms[competitor_ligand_index]["xyz"], atoms[competitor_index]["xyz"])

        if str(atoms[intended_index]["name"]) != record.get("intended_water_atom"):
            raise SystemExit(f"Intended water atom reconstruction differs: {orientation_id}")
        if abs(intended_distance - float(record.get("intended_distance_angstrom_reconstructed"))) > 1e-12:
            raise SystemExit(f"Intended distance reconstruction differs: {orientation_id}")
        if abs(competitor_distance - float(record.get("nearest_competing_distance_angstrom_reconstructed"))) > 1e-12:
            raise SystemExit(f"Competitor distance reconstruction differs: {orientation_id}")
        if str(atoms[competitor_index]["element"]) != row["nearest_non_target_water_atom_at_2p0A"] and str(atoms[competitor_index]["name"]) != row["nearest_non_target_water_atom_at_2p0A"]:
            raise SystemExit(f"Competitor water identity reconstruction differs: {orientation_id}")
        reconstructed.append(
            {
                "orientation_id": orientation_id,
                "intended_water_atom": atoms[intended_index]["name"],
                "intended_distance_angstrom": intended_distance,
                "competitor_water_atom": atoms[competitor_index]["name"],
                "competitor_distance_angstrom": competitor_distance,
                "panel": artifact(panel_path),
            }
        )

    if seen != set(row_by_id):
        raise SystemExit("Rendering orientation universe differs from the source table")
    sheet_path = verify_record(rendering.get("contact_sheet"), "rendering.contact_sheet")
    with Image.open(sheet_path) as sheet:
        sheet.verify()
    with Image.open(sheet_path) as sheet:
        if (sheet.width, sheet.height) != (1200, 2440):
            raise SystemExit("Contact-sheet dimensions differ from the frozen 2x10 layout")

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_review_rendering_independent_validation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_rendering_independent_geometry_reconstruction",
        "expected_role": args.expected_role,
        "production_approved": False,
        "rendering_report": artifact(rendering_path),
        "adjudication_table": artifact(table_path),
        "contact_sheet": artifact(sheet_path),
        "orientation_count": len(reconstructed),
        "orientations": reconstructed,
        "automatic_decision_applied": False,
        "release_boundary": (
            "This report validates image artifacts and independently reconstructs the distances represented by the "
            "panels. It does not validate a human/chemical judgment and does not authorize water-interaction QM."
        ),
    }
    atomic_json(args.report.resolve(), report)
    print(json.dumps({"status": report["status"], "report": str(args.report.resolve()), "sha256": sha256(args.report.resolve())}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
