#!/usr/bin/env python3
"""Independently reconstruct an O6U FFParam water-geometry audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import runpy
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_ORIENTATIONS = 20
EXPECTED_GRID = [value / 100 for value in range(150, 301, 5)]
REPRESENTATIVE_DISTANCE = 2.0
INTENDED_TOLERANCE = 2.0e-4
COLLISION_DISTANCE = 1.2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def parse_block(block: str) -> list[tuple[str, tuple[float, float, float]]]:
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    if not lines or len(lines[0].split()) != 2:
        raise SystemExit("Coordinate block lacks charge and multiplicity")
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 4:
            raise SystemExit(f"Invalid coordinate line: {line}")
        atoms.append((fields[0].upper(), tuple(float(value) for value in fields[1:])))
    return atoms


def distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right, strict=True)))


def rtf_atoms(path: Path) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] == "ATOM":
            result.append({"name": fields[1], "element": fields[1][0].upper()})
    if len(result) != 76:
        raise SystemExit(f"Expected 76 RTF atoms, found {len(result)}")
    return result


def close(observed: float, expected: float, tolerance: float = 1.0e-9) -> bool:
    return math.isfinite(observed) and abs(observed - expected) <= tolerance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-report", required=True, type=Path)
    parser.add_argument("--generation-report", required=True, type=Path)
    parser.add_argument("--rtf", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    audit_path = args.audit_report.resolve()
    generation_path = args.generation_report.resolve()
    rtf_path = args.rtf.resolve()
    report_path = args.report.resolve()
    for path in (audit_path, generation_path, rtf_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty input: {path}")
    if report_path.exists():
        raise SystemExit(f"Refusing to overwrite report: {report_path}")

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    generation = json.loads(generation_path.read_text(encoding="utf-8"))
    if audit.get("status") != "pass_geometry_integrity_visual_review_required":
        raise SystemExit("Geometry audit did not pass its integrity-only gate")
    if audit.get("production_approved") is not False or audit.get("orientation_count") != EXPECTED_ORIENTATIONS:
        raise SystemExit("Geometry audit improperly claims approval or has the wrong count")
    if generation.get("status") != "pass_generation_only_visual_review_required":
        raise SystemExit("Input-generation report did not pass its generation-only gate")
    if generation.get("production_approved") is not False or generation.get("orientation_count") != EXPECTED_ORIENTATIONS:
        raise SystemExit("Input-generation report improperly claims approval or has the wrong count")
    audit_generation = audit.get("inputs", {}).get("generation_report", {})
    audit_rtf = audit.get("inputs", {}).get("rtf", {})
    if audit_generation.get("sha256") != sha256(generation_path):
        raise SystemExit("Audit does not bind the supplied generation report")
    if audit_rtf.get("sha256") != sha256(rtf_path):
        raise SystemExit("Audit does not bind the supplied RTF")
    if not close(float(audit.get("representative_distance_angstrom", math.nan)), REPRESENTATIVE_DISTANCE):
        raise SystemExit("Audit representative distance differs")
    if not close(float(audit.get("distance_integrity_tolerance_angstrom", math.nan)), INTENDED_TOLERANCE):
        raise SystemExit("Audit intended-distance tolerance differs")
    if not close(float(audit.get("sanity_collision_distance_angstrom", math.nan)), COLLISION_DISTANCE):
        raise SystemExit("Audit sanity-collision distance differs")

    topology = rtf_atoms(rtf_path)
    generated = {item["orientation_id"]: item for item in generation["generated_pairs"]}
    audited = {item["orientation_id"]: item for item in audit["orientations"]}
    if len(generated) != EXPECTED_ORIENTATIONS or set(generated) != set(audited):
        raise SystemExit("Generation and audit orientation identities differ")

    reconstructed: list[dict[str, object]] = []
    for orientation_id in sorted(generated):
        source = generated[orientation_id]
        observed = audited[orientation_id]
        if source.get("source_definition") != observed.get("source_definition"):
            raise SystemExit(f"Source definition differs for {orientation_id}")
        fields = str(source["source_definition"]).split()
        target_name = fields[2]
        target_indices = [index for index, atom in enumerate(topology) if atom["name"] == target_name]
        if len(target_indices) != 1 or observed.get("target_atom") != target_name:
            raise SystemExit(f"Target atom identity differs for {orientation_id}")
        target_index = target_indices[0]

        coordinate_path = Path(str(source["coordinate_file"]["path"])).resolve()
        if not coordinate_path.is_file() or sha256(coordinate_path) != source["coordinate_file"]["sha256"]:
            raise SystemExit(f"Coordinate module integrity differs for {orientation_id}")
        if observed.get("coordinate_module", {}).get("sha256") != sha256(coordinate_path):
            raise SystemExit(f"Audit coordinate-module binding differs for {orientation_id}")
        module = runpy.run_path(str(coordinate_path))
        ligand = parse_block(str(module.get("basecoor", "")))
        if len(ligand) != len(topology):
            raise SystemExit(f"Ligand atom count differs for {orientation_id}")
        if [atom[0] for atom in ligand] != [atom["element"] for atom in topology]:
            raise SystemExit(f"Ligand element order differs for {orientation_id}")
        labels = list(module.get("intrange", []))
        grid = [float(str(label).replace("_", ".")) for label in labels]
        if grid != EXPECTED_GRID or list(source.get("distance_grid_angstrom", [])) != EXPECTED_GRID:
            raise SystemExit(f"Distance grid differs for {orientation_id}")

        target_xyz = ligand[target_index][1]
        scan_rows: list[dict[str, object]] = []
        representative: dict[str, object] | None = None
        for label, scan_distance in zip(labels, grid, strict=True):
            water = parse_block(str(module.get(f"interaction_{label}", "")))
            if len(water) != 3 or sorted(atom[0] for atom in water) != ["H", "H", "O"]:
                raise SystemExit(f"Water identity differs for {orientation_id} at {scan_distance}")
            if target_name.startswith("H"):
                intended_index = next(index for index, atom in enumerate(water) if atom[0] == "O")
            else:
                h_indices = [index for index, atom in enumerate(water) if atom[0] == "H"]
                intended_index = min(h_indices, key=lambda index: distance(target_xyz, water[index][1]))
            intended_distance = distance(target_xyz, water[intended_index][1])
            if abs(intended_distance - scan_distance) > INTENDED_TOLERANCE:
                raise SystemExit(f"Intended distance differs for {orientation_id} at {scan_distance}")
            contacts = []
            for ligand_index, ligand_atom in enumerate(ligand):
                if ligand_index == target_index:
                    continue
                for water_index, water_atom in enumerate(water):
                    contacts.append((distance(ligand_atom[1], water_atom[1]), ligand_index, water_index))
            nearest, ligand_index, water_index = min(contacts)
            row = {
                "scan_distance_angstrom": scan_distance,
                "intended_distance_angstrom": intended_distance,
                "nearest_non_target_distance_angstrom": nearest,
                "nearest_non_target_ligand_atom": topology[ligand_index]["name"],
                "nearest_non_target_water_atom": water[water_index][0],
                "sanity_collision_below_1p2A": nearest < COLLISION_DISTANCE,
            }
            scan_rows.append(row)
            if scan_distance == REPRESENTATIVE_DISTANCE:
                representative = row
        if representative is None:
            raise SystemExit(f"Representative scan point missing for {orientation_id}")

        observed_scan = observed.get("scan_geometry")
        if not isinstance(observed_scan, list) or len(observed_scan) != len(scan_rows):
            raise SystemExit(f"Audit scan row count differs for {orientation_id}")
        for expected_row, observed_row in zip(scan_rows, observed_scan, strict=True):
            for key in ("scan_distance_angstrom", "intended_distance_angstrom", "nearest_non_target_distance_angstrom"):
                if not close(float(observed_row.get(key, math.nan)), float(expected_row[key]), 1.0e-10):
                    raise SystemExit(f"Audit numeric reconstruction differs for {orientation_id}: {key}")
            for key in ("nearest_non_target_ligand_atom", "nearest_non_target_water_atom", "sanity_collision_below_1p2A"):
                if observed_row.get(key) != expected_row[key]:
                    raise SystemExit(f"Audit categorical reconstruction differs for {orientation_id}: {key}")
        expected_minimum = min(float(row["nearest_non_target_distance_angstrom"]) for row in scan_rows)
        expected_collision = any(bool(row["sanity_collision_below_1p2A"]) for row in scan_rows)
        if not close(float(observed.get("minimum_non_target_distance_over_scan_angstrom", math.nan)), expected_minimum, 1.0e-10):
            raise SystemExit(f"Audit scan minimum differs for {orientation_id}")
        if observed.get("sanity_collision_anywhere_in_scan") is not expected_collision:
            raise SystemExit(f"Audit collision summary differs for {orientation_id}")
        if observed.get("disposition") != "pending_geometry_specific_visual_review":
            raise SystemExit(f"Audit prematurely assigned a disposition for {orientation_id}")

        pdb_record = observed.get("representative_2p0A", {}).get("pdb", {})
        pdb_path = Path(str(pdb_record.get("path", ""))).resolve()
        if not pdb_path.is_file() or sha256(pdb_path) != pdb_record.get("sha256"):
            raise SystemExit(f"Representative PDB integrity differs for {orientation_id}")
        hetatm_count = sum(line.startswith("HETATM") for line in pdb_path.read_text(encoding="ascii").splitlines())
        if hetatm_count != len(topology) + 3:
            raise SystemExit(f"Representative PDB atom count differs for {orientation_id}")
        reconstructed.append(
            {
                "orientation_id": orientation_id,
                "target_atom": target_name,
                "minimum_non_target_distance_over_scan_angstrom": expected_minimum,
                "sanity_collision_anywhere_in_scan": expected_collision,
                "representative_pdb_sha256": sha256(pdb_path),
            }
        )

    report = {
        "schema_version": "1.0",
        "report_type": "independent_o6u_ffparam_water_geometry_audit_validation",
        "status": "pass_geometry_audit_independently_reconstructed",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_report": artifact(audit_path),
        "generation_report": artifact(generation_path),
        "rtf": artifact(rtf_path),
        "orientation_count": len(reconstructed),
        "reconstructed_orientations": reconstructed,
        "release_boundary": (
            "This independently validates geometry bookkeeping only. All dispositions remain pending, and no "
            "water-interaction QM, fitted parameter, CHARMM-GUI build, or MD stage is authorized."
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "orientation_count": len(reconstructed), "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
