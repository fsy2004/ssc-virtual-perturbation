#!/usr/bin/env python3
"""Audit generated O6U-water scan geometry without selecting or running probes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import runpy
from datetime import datetime, timezone
from pathlib import Path


REPRESENTATIVE_DISTANCE = 2.0
DISTANCE_TOLERANCE = 2.0e-4
SANITY_COLLISION_DISTANCE = 1.2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def parse_xyz_block(block: str) -> list[tuple[str, tuple[float, float, float]]]:
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    if not lines or len(lines[0].split()) != 2:
        raise ValueError("coordinate block lacks charge and multiplicity")
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"invalid coordinate line: {line}")
        atoms.append((fields[0], tuple(float(value) for value in fields[1:])))
    return atoms


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def rtf_atoms(path: Path) -> list[dict[str, str]]:
    atoms: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] == "ATOM":
            atoms.append({"name": fields[1], "type": fields[2], "element": fields[1][0].upper()})
    if not atoms:
        raise SystemExit("RTF contains no ATOM records")
    return atoms


def pdb_atom_line(serial: int, name: str, residue: str, chain: str, residue_id: int,
                  xyz: tuple[float, float, float], element: str) -> str:
    x, y, z = xyz
    return (
        f"HETATM{serial:5d} {name:>4s} {residue:>3s} {chain}{residue_id:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generation-report", required=True, type=Path)
    parser.add_argument("--rtf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    report_path = args.generation_report.resolve()
    rtf_path = args.rtf.resolve()
    if not report_path.is_file() or not rtf_path.is_file():
        raise SystemExit("Missing generation report or RTF")
    generation = json.loads(report_path.read_text(encoding="utf-8"))
    if generation.get("status") != "pass_generation_only_visual_review_required":
        raise SystemExit("Generation report has not passed its generation-only gate")
    if generation.get("production_approved") is not False:
        raise SystemExit("Generation report improperly claims production approval")
    pairs = generation.get("generated_pairs")
    if not isinstance(pairs, list) or len(pairs) != 20:
        raise SystemExit("Expected exactly 20 prescreen-retained orientations")

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise SystemExit(f"Refusing to reuse output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    pdb_dir = output_dir / "representative_2p0A_pdb"
    pdb_dir.mkdir()

    topology = rtf_atoms(rtf_path)
    audited: list[dict[str, object]] = []
    seen_ids: set[str] = set()
    for pair in pairs:
        orientation_id = pair.get("orientation_id")
        definition = pair.get("source_definition")
        if not isinstance(orientation_id, str) or not isinstance(definition, str):
            raise SystemExit("Orientation record lacks identity or source definition")
        if orientation_id in seen_ids:
            raise SystemExit(f"Duplicate orientation identity: {orientation_id}")
        seen_ids.add(orientation_id)
        definition_fields = definition.split()
        if len(definition_fields) < 3:
            raise SystemExit(f"Cannot parse orientation definition: {definition}")
        target_name = definition_fields[2]
        target_indices = [i for i, atom in enumerate(topology) if atom["name"] == target_name]
        if len(target_indices) != 1:
            raise SystemExit(f"Target atom is not unique in RTF: {target_name}")
        target_index = target_indices[0]

        coordinate_path = Path(pair["coordinate_file"]["path"]).resolve()
        if not coordinate_path.is_file() or sha256(coordinate_path) != pair["coordinate_file"]["sha256"]:
            raise SystemExit(f"Coordinate module integrity differs: {coordinate_path}")
        module = runpy.run_path(str(coordinate_path))
        ligand = parse_xyz_block(module.get("basecoor", ""))
        if len(ligand) != len(topology):
            raise SystemExit(f"Ligand atom count differs for {orientation_id}")
        for index, ((element, _), atom) in enumerate(zip(ligand, topology, strict=True)):
            if element.upper() != atom["element"]:
                raise SystemExit(f"RTF-coordinate element mismatch at atom {index + 1} for {orientation_id}")

        labels = list(module.get("intrange", []))
        expected_grid = list(pair.get("distance_grid_angstrom", []))
        observed_grid = [float(label.replace("_", ".")) for label in labels]
        if observed_grid != expected_grid or REPRESENTATIVE_DISTANCE not in observed_grid:
            raise SystemExit(f"Distance grid/report mismatch for {orientation_id}")

        scan_records: list[dict[str, object]] = []
        representative: dict[str, object] | None = None
        target_xyz = ligand[target_index][1]
        for label, scan_distance in zip(labels, observed_grid, strict=True):
            water = parse_xyz_block(module.get(f"interaction_{label}", ""))
            if len(water) != 3 or sorted(element.upper() for element, _ in water) != ["H", "H", "O"]:
                raise SystemExit(f"Water geometry differs for {orientation_id} at {scan_distance:.2f} A")
            if target_name.startswith("H"):
                intended_water_index = next(i for i, (element, _) in enumerate(water) if element.upper() == "O")
            else:
                hydrogen_indices = [i for i, (element, _) in enumerate(water) if element.upper() == "H"]
                intended_water_index = min(hydrogen_indices, key=lambda i: distance(target_xyz, water[i][1]))
            intended_distance = distance(target_xyz, water[intended_water_index][1])
            if abs(intended_distance - scan_distance) > DISTANCE_TOLERANCE:
                raise SystemExit(
                    f"Intended target-water distance differs for {orientation_id}: "
                    f"{intended_distance:.6f} vs {scan_distance:.6f} A"
                )

            non_target = []
            for ligand_index, (_, ligand_xyz) in enumerate(ligand):
                if ligand_index == target_index:
                    continue
                for water_index, (_, water_xyz) in enumerate(water):
                    non_target.append((distance(ligand_xyz, water_xyz), ligand_index, water_index))
            nearest_distance, nearest_ligand_index, nearest_water_index = min(non_target)
            row = {
                "scan_distance_angstrom": scan_distance,
                "intended_distance_angstrom": intended_distance,
                "nearest_non_target_distance_angstrom": nearest_distance,
                "nearest_non_target_ligand_atom": topology[nearest_ligand_index]["name"],
                "nearest_non_target_water_atom": water[nearest_water_index][0].upper(),
                "sanity_collision_below_1p2A": nearest_distance < SANITY_COLLISION_DISTANCE,
            }
            scan_records.append(row)
            if scan_distance == REPRESENTATIVE_DISTANCE:
                pdb_lines = [
                    pdb_atom_line(i + 1, atom["name"], "O6U", "A", 1, ligand[i][1], atom["element"])
                    for i, atom in enumerate(topology)
                ]
                for water_index, (element, xyz) in enumerate(water, start=1):
                    water_name = "OH2" if element.upper() == "O" else f"H{water_index}"
                    pdb_lines.append(pdb_atom_line(len(pdb_lines) + 1, water_name, "TIP", "W", 2, xyz, element.upper()))
                pdb_lines.extend(["TER", "END"])
                pdb_path = pdb_dir / f"{orientation_id}_2p0A.pdb"
                pdb_path.write_text("\n".join(pdb_lines) + "\n", encoding="ascii", newline="\n")
                representative = {**row, "pdb": artifact(pdb_path)}
        if representative is None:
            raise SystemExit(f"Representative geometry missing for {orientation_id}")
        audited.append(
            {
                "orientation_id": orientation_id,
                "source_definition": definition,
                "target_atom": target_name,
                "coordinate_module": artifact(coordinate_path),
                "representative_2p0A": representative,
                "minimum_non_target_distance_over_scan_angstrom": min(
                    row["nearest_non_target_distance_angstrom"] for row in scan_records
                ),
                "sanity_collision_anywhere_in_scan": any(row["sanity_collision_below_1p2A"] for row in scan_records),
                "scan_geometry": scan_records,
                "disposition": "pending_geometry_specific_visual_review",
            }
        )

    result = {
        "schema_version": "1.0",
        "report_type": "o6u_ffparam_water_geometry_audit",
        "status": "pass_geometry_integrity_visual_review_required",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {"generation_report": artifact(report_path), "rtf": artifact(rtf_path)},
        "orientation_count": len(audited),
        "representative_distance_angstrom": REPRESENTATIVE_DISTANCE,
        "distance_integrity_tolerance_angstrom": DISTANCE_TOLERANCE,
        "sanity_collision_distance_angstrom": SANITY_COLLISION_DISTANCE,
        "sanity_collision_note": (
            "This threshold is a file/geometry sanity alarm only. It does not select, exclude, or rank an orientation."
        ),
        "orientations": audited,
        "release_boundary": (
            "All orientations remain pending. Visual inspection and chemical-role adjudication must be recorded; "
            "this audit cannot authorize QM execution."
        ),
    }
    result_path = output_dir / "O6U_FFPARAM_WATER_GEOMETRY_AUDIT.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": result["status"], "orientations": len(audited), "report": str(result_path), "sha256": sha256(result_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
