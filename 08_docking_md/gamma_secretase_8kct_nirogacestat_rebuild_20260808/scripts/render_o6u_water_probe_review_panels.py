#!/usr/bin/env python3
"""Render hash-bound O6U water-probe geometry panels for direct review."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from rdkit import Chem


COLORS = {
    "C": "#7f8c8d",
    "H": "#d9d9d9",
    "N": "#2f6bff",
    "O": "#e53935",
    "F": "#43a047",
}


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


def verify_recorded_artifact(record: object, label: str) -> Path:
    if not isinstance(record, dict):
        raise RuntimeError(f"Missing artifact record: {label}")
    path = Path(str(record.get("path", ""))).resolve()
    if (
        not path.is_file()
        or path.stat().st_size != record.get("size_bytes")
        or sha256(path) != record.get("sha256")
    ):
        raise RuntimeError(f"Artifact failed size/hash verification: {label}")
    return path


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        atoms.append(
            {
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": int(line[22:26]),
                "coord": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
            }
        )
    return atoms


def deterministic_pca(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = coords.mean(axis=0)
    _, _, vh = np.linalg.svd(coords - center, full_matrices=False)
    rotation = vh.T
    projected = (coords - center) @ rotation
    for column in range(3):
        anchor = int(np.argmax(np.abs(projected[:, column])))
        if projected[anchor, column] < 0:
            rotation[:, column] *= -1
            projected[:, column] *= -1
    if np.linalg.det(rotation) < 0:
        rotation[:, 2] *= -1
    return center, rotation


def draw_panel(
    axis: plt.Axes,
    coords: np.ndarray,
    atoms: list[dict[str, object]],
    bonds: list[tuple[int, int]],
    shown: set[int],
    dims: tuple[int, int],
    target_index: int,
    intended_water_index: int,
    competitor_ligand_index: int,
    competitor_water_index: int,
    title: str,
) -> None:
    xdim, ydim = dims
    for left, right in bonds:
        if left in shown and right in shown:
            axis.plot(
                [coords[left, xdim], coords[right, xdim]],
                [coords[left, ydim], coords[right, ydim]],
                color="#a7a7a7",
                linewidth=0.8,
                zorder=1,
            )
    axis.plot(
        [coords[target_index, xdim], coords[intended_water_index, xdim]],
        [coords[target_index, ydim], coords[intended_water_index, ydim]],
        color="#7b1fa2",
        linewidth=2.0,
        linestyle="--",
        zorder=2,
    )
    axis.plot(
        [coords[competitor_ligand_index, xdim], coords[competitor_water_index, xdim]],
        [coords[competitor_ligand_index, ydim], coords[competitor_water_index, ydim]],
        color="#fb8c00",
        linewidth=1.5,
        linestyle=":",
        zorder=2,
    )
    for index in sorted(shown):
        atom = atoms[index]
        element = str(atom["element"])
        emphasized = index in {target_index, intended_water_index, competitor_ligand_index, competitor_water_index}
        axis.scatter(
            coords[index, xdim],
            coords[index, ydim],
            s=55 if emphasized else (24 if element != "H" else 11),
            c=COLORS.get(element, "#ab47bc"),
            edgecolors="#111111" if emphasized else "none",
            linewidths=0.7,
            zorder=3,
        )
        if emphasized:
            axis.annotate(
                f"{atom['name']}:{atom['resname']}",
                (coords[index, xdim], coords[index, ydim]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
            )
    axis.set_aspect("equal", adjustable="datalim")
    axis.set_title(title, fontsize=8)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.spines[:].set_visible(False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-report", required=True, type=Path)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", required=True, choices=("canary_geometry_review", "formal_geometry_review"))
    args = parser.parse_args()

    template_path = args.template_report.resolve()
    source_sdf = args.source_sdf.resolve()
    output_dir = args.output_dir.resolve()
    if not template_path.is_file() or not source_sdf.is_file():
        raise SystemExit("Template report and source SDF must exist")
    if output_dir.exists():
        raise SystemExit("Output directory already exists; refusing reuse")

    template = load_json(template_path)
    expected_template_role = "canary_template" if args.role == "canary_geometry_review" else "formal_mp2_template"
    if (
        template.get("status") != "pending_visual_adjudication_no_qm_authorized"
        or template.get("role") != expected_template_role
        or template.get("production_approved") is not False
    ):
        raise SystemExit("Template report does not pass its exact role/status gate")
    table_path = verify_recorded_artifact(template.get("adjudication_table"), "template.adjudication_table")
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if len(rows) != 20 or any(row.get("review_decision") != "PENDING" for row in rows):
        raise SystemExit("Expected the exact unsigned 20-row pending adjudication table")

    supplier = Chem.SDMolSupplier(str(source_sdf), removeHs=False, sanitize=False)
    molecule = supplier[0] if len(supplier) else None
    if molecule is None or molecule.GetNumAtoms() != 76:
        raise SystemExit("Source SDF must contain the frozen 76-atom O6U molecule")
    sdf_elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    ligand_bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in molecule.GetBonds()]

    output_dir.mkdir(parents=True, exist_ok=False)
    panel_dir = output_dir / "panels"
    panel_dir.mkdir()
    rendered: list[dict[str, object]] = []
    for row in rows:
        orientation_id = row["orientation_id"]
        pdb_path = Path(row["representative_pdb_path"]).resolve()
        if not pdb_path.is_file() or sha256(pdb_path) != row["representative_pdb_sha256"]:
            raise SystemExit(f"Representative PDB failed hash verification: {orientation_id}")
        atoms = parse_pdb(pdb_path)
        if len(atoms) != 79 or [str(atom["element"]) for atom in atoms[:76]] != sdf_elements:
            raise SystemExit(f"PDB/SDF atom identity differs: {orientation_id}")
        if any(atom["resname"] != "O6U" for atom in atoms[:76]) or any(atom["resname"] != "TIP" for atom in atoms[76:]):
            raise SystemExit(f"PDB ligand/water component identity differs: {orientation_id}")

        names = [str(atom["name"]) for atom in atoms[:76]]
        if len(names) != len(set(names)):
            raise SystemExit(f"Ligand atom names are not unique: {orientation_id}")
        target_index = names.index(row["target_atom"])
        competitor_ligand_index = names.index(row["nearest_non_target_ligand_atom_at_2p0A"])
        water_names = [str(atom["name"]) for atom in atoms[76:]]
        recorded_competitor_water = row["nearest_non_target_water_atom_at_2p0A"]
        competitor_water_candidates = [
            index
            for index in range(76, 79)
            if atoms[index]["name"] == recorded_competitor_water
            or atoms[index]["element"] == recorded_competitor_water
        ]
        if not competitor_water_candidates:
            raise SystemExit(
                f"Recorded competitor water atom cannot be mapped by name or element: {orientation_id}"
            )
        competitor_ligand_coord = np.asarray(atoms[competitor_ligand_index]["coord"])
        competitor_water_index = min(
            competitor_water_candidates,
            key=lambda index: float(
                np.linalg.norm(np.asarray(atoms[index]["coord"]) - competitor_ligand_coord)
            ),
        )
        target_coord = np.asarray(atoms[target_index]["coord"])
        water_indices = list(range(76, 79))
        intended_water_index = min(water_indices, key=lambda index: float(np.linalg.norm(np.asarray(atoms[index]["coord"]) - target_coord)))

        raw_coords = np.vstack([np.asarray(atom["coord"]) for atom in atoms])
        center, rotation = deterministic_pca(raw_coords[:76][np.array(sdf_elements) != "H"])
        coords = (raw_coords - center) @ rotation
        water_bonds = [(76, 77), (77, 78)] if water_names[1] == "OH2" else [(76, 78), (77, 78)]
        bonds = ligand_bonds + water_bonds

        whole = {index for index, atom in enumerate(atoms) if atom["element"] != "H" or index >= 76}
        whole.update({target_index, competitor_ligand_index, competitor_water_index, intended_water_index})
        local = {
            index
            for index in range(79)
            if min(
                np.linalg.norm(raw_coords[index] - raw_coords[target_index]),
                np.linalg.norm(raw_coords[index] - raw_coords[intended_water_index]),
            ) <= 4.0
        }
        local.update({target_index, competitor_ligand_index, competitor_water_index, intended_water_index})

        intended_distance = float(np.linalg.norm(raw_coords[target_index] - raw_coords[intended_water_index]))
        competitor_distance = float(np.linalg.norm(raw_coords[competitor_ligand_index] - raw_coords[competitor_water_index]))
        figure, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=180)
        draw_panel(axes[0], coords, atoms, bonds, whole, (0, 1), target_index, intended_water_index, competitor_ligand_index, competitor_water_index, "Whole molecule PCA XY")
        draw_panel(axes[1], coords, atoms, bonds, local, (0, 1), target_index, intended_water_index, competitor_ligand_index, competitor_water_index, "Local PCA XY (4 A)")
        draw_panel(axes[2], coords, atoms, bonds, local, (0, 2), target_index, intended_water_index, competitor_ligand_index, competitor_water_index, "Local PCA XZ (4 A)")
        figure.suptitle(
            f"{orientation_id} | {row['source_definition']}\n"
            f"purple intended={intended_distance:.3f} A; orange nearest competitor={competitor_distance:.3f} A; "
            f"scan collision annotation={row['sanity_collision_anywhere_in_scan']}",
            fontsize=9,
        )
        figure.tight_layout(rect=(0, 0, 1, 0.88))
        png = panel_dir / f"{orientation_id}.png"
        figure.savefig(png, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        rendered.append(
            {
                "orientation_id": orientation_id,
                "panel": artifact(png),
                "representative_pdb": artifact(pdb_path),
                "target_atom": row["target_atom"],
                "intended_water_atom": atoms[intended_water_index]["name"],
                "intended_distance_angstrom_reconstructed": intended_distance,
                "nearest_competing_ligand_atom": row["nearest_non_target_ligand_atom_at_2p0A"],
                "nearest_competing_water_atom": row["nearest_non_target_water_atom_at_2p0A"],
                "nearest_competing_distance_angstrom_reconstructed": competitor_distance,
                "collision_annotation_only": row["sanity_collision_anywhere_in_scan"].lower() == "true",
            }
        )

    thumbs: list[Image.Image] = []
    width, height = 600, 220
    for record in rendered:
        image = Image.open(record["panel"]["path"]).convert("RGB")
        image.thumbnail((width, height))
        canvas = Image.new("RGB", (width, height + 24), "white")
        canvas.paste(image, ((width - image.width) // 2, 24))
        ImageDraw.Draw(canvas).text((6, 5), str(record["orientation_id"]), fill="black")
        thumbs.append(canvas)
    sheet = Image.new("RGB", (width * 2, (height + 24) * 10), "white")
    for index, thumb in enumerate(thumbs):
        sheet.paste(thumb, ((index % 2) * width, (index // 2) * (height + 24)))
    sheet_path = output_dir / "O6U_WATER_PROBE_REVIEW_CONTACT_SHEET.png"
    sheet.save(sheet_path, dpi=(180, 180))

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_geometry_review_rendering",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_rendering_only_direct_review_required",
        "role": args.role,
        "production_approved": False,
        "template_report": artifact(template_path),
        "adjudication_table": artifact(table_path),
        "source_sdf": artifact(source_sdf),
        "orientation_count": len(rendered),
        "rendering": {
            "engine": "RDKit topology plus Matplotlib 2D projection",
            "orientation": "deterministic heavy-atom PCA",
            "panels": ["whole PCA XY", "local 4 A PCA XY", "local 4 A PCA XZ"],
            "distance_line_policy": "purple intended nearest target-water atom; orange recorded nearest non-target contact",
        },
        "contact_sheet": artifact(sheet_path),
        "orientations": rendered,
        "automatic_decision_applied": False,
        "release_boundary": (
            "These panels support direct geometry-specific review only. They do not decide RUN_QM versus exclusion, "
            "do not authorize QM, and do not replace inspection of the hash-bound PDB coordinates."
        ),
    }
    report_path = output_dir / "O6U_WATER_PROBE_GEOMETRY_REVIEW_RENDERING.json"
    atomic_json = report_path.with_suffix(report_path.suffix + ".tmp")
    atomic_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(atomic_json, report_path)
    print(json.dumps({"status": report["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
