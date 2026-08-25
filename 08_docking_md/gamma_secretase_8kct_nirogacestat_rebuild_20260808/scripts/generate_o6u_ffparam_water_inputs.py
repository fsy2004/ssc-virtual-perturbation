#!/usr/bin/env python3
"""Generate, but never execute, the frozen O6U FFParam water-interaction inputs."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import runpy
from datetime import datetime, timezone
from pathlib import Path

from ffparam.script_core.CoordinateWriter import crdformat
from ffparam.script_core.Psi4InputWriter import Psi4Writer
from ffparam.script_core.moleculereader import Molecule
from ffparam.script_core.rtftopsf import rtftopsf

from prepare_o6u_crest_input import load_single_sdf, validate_identity
from validate_o6u_crest_ensemble import read_xyz_ensemble


EXPECTED = {
    "source_sdf": "2cb9d769cde4157181a6199b83294cad56cade14ab34a5e86a6deb6790fc28d5",
    "correspondence_tsv": "62b5a9500a0c5e0c2d85eb3fa51fb4e4cb82881dafd3d0b2b38f02231d6935f5",
    "rtf": "1cafc443b2b5a7c2b43a8ae1a5c4c5e6e8d0cfdb81971a0ba0c511a3664c1def",
    "coordinate_template": "9d18e691f3e3afb29d5bb18584a6245b7e1cf3aefb2e660b2bcf08dee6857167",
    "orientation_da": "5ea7e12c464750b8f35d9fa3feed875bb6b69a4091867b8fb6ddf2ad3c6272a3",
    "prescreen": "976211319f88c7c6ef19779d2e0e4ac011eed443ece525d7fb2c030556deb96a",
    "policy": "4e65a9ecb9d90ec9d0c9ee849c9c05f1fd0aae09c7e16a74bb59d8fb5043fccb",
}
EXPECTED_RETAINED = 20
EXPECTED_DISTANCE_GRID = [value / 100 for value in range(150, 301, 5)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def checked(path: Path, label: str, expected_sha256: str | None = None) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SystemExit(f"Missing or empty {label}: {resolved}")
    if expected_sha256 is not None and sha256(resolved) != expected_sha256:
        raise SystemExit(f"{label} hash differs")
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sdf", required=True, type=Path)
    parser.add_argument("--correspondence-tsv", required=True, type=Path)
    parser.add_argument("--rtf", required=True, type=Path)
    parser.add_argument("--coordinate-template", required=True, type=Path)
    parser.add_argument("--geometry-xyz", required=True, type=Path)
    parser.add_argument("--orientation-da", required=True, type=Path)
    parser.add_argument("--prescreen", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--role", choices=("generation_canary", "formal_mp2_target"), required=True)
    parser.add_argument("--ensemble-report", type=Path)
    args = parser.parse_args()

    source = checked(args.source_sdf, "source SDF", EXPECTED["source_sdf"])
    correspondence = checked(args.correspondence_tsv, "atom correspondence", EXPECTED["correspondence_tsv"])
    rtf = checked(args.rtf, "RTF", EXPECTED["rtf"])
    coordinate_template = checked(args.coordinate_template, "coordinate template", EXPECTED["coordinate_template"])
    geometry = checked(args.geometry_xyz, "geometry XYZ")
    orientation_da = checked(args.orientation_da, "orientation plan", EXPECTED["orientation_da"])
    prescreen_path = checked(args.prescreen, "chemical-role prescreen", EXPECTED["prescreen"])
    policy_path = checked(args.policy, "disposition policy", EXPECTED["policy"])
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"Refusing non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    ensemble_path: Path | None = None
    if args.role == "formal_mp2_target":
        if args.ensemble_report is None:
            raise SystemExit("Formal MP2 target generation requires an ensemble report")
        ensemble_path = checked(args.ensemble_report, "MP2 ensemble report")
        ensemble = json.loads(ensemble_path.read_text(encoding="utf-8"))
        if (
            ensemble.get("status") != "pass_five_member_ensemble_independently_reconstructed"
            or ensemble.get("production_approved") is not False
        ):
            raise SystemExit("MP2 ensemble report is not a valid geometry release")
        if ensemble.get("charge_water_target_optimized_xyz", {}).get("sha256") != sha256(geometry):
            raise SystemExit("Geometry XYZ differs from the frozen charge/water target")
    elif args.ensemble_report is not None:
        raise SystemExit("Generation canary must not claim an MP2 ensemble release")

    molecule = load_single_sdf(source)
    identity = validate_identity(molecule)
    expected_elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    frames = read_xyz_ensemble(geometry)
    if len(frames) != 1 or frames[0]["elements"] != expected_elements:
        raise SystemExit("Geometry XYZ atom identity/order differs from frozen O6U")

    import csv

    with correspondence.open("r", encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    if len(mapping) != 76 or len({row["cgenff_atom_name"] for row in mapping}) != 76:
        raise SystemExit("Atom correspondence is not a unique 76-row mapping")
    source_names = [row["cgenff_atom_name"] for row in mapping]

    topology = rtftopsf(topparlist=[str(rtf)], resilist=["o6u"], psf=False).topparPsfs[0]
    if len(topology["atomnames"]) != 76 or len(topology["bondnames"]) != 78:
        raise SystemExit("FFParam topology identity differs")
    template_molecule = Molecule()
    if template_molecule.readcoor(str(coordinate_template), center=False) is None:
        raise SystemExit("FFParam could not read the high-precision CRD template")
    from ffparam.script_core import toppario

    mapped_indices = toppario.createmap(topology["atomnames"])["mapindices"]
    formatted = crdformat(template_molecule.mergepos(), mapped_indices)
    if len(formatted) != 76:
        raise SystemExit("CRD template did not produce 76 mapped coordinates")
    xyz_by_name = {
        name: [float(value) for value in coordinate]
        for name, coordinate in zip(source_names, frames[0]["coordinates"], strict=True)
    }
    formatted_names = [row[1] for row in formatted]
    if set(formatted_names) != set(source_names):
        raise SystemExit("FFParam CRD/topology names differ from the frozen atom mapping")
    for row in formatted:
        row[2:5] = xyz_by_name[row[1]]

    prescreen = json.loads(prescreen_path.read_text(encoding="utf-8"))
    if (
        prescreen.get("status") != "pass_chemical_role_prescreen_visual_review_required"
        or prescreen.get("orientation_count") != 70
        or prescreen.get("retained_for_visual_review_count") != EXPECTED_RETAINED
        or prescreen.get("production_approved") is not False
    ):
        raise SystemExit("Chemical-role prescreen is invalid")
    full_lines = [line.strip() for line in orientation_da.read_text(encoding="utf-8").splitlines() if line.strip()]
    retained_rows = [row for row in prescreen["orientations"] if row.get("prescreen_suggestion") == "retain_for_visual_review"]
    retained_rows.sort(key=lambda row: int(row["source_line_number"]))
    retained_lines = [str(row["source_definition"]) for row in retained_rows]
    if len(retained_lines) != EXPECTED_RETAINED:
        raise SystemExit("Chemical-role prescreen did not retain exactly 20 rows")
    for row, line in zip(retained_rows, retained_lines, strict=True):
        line_number = int(row["source_line_number"])
        if not 1 <= line_number <= len(full_lines) or full_lines[line_number - 1] != line:
            raise SystemExit(f"Prescreen row {row.get('orientation_id')} differs from the frozen DA file")
    subset_da = output_dir / "O6U_WATER_PROBE_POLAR_SUBSET.da"
    subset_da.write_text("\n".join(retained_lines) + "\n", encoding="utf-8", newline="\n")

    qmtop = "nproc=1 mem=2GB lot=HF basis=6-31G(d) keyword={WATER NOSYMM}"
    writer = Psi4Writer(topology, qmtop=qmtop, outpath=str(output_dir), resn="O6U", rescharge=0, multiplicity=1)
    created = writer.run(
        "interaction",
        outpath=str(output_dir),
        coor=formatted,
        da=str(subset_da),
        water=True,
        bsse=False,
        da_atoms=[],
    )
    if not isinstance(created, list) or len(created) != EXPECTED_RETAINED:
        raise SystemExit(f"FFParam created an unexpected orientation count: {created!r}")

    generated_pairs: list[dict[str, object]] = []
    seen_paths: set[Path] = set()
    for retained, pair in zip(retained_rows, created, strict=True):
        if not isinstance(pair, list) or len(pair) != 2:
            raise SystemExit("FFParam returned an invalid input/coordinate pair")
        run_file = checked(Path(pair[0]), "generated Psi4 run file")
        coordinate_file = checked(Path(pair[1]), "generated Psi4 coordinate module")
        if run_file in seen_paths or coordinate_file in seen_paths:
            raise SystemExit("FFParam generated duplicate file paths")
        seen_paths.update((run_file, coordinate_file))
        coordinate_data = runpy.run_path(str(coordinate_file))
        labels = list(coordinate_data.get("intrange", []))
        try:
            distance_grid = [float(label.replace("_", ".")) for label in labels]
        except (AttributeError, ValueError) as exc:
            raise SystemExit(f"Distance scan grid is not numeric for {coordinate_file}") from exc
        if len(distance_grid) != len(EXPECTED_DISTANCE_GRID) or any(
            abs(observed - expected) > 1e-10
            for observed, expected in zip(distance_grid, EXPECTED_DISTANCE_GRID, strict=True)
        ):
            raise SystemExit(f"Distance scan grid differs for {coordinate_file}")
        coordinate_blocks = [key for key in coordinate_data if key.startswith("interaction_")]
        if len(coordinate_blocks) != len(EXPECTED_DISTANCE_GRID):
            raise SystemExit(f"Distance-coordinate count differs for {coordinate_file}")
        run_text = run_file.read_text(encoding="utf-8", errors="replace")
        if "opttheory=\"hf\"" not in run_text.lower() or "optbasis=\"6-31g(d)\"" not in run_text.lower():
            raise SystemExit(f"Generated model chemistry differs in {run_file}")
        if "psi4.energy" not in run_text or "psi4.optimize" in run_text:
            raise SystemExit(f"Generated water input is not a fixed-orientation distance energy scan: {run_file}")
        generated_pairs.append(
            {
                "orientation_id": retained["orientation_id"],
                "source_line_number": retained["source_line_number"],
                "source_definition": retained["source_definition"],
                "run_file": artifact(run_file),
                "coordinate_file": artifact(coordinate_file),
                "distance_grid_angstrom": distance_grid,
            }
        )

    parent_files = sorted(output_dir.glob("monomer*_hf_6-31g(d).py"))
    if len(parent_files) != 2:
        raise SystemExit(f"Expected two FFParam monomer files, found {len(parent_files)}")
    forbidden_outputs = sorted(output_dir.glob("*.out")) + sorted(output_dir.glob("*.log"))
    if forbidden_outputs:
        raise SystemExit("Generation-only stage unexpectedly produced QM outputs")

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_ffparam_water_input_generation",
        "status": "pass_generation_only_visual_review_required",
        "role": args.role,
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "identity": identity,
        "ffparam_distribution_version": importlib.metadata.version("ffparam"),
        "model_chemistry": "HF/6-31G(d) fixed-orientation distance scan, neutral singlet, no BSSE",
        "inputs": {
            "source_sdf": artifact(source),
            "correspondence_tsv": artifact(correspondence),
            "rtf": artifact(rtf),
            "coordinate_template": artifact(coordinate_template),
            "geometry_xyz": artifact(geometry),
            "orientation_da": artifact(orientation_da),
            "prescreen": artifact(prescreen_path),
            "policy": artifact(policy_path),
            "ensemble_report": artifact(ensemble_path) if ensemble_path is not None else None,
        },
        "subset_da": artifact(subset_da),
        "orientation_count": len(generated_pairs),
        "distance_grid_start_angstrom": 1.5,
        "distance_grid_end_angstrom": 3.0,
        "distance_grid_increment_angstrom": 0.05,
        "generated_pairs": generated_pairs,
        "monomer_files": [artifact(path) for path in parent_files],
        "release_boundary": (
            "Generation canary only; archive these inputs and regenerate from the frozen MP2 charge target."
            if args.role == "generation_canary"
            else "Formal input generation only; geometry-specific visual review must freeze the final subset before any QM execution."
        ),
    }
    report_path = output_dir / "O6U_FFPARAM_WATER_INPUT_GENERATION.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "role": args.role, "orientation_count": len(generated_pairs), "report": str(report_path), "report_sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
