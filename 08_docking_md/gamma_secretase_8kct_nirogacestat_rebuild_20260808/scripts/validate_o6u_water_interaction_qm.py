#!/usr/bin/env python3
"""Independently reconstruct O6U HF/6-31G(d) water-interaction curves."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.5094740631
BATCH_STATUS = "pass_raw_water_qm_execution_outputs_present_validation_required"
GRID = [round(1.5 + 0.05 * index, 2) for index in range(31)]


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


def parse_monomer_energy(text: str) -> float:
    matches = re.findall(r"INTERACTION MONOMER ENERGY is\s*:\s*([-+0-9.eE]+)", text)
    if len(matches) != 1:
        raise RuntimeError("Expected exactly one monomer-energy marker")
    value = float(matches[0])
    if not math.isfinite(value):
        raise RuntimeError("Monomer energy is not finite")
    return value


def parse_interaction_table(text: str) -> tuple[list[tuple[float, float]], float, float]:
    start_marker = "INTERACTION TABLE NOBSSE START"
    end_marker = "INTERACTION TABLE NOBSSE END"
    if text.count(start_marker) != 1 or text.count(end_marker) != 1:
        raise RuntimeError("Expected exactly one interaction-table marker pair")
    section = text.split(start_marker, 1)[1].split(end_marker, 1)[0]
    rows: list[tuple[float, float]] = []
    for line in section.splitlines():
        fields = line.strip().split()
        if len(fields) != 2:
            continue
        try:
            distance, energy = map(float, fields)
        except ValueError:
            continue
        if not math.isfinite(distance) or not math.isfinite(energy):
            raise RuntimeError("Interaction table contains a non-finite value")
        rows.append((distance, energy))
    if len(rows) != 31 or [round(distance, 2) for distance, _ in rows] != GRID:
        raise RuntimeError("Interaction table does not contain the exact ordered 31-point grid")
    if len({distance for distance, _ in rows}) != 31:
        raise RuntimeError("Interaction table contains duplicate distances")
    reported = re.findall(r"INTERACTION DISTANCE and ENERGY are:\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)", text)
    if len(reported) != 1:
        raise RuntimeError("Expected exactly one reported interaction minimum")
    reported_distance, reported_energy = map(float, reported[0])
    minimum_distance, minimum_energy = min(rows, key=lambda item: item[1])
    if abs(reported_distance - minimum_distance) > 1e-12 or abs(reported_energy - minimum_energy) > 1e-12:
        raise RuntimeError("Reported interaction minimum differs from independent table reconstruction")
    return rows, minimum_distance, minimum_energy


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def atomic_json(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def self_test() -> int:
    monomer = "header\nINTERACTION MONOMER ENERGY is : -76.000000000000\n"
    rows = [(distance, -100.0 - (distance - 2.1) ** 2) for distance in GRID]
    minimum_distance, minimum_energy = min(rows, key=lambda item: item[1])
    text = (
        f"INTERACTION DISTANCE and ENERGY are: {minimum_distance} {minimum_energy}\n"
        "INTERACTION TABLE NOBSSE START\n"
        + "\n".join(f"{distance} {energy}" for distance, energy in rows)
        + "\nINTERACTION TABLE NOBSSE END\n"
    )
    assert parse_monomer_energy(monomer) == -76.0
    parsed, parsed_distance, parsed_energy = parse_interaction_table(text)
    assert parsed == rows and parsed_distance == minimum_distance and parsed_energy == minimum_energy
    print(json.dumps({"status": "pass_synthetic_parser_self_test_no_scientific_data"}, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-report", type=Path)
    parser.add_argument("--generation-report", type=Path)
    parser.add_argument("--authorization-report", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        if any(value is not None for value in (args.batch_report, args.generation_report, args.authorization_report, args.output_dir)):
            raise SystemExit("--self-test cannot be combined with file inputs")
        return self_test()
    if any(value is None for value in (args.batch_report, args.generation_report, args.authorization_report, args.output_dir)):
        raise SystemExit("Formal validation requires all four file/directory arguments")

    batch_path = args.batch_report.resolve()
    generation_path = args.generation_report.resolve()
    authorization_path = args.authorization_report.resolve()
    output_dir = args.output_dir.resolve()
    if not batch_path.is_file() or not generation_path.is_file() or not authorization_path.is_file():
        raise SystemExit("Batch, generation, and authorization reports must exist")
    if output_dir.exists():
        raise SystemExit("Output directory already exists; refusing reuse")
    output_dir.mkdir(parents=True, exist_ok=False)

    batch = load_json(batch_path)
    generation = load_json(generation_path)
    authorization = load_json(authorization_path)
    if (
        batch.get("status") != BATCH_STATUS
        or batch.get("role") != "formal_execution"
        or batch.get("water_interaction_qm_executed") is not True
        or batch.get("parameter_fitting_authorized") is not False
        or batch.get("production_approved") is not False
    ):
        raise SystemExit("Batch report does not pass its exact formal raw-output gate")
    if verify_record(batch.get("generation_report"), "batch.generation_report") != generation_path:
        raise SystemExit("Batch is bound to a different generation report")
    if verify_record(batch.get("authorization_report"), "batch.authorization_report") != authorization_path:
        raise SystemExit("Batch is bound to a different authorization report")
    selected_ids = authorization.get("run_qm_orientation_ids")
    if not isinstance(selected_ids, list) or not selected_ids or selected_ids != batch.get("selected_orientation_ids"):
        raise SystemExit("Batch selected orientations differ from the frozen authorization")

    jobs = batch.get("jobs")
    if not isinstance(jobs, list) or len(jobs) != len(selected_ids) + 2:
        raise SystemExit("Batch job universe is incomplete")
    job_by_label: dict[str, dict[str, object]] = {}
    for job in jobs:
        if not isinstance(job, dict):
            raise SystemExit("Batch job record is not an object")
        label = str(job.get("label", ""))
        if not label or label in job_by_label or job.get("returncode") != 0:
            raise SystemExit("Batch contains a blank, duplicate, or failed job")
        job_by_label[label] = job
    if set(job_by_label) != set(selected_ids) | {"monomer1", "monomer2"}:
        raise SystemExit("Batch job labels differ from the frozen selected universe plus two monomers")

    monomer_energies: dict[str, float] = {}
    for label in ("monomer1", "monomer2"):
        raw = verify_record(job_by_label[label].get("raw_output"), f"{label}.raw_output")
        monomer_energies[label] = parse_monomer_energy(raw.read_text(encoding="utf-8", errors="replace"))
    monomer_sum = monomer_energies["monomer1"] + monomer_energies["monomer2"]

    curve_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    raw_outputs: list[dict[str, object]] = []
    for orientation_id in selected_ids:
        raw = verify_record(job_by_label[orientation_id].get("raw_output"), f"{orientation_id}.raw_output")
        table, minimum_distance, minimum_complex_energy = parse_interaction_table(raw.read_text(encoding="utf-8", errors="replace"))
        raw_outputs.append({"orientation_id": orientation_id, "raw_output": artifact(raw)})
        for distance_value, complex_energy in table:
            curve_rows.append(
                {
                    "orientation_id": orientation_id,
                    "distance_angstrom": distance_value,
                    "complex_energy_hartree": complex_energy,
                    "monomer_sum_hartree": monomer_sum,
                    "interaction_energy_hartree": complex_energy - monomer_sum,
                    "interaction_energy_kcal_mol": (complex_energy - monomer_sum) * HARTREE_TO_KCAL_MOL,
                }
            )
        summary_rows.append(
            {
                "orientation_id": orientation_id,
                "minimum_distance_angstrom": minimum_distance,
                "minimum_complex_energy_hartree": minimum_complex_energy,
                "minimum_interaction_energy_hartree": minimum_complex_energy - monomer_sum,
                "minimum_interaction_energy_kcal_mol": (minimum_complex_energy - monomer_sum) * HARTREE_TO_KCAL_MOL,
            }
        )

    curve_path = output_dir / "O6U_WATER_INTERACTION_CURVES_RAW.tsv"
    summary_path = output_dir / "O6U_WATER_INTERACTION_MINIMA_RAW.tsv"
    write_tsv(curve_path, list(curve_rows[0]), curve_rows)
    write_tsv(summary_path, list(summary_rows[0]), summary_rows)
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_interaction_qm_independent_numerical_validation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass_raw_water_qm_independent_numerical_reconstruction",
        "production_approved": False,
        "parameter_fitting_authorized": False,
        "batch_report": artifact(batch_path),
        "generation_report": artifact(generation_path),
        "authorization_report": artifact(authorization_path),
        "monomer_energies_hartree": monomer_energies,
        "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        "selected_orientation_count": len(selected_ids),
        "selected_orientation_ids": selected_ids,
        "curve_point_count": len(curve_rows),
        "curve_table": artifact(curve_path),
        "minimum_table": artifact(summary_path),
        "raw_outputs": raw_outputs,
        "automatic_scientific_classification_applied": False,
        "release_boundary": (
            "This report reconstructs raw no-BSSE HF/6-31G(d) curves and interaction energies only. Every selected "
            "orientation still requires a signed post-QM disposition (applicable, weak, or unfavourable), failed "
            "attempts remain visible, and no parameter fitting or downstream structure/MD is authorized."
        ),
    }
    report_path = output_dir / "O6U_WATER_INTERACTION_QM_INDEPENDENT_VALIDATION.json"
    atomic_json(report_path, report)
    print(json.dumps({"status": report["status"], "report": str(report_path), "sha256": sha256(report_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
