#!/usr/bin/env python3
"""Fit a coordinate-level Fourier diagnostic to closed-scan QM-MM residuals.

This diagnostic selects a small harmonic order by leave-one-out error.  Its
coefficients are coordinate-level residual descriptors, not CHARMM atom-type
parameters and not a promoted force-field candidate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ALLOWED_ROTORS = {"ROT_C09_N04", "ROT_C15_N05", "ROT_C17_C15", "ROT_C24_C14"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def design_matrix(angles_deg: np.ndarray, reference_deg: float, order: int) -> np.ndarray:
    angles = np.deg2rad(np.asarray(angles_deg, dtype=float))
    reference = math.radians(reference_deg)
    columns = []
    for harmonic in range(1, order + 1):
        columns.append(np.cos(harmonic * angles) - math.cos(harmonic * reference))
        columns.append(np.sin(harmonic * angles) - math.sin(harmonic * reference))
    return np.column_stack(columns)


def coefficient_terms(coefficients: list[float] | np.ndarray) -> list[dict[str, float | int]]:
    values = np.asarray(coefficients, dtype=float)
    terms = []
    for index in range(0, len(values), 2):
        cosine = float(values[index])
        sine = float(values[index + 1])
        terms.append({
            "harmonic": index // 2 + 1,
            "cosine_coefficient_kcal_mol": cosine,
            "sine_coefficient_kcal_mol": sine,
            "amplitude_kcal_mol": math.hypot(cosine, sine),
            "phase_deg": math.degrees(math.atan2(sine, cosine)),
        })
    return terms


def fit_order(angles: np.ndarray, residual: np.ndarray, reference: float, order: int) -> dict:
    matrix = design_matrix(angles, reference, order)
    coefficients, _, rank, singular = np.linalg.lstsq(matrix, residual, rcond=None)
    prediction = matrix @ coefficients
    training_rmse = float(np.sqrt(np.mean((prediction - residual) ** 2)))
    loo_errors = []
    for held_out in range(len(residual)):
        keep = np.arange(len(residual)) != held_out
        train_matrix = matrix[keep]
        train_residual = residual[keep]
        if np.linalg.matrix_rank(train_matrix) < matrix.shape[1]:
            loo_errors.append(float("nan"))
            continue
        loo_coefficients, _, _, _ = np.linalg.lstsq(train_matrix, train_residual, rcond=None)
        loo_errors.append(float(matrix[held_out] @ loo_coefficients - residual[held_out]))
    finite_loo = np.asarray([value for value in loo_errors if math.isfinite(value)])
    loo_rmse = float(np.sqrt(np.mean(finite_loo**2))) if len(finite_loo) == len(loo_errors) else None
    condition = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0 else None
    return {
        "order": order,
        "coefficients": coefficients.tolist(),
        "terms": coefficient_terms(coefficients),
        "rank": int(rank),
        "condition_number": condition,
        "training_rmse_kcal_mol": training_rmse,
        "leave_one_out_rmse_kcal_mol": loo_rmse,
        "leave_one_out_errors_kcal_mol": loo_errors,
        "predicted_residual_kcal_mol": prediction.tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consolidated-input", required=True)
    parser.add_argument("--rotor-id", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    if args.rotor_id not in ALLOWED_ROTORS:
        raise ValueError(f"Not an authorized adaptive rotor: {args.rotor_id}")
    input_path = Path(args.consolidated_input).resolve()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    source = json.loads(input_path.read_text())
    if source.get("status") != "pass_consolidated_adaptive_scan_initial_cgenff_fitting_inputs":
        raise ValueError("Consolidated input did not pass")
    points = [row for row in source["points"] if row["rotor_id"] == args.rotor_id]
    if len(points) < 4:
        raise ValueError("At least four residual points are required")
    points.sort(key=lambda row: row["signed_step_index"])
    references = [row["reference_dihedral_deg"] for row in points]
    if max(references) - min(references) > 1e-6:
        raise ValueError("Inconsistent inferred reference dihedral")
    reference = float(np.mean(references))
    angles = np.asarray([row["final_dihedral_deg"] for row in points], dtype=float)
    residual = np.asarray([row["qm_minus_initial_cgenff_delta_kcal_mol"] for row in points], dtype=float)
    max_order = min(2, max(1, (len(points) - 2) // 2))
    candidates = [fit_order(angles, residual, reference, order) for order in range(1, max_order + 1)]
    eligible = [row for row in candidates if row["leave_one_out_rmse_kcal_mol"] is not None]
    if not eligible:
        raise RuntimeError("No full-rank leave-one-out candidate")
    selected = min(eligible, key=lambda row: (row["leave_one_out_rmse_kcal_mol"], row["order"]))
    report = {
        "schema_version": "1.0",
        "report_type": "o6u_coordinate_fourier_residual_diagnostic",
        "status": "pass_coordinate_fourier_residual_diagnostic",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rotor_id": args.rotor_id,
        "input": {"path": str(input_path), "sha256": sha256(input_path), "size_bytes": input_path.stat().st_size},
        "point_count": len(points),
        "reference_dihedral_deg": reference,
        "points": points,
        "candidate_models": candidates,
        "selected_model": selected,
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__},
        "parameter_mutation": False,
        "production_md_approved": False,
        "interpretation_boundary": (
            "Coordinate-level Fourier residual diagnostic only. Coefficients are not directly "
            "transferable CHARMM atom-type torsion parameters and are not a promoted parameter set."
        ),
    }
    output_dir.mkdir(parents=True)
    report_path = output_dir / "O6U_COORDINATE_FOURIER_RESIDUAL_DIAGNOSTIC.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(report_path)
    print(json.dumps({
        "status": report["status"],
        "rotor_id": args.rotor_id,
        "selected_order": selected["order"],
        "leave_one_out_rmse_kcal_mol": selected["leave_one_out_rmse_kcal_mol"],
        "sha256": sha256(report_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
