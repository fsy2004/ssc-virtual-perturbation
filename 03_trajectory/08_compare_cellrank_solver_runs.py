#!/usr/bin/env python3
"""Compare two complete CellRank runs for numerical solver stability."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ARMS = (
    "baseline",
    "remove_hes1",
    "remove_collectri_hes1_targets",
)
METADATA_COLUMNS = (
    "cell_id",
    "arm",
    "donor_id",
    "cohort",
    "condition",
    "sample",
    "fib_subtype",
)
NUMERIC_COLUMNS = ("pseudotime", "myo_fate")
CELL_TOLERANCE = 1e-6
SUMMARY_TOLERANCE = 1e-6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--candidate-dir", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare_cells(reference_dir: Path, candidate_dir: Path) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        name = f"Figure3_{arm}_per_cell.csv.gz"
        reference = pd.read_csv(reference_dir / name)
        candidate = pd.read_csv(candidate_dir / name)
        if reference.shape != candidate.shape:
            raise RuntimeError(
                f"{arm}: shape differs: {reference.shape} versus {candidate.shape}"
            )
        metadata_equal = all(
            reference[column].equals(candidate[column])
            for column in METADATA_COLUMNS
        )
        if not metadata_equal:
            raise RuntimeError(f"{arm}: cell order or metadata differs")
        for column in NUMERIC_COLUMNS:
            difference = np.abs(
                reference[column].to_numpy(float)
                - candidate[column].to_numpy(float)
            )
            rows.append(
                {
                    "arm": arm,
                    "column": column,
                    "n_cells": len(reference),
                    "metadata_equal": metadata_equal,
                    "maximum_absolute_difference": float(difference.max()),
                    "mean_absolute_difference": float(difference.mean()),
                    "tolerance": CELL_TOLERANCE,
                    "passes_tolerance": bool(
                        difference.max() <= CELL_TOLERANCE
                    ),
                    "reference_sha256": sha256(reference_dir / name),
                    "candidate_sha256": sha256(candidate_dir / name),
                }
            )
    return pd.DataFrame(rows)


def compare_table(
    reference_dir: Path,
    candidate_dir: Path,
    filename: str,
    keys: list[str],
    values: list[str],
) -> pd.DataFrame:
    reference = pd.read_csv(reference_dir / filename)
    candidate = pd.read_csv(candidate_dir / filename)
    merged = reference.merge(
        candidate,
        on=keys,
        how="outer",
        suffixes=("_reference", "_candidate"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        raise RuntimeError(f"{filename}: comparison keys differ")
    rows = []
    for column in values:
        difference = np.abs(
            merged[f"{column}_reference"].to_numpy(float)
            - merged[f"{column}_candidate"].to_numpy(float)
        )
        rows.append(
            {
                "table": filename,
                "column": column,
                "n_rows": len(merged),
                "maximum_absolute_difference": float(difference.max()),
                "mean_absolute_difference": float(difference.mean()),
                "tolerance": SUMMARY_TOLERANCE,
                "passes_tolerance": bool(
                    difference.max() <= SUMMARY_TOLERANCE
                ),
            }
        )
    return pd.DataFrame(rows)


def self_test() -> None:
    values = np.array([0.0, 1e-8, 2e-8])
    assert float(values.max()) <= CELL_TOLERANCE
    print("self-test passed")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if args.reference_dir is None or args.candidate_dir is None or args.outdir is None:
        raise SystemExit(
            "--reference-dir, --candidate-dir and --outdir are required"
        )
    args.outdir.mkdir(parents=True, exist_ok=True)
    cell_comparison = compare_cells(args.reference_dir, args.candidate_dir)
    summary_comparison = compare_table(
        args.reference_dir,
        args.candidate_dir,
        "Figure3_leave_feature_out_summary.csv",
        ["arm", "activity"],
        [
            "median_rho",
            "bootstrap_95_ci_low",
            "bootstrap_95_ci_high",
            "positive_fraction",
            "wilcoxon_p",
            "bh_q_within_arm",
            "holm_p_within_arm",
        ],
    )
    paired_comparison = compare_table(
        args.reference_dir,
        args.candidate_dir,
        "Figure3_leave_feature_out_paired_comparison.csv",
        ["comparison_arm", "activity"],
        [
            "median_delta_rho",
            "bootstrap_95_ci_low",
            "bootstrap_95_ci_high",
            "paired_wilcoxon_p",
            "bh_q_within_comparison",
            "holm_p_within_comparison",
        ],
    )
    tables = pd.concat(
        [summary_comparison, paired_comparison],
        ignore_index=True,
    )
    cell_path = args.outdir / "Figure3_solver_cell_comparison.csv"
    table_path = args.outdir / "Figure3_solver_summary_comparison.csv"
    cell_comparison.to_csv(cell_path, index=False)
    tables.to_csv(table_path, index=False)
    passed = bool(
        cell_comparison["passes_tolerance"].all()
        and tables["passes_tolerance"].all()
    )
    manifest = {
        "purpose": "CellRank PETSc execution-mode numerical stability",
        "reference_dir": str(args.reference_dir),
        "candidate_dir": str(args.candidate_dir),
        "cell_tolerance": CELL_TOLERANCE,
        "summary_tolerance": SUMMARY_TOLERANCE,
        "all_checks_passed": passed,
        "reference_manifest_sha256": sha256(
            args.reference_dir / "Figure3_leave_feature_out_manifest.json"
        ),
        "candidate_manifest_sha256": sha256(
            args.candidate_dir / "Figure3_leave_feature_out_manifest.json"
        ),
        "script_sha256": sha256(Path(__file__).resolve()),
        "outputs": {
            str(cell_path): sha256(cell_path),
            str(table_path): sha256(table_path),
        },
    }
    (args.outdir / "Figure3_solver_stability_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    if not passed:
        raise SystemExit("solver-run comparison exceeded numerical tolerance")
    print(f"completed: {args.outdir}")


if __name__ == "__main__":
    main()
