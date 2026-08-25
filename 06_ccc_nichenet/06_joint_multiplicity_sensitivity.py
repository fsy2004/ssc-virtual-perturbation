#!/usr/bin/env python3
"""Joint permutation multiplicity sensitivity for four Notch routes."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASE = (
    ROOT
    / "04_manuscript"
    / "revision_20260722"
    / "communication_spatial_extension_20260729"
)
DEFAULT_LIANA = (
    DEFAULT_BASE
    / "liana_multisample"
    / "formal_per_donor_20260729"
    / "liana_primary_donor_scores.csv"
)
DEFAULT_VISIUM = (
    DEFAULT_BASE
    / "commot"
    / "formal_full_gene_normalized_20260729"
    / "commot_visium_sample_scores.csv"
)
DEFAULT_OUT = DEFAULT_BASE / "joint_multiplicity_sensitivity_20260730"
ROUTES = (
    ("Pericyte", "JAG1", "NOTCH2"),
    ("Pericyte", "JAG1", "NOTCH3"),
    ("Endothelial", "DLL4", "NOTCH2"),
    ("Endothelial", "DLL4", "NOTCH3"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--liana-scores", type=Path, default=DEFAULT_LIANA)
    parser.add_argument("--visium-scores", type=Path, default=DEFAULT_VISIUM)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--liana-permutations", type=int, default=9_999)
    parser.add_argument("--seed", type=int, default=20260730)
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def route_label(source: str, ligand: str, receptor: str) -> str:
    return f"{source}:{ligand}-{receptor}"


def robust_stat(
    score: np.ndarray,
    condition: np.ndarray,
    cohort: np.ndarray | None = None,
) -> tuple[float, float]:
    columns = [np.ones(len(score), dtype=float), condition.astype(float)]
    if cohort is not None:
        dummies = pd.get_dummies(
            pd.Series(cohort, dtype="string"), drop_first=True, dtype=float
        )
        if dummies.shape[1]:
            columns.extend(
                dummies[column].to_numpy(dtype=float) for column in dummies.columns
            )
    design = np.column_stack(columns)
    fit = sm.OLS(score.astype(float), design).fit(cov_type="HC3")
    return float(fit.params[1]), float(fit.tvalues[1])


def adjusted_pvalues(
    observed_t: np.ndarray,
    permuted_t: np.ndarray,
    monte_carlo: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    observed = np.abs(observed_t)
    permuted = np.abs(permuted_t)

    def probability(events: np.ndarray) -> float:
        count = int(np.count_nonzero(events))
        if monte_carlo:
            return (count + 1.0) / (len(events) + 1.0)
        return count / len(events)

    raw = np.asarray(
        [probability(permuted[:, index] >= observed[index]) for index in range(4)]
    )
    row_max = permuted.max(axis=1)
    single_step = np.asarray(
        [probability(row_max >= observed[index]) for index in range(4)]
    )

    order = np.argsort(-observed)
    stepdown_ordered = np.empty(4, dtype=float)
    for rank, route_index in enumerate(order):
        remaining_max = permuted[:, order[rank:]].max(axis=1)
        stepdown_ordered[rank] = probability(
            remaining_max >= observed[route_index]
        )
    stepdown_ordered = np.maximum.accumulate(stepdown_ordered)
    stepdown = np.empty(4, dtype=float)
    stepdown[order] = stepdown_ordered

    positive_raw = np.asarray(
        [probability(permuted_t[:, index] >= observed_t[index]) for index in range(4)]
    )
    positive_row_max = permuted_t.max(axis=1)
    positive_single_step = np.asarray(
        [
            probability(positive_row_max >= observed_t[index])
            for index in range(4)
        ]
    )
    positive_order = np.argsort(-observed_t)
    positive_stepdown_ordered = np.empty(4, dtype=float)
    for rank, route_index in enumerate(positive_order):
        remaining_max = permuted_t[:, positive_order[rank:]].max(axis=1)
        positive_stepdown_ordered[rank] = probability(
            remaining_max >= observed_t[route_index]
        )
    positive_stepdown_ordered = np.maximum.accumulate(
        positive_stepdown_ordered
    )
    positive_stepdown = np.empty(4, dtype=float)
    positive_stepdown[positive_order] = positive_stepdown_ordered
    return (
        raw,
        single_step,
        stepdown,
        positive_raw,
        positive_single_step,
        positive_stepdown,
    )


def global_results(
    observed_t: np.ndarray,
    permuted_t: np.ndarray,
    monte_carlo: bool,
) -> dict[str, float]:
    observed_max = float(np.max(np.abs(observed_t)))
    permuted_max = np.max(np.abs(permuted_t), axis=1)
    observed_mean = float(np.mean(observed_t))
    permuted_mean = np.mean(permuted_t, axis=1)

    def probability(events: np.ndarray) -> float:
        count = int(np.count_nonzero(events))
        if monte_carlo:
            return (count + 1.0) / (len(events) + 1.0)
        return count / len(events)

    return {
        "observed_max_abs_hc3_t": observed_max,
        "global_maxT_two_sided_p": probability(permuted_max >= observed_max),
        "observed_mean_hc3_t": observed_mean,
        "global_concordant_positive_p": probability(
            permuted_mean >= observed_mean
        ),
        "global_mean_t_two_sided_p": probability(
            np.abs(permuted_mean) >= abs(observed_mean)
        ),
    }


def analyse_visium(
    path: Path,
) -> tuple[pd.DataFrame, dict[str, float], int]:
    scores = pd.read_csv(path)
    scores["route"] = [
        route_label(
            "Pericyte" if ligand == "JAG1" else "Endothelial",
            ligand,
            receptor,
        )
        for ligand, receptor in zip(scores["ligand"], scores["receptor"])
    ]
    samples = (
        scores[["sample", "condition"]]
        .drop_duplicates()
        .sort_values("sample")
        .reset_index(drop=True)
    )
    if len(samples) != 14 or int(samples["condition"].eq("HC").sum()) != 4:
        raise ValueError("Visium contract requires 14 samples with four HC")
    sample_position = {sample: index for index, sample in enumerate(samples["sample"])}
    route_data = []
    observed_t = []
    observed_beta = []
    for source, ligand, receptor in ROUTES:
        route = scores[
            scores["ligand"].eq(ligand) & scores["receptor"].eq(receptor)
        ].copy()
        route = route.sort_values("sample")
        positions = route["sample"].map(sample_position).to_numpy(dtype=int)
        condition = samples.loc[positions, "condition"].eq("SSc").to_numpy(dtype=float)
        beta, t_value = robust_stat(
            route["primary_score"].to_numpy(dtype=float), condition
        )
        route_data.append((route, positions))
        observed_beta.append(beta)
        observed_t.append(t_value)

    assignments = list(itertools.combinations(range(len(samples)), 4))
    permuted_t = np.empty((len(assignments), 4), dtype=float)
    for permutation_index, hc_positions in enumerate(assignments):
        permuted_condition = np.ones(len(samples), dtype=float)
        permuted_condition[list(hc_positions)] = 0.0
        for route_index, (route, positions) in enumerate(route_data):
            _, permuted_t[permutation_index, route_index] = robust_stat(
                route["primary_score"].to_numpy(dtype=float),
                permuted_condition[positions],
            )

    observed_t_array = np.asarray(observed_t)
    (
        raw,
        single_step,
        stepdown,
        positive_raw,
        positive_single_step,
        positive_stepdown,
    ) = adjusted_pvalues(
        observed_t_array, permuted_t, monte_carlo=False
    )
    rows = []
    for index, (source, ligand, receptor) in enumerate(ROUTES):
        rows.append(
            {
                "platform": "Visium_COMMOT",
                "source": source,
                "ligand": ligand,
                "receptor": receptor,
                "n_samples": 14,
                "observed_difference": observed_beta[index],
                "observed_hc3_t": observed_t[index],
                "joint_raw_permutation_p": raw[index],
                "westfall_young_single_step_p": single_step[index],
                "westfall_young_stepdown_p": stepdown[index],
                "positive_raw_permutation_p": positive_raw[index],
                "positive_westfall_young_single_step_p": positive_single_step[index],
                "positive_westfall_young_stepdown_p": positive_stepdown[index],
            }
        )
    return (
        pd.DataFrame(rows),
        global_results(observed_t_array, permuted_t, monte_carlo=False),
        len(assignments),
    )


def analyse_liana(
    path: Path,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    scores = pd.read_csv(path)
    scores["route"] = [
        route_label(source, ligand, receptor)
        for source, ligand, receptor in zip(
            scores["source"], scores["ligand"], scores["receptor"]
        )
    ]
    donor_meta = (
        scores[["donor_id", "cohort", "condition"]]
        .drop_duplicates()
        .sort_values("donor_id")
        .reset_index(drop=True)
    )
    if donor_meta["donor_id"].duplicated().any():
        raise ValueError("Donor metadata are not unique")
    donor_position = {
        donor: index for index, donor in enumerate(donor_meta["donor_id"])
    }
    observed_condition = donor_meta["condition"].eq("SSc").to_numpy(dtype=float)
    cohort_values = donor_meta["cohort"].astype(str).to_numpy()
    cohort_positions = [
        np.flatnonzero(cohort_values == cohort)
        for cohort in sorted(donor_meta["cohort"].astype(str).unique())
    ]

    route_data = []
    observed_t = []
    observed_beta = []
    for source, ligand, receptor in ROUTES:
        route = scores[
            scores["source"].eq(source)
            & scores["ligand"].eq(ligand)
            & scores["receptor"].eq(receptor)
        ].copy()
        route = route.sort_values("donor_id")
        positions = route["donor_id"].map(donor_position).to_numpy(dtype=int)
        beta, t_value = robust_stat(
            route["magnitude_score"].to_numpy(dtype=float),
            observed_condition[positions],
            cohort_values[positions],
        )
        route_data.append((route, positions))
        observed_beta.append(beta)
        observed_t.append(t_value)

    rng = np.random.default_rng(seed)
    permuted_t = np.empty((permutations, 4), dtype=float)
    for permutation_index in range(permutations):
        permuted_condition = observed_condition.copy()
        for positions in cohort_positions:
            permuted_condition[positions] = rng.permutation(
                observed_condition[positions]
            )
        for route_index, (route, positions) in enumerate(route_data):
            _, permuted_t[permutation_index, route_index] = robust_stat(
                route["magnitude_score"].to_numpy(dtype=float),
                permuted_condition[positions],
                cohort_values[positions],
            )

    observed_t_array = np.asarray(observed_t)
    (
        raw,
        single_step,
        stepdown,
        positive_raw,
        positive_single_step,
        positive_stepdown,
    ) = adjusted_pvalues(
        observed_t_array, permuted_t, monte_carlo=True
    )
    rows = []
    for index, (source, ligand, receptor) in enumerate(ROUTES):
        rows.append(
            {
                "platform": "LIANA_donor",
                "source": source,
                "ligand": ligand,
                "receptor": receptor,
                "n_donors": len(route_data[index][0]),
                "observed_difference": observed_beta[index],
                "observed_hc3_t": observed_t[index],
                "joint_raw_permutation_p": raw[index],
                "westfall_young_single_step_p": single_step[index],
                "westfall_young_stepdown_p": stepdown[index],
                "positive_raw_permutation_p": positive_raw[index],
                "positive_westfall_young_single_step_p": positive_single_step[index],
                "positive_westfall_young_stepdown_p": positive_stepdown[index],
            }
        )
    return (
        pd.DataFrame(rows),
        global_results(observed_t_array, permuted_t, monte_carlo=True),
    )


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    visium_routes, visium_global, visium_assignments = analyse_visium(
        args.visium_scores
    )
    liana_routes, liana_global = analyse_liana(
        args.liana_scores, args.liana_permutations, args.seed
    )
    route_results = pd.concat([visium_routes, liana_routes], ignore_index=True)
    for alpha in (0.05, 0.10):
        suffix = str(alpha).replace(".", "_")
        route_results[f"raw_two_sided_lt_{suffix}"] = (
            route_results["joint_raw_permutation_p"] < alpha
        )
        route_results[f"stepdown_two_sided_lt_{suffix}"] = (
            route_results["westfall_young_stepdown_p"] < alpha
        )
        route_results[f"raw_positive_lt_{suffix}"] = (
            route_results["positive_raw_permutation_p"] < alpha
        )
        route_results[f"stepdown_positive_lt_{suffix}"] = (
            route_results["positive_westfall_young_stepdown_p"] < alpha
        )
    route_results.to_csv(args.outdir / "joint_route_results.csv", index=False)

    global_table = pd.DataFrame(
        [
            {
                "platform": "Visium_COMMOT",
                "permutations": visium_assignments,
                "permutation_type": "all assignments of four HC labels",
                **visium_global,
            },
            {
                "platform": "LIANA_donor",
                "permutations": args.liana_permutations,
                "permutation_type": "Monte Carlo labels within cohort",
                **liana_global,
            },
        ]
    )
    global_table.to_csv(args.outdir / "joint_global_results.csv", index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_status": (
            "posthoc multiplicity sensitivity; primary BH and Holm results retained"
        ),
        "parameters": {
            "liana_permutations": args.liana_permutations,
            "seed": args.seed,
        },
        "inputs": {
            "liana_scores": {
                "path": str(args.liana_scores.resolve()),
                "bytes": args.liana_scores.stat().st_size,
                "sha256": sha256(args.liana_scores),
            },
            "visium_scores": {
                "path": str(args.visium_scores.resolve()),
                "bytes": args.visium_scores.stat().st_size,
                "sha256": sha256(args.visium_scores),
            },
        },
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "statsmodels": sm.__version__,
        },
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    (args.outdir / "JOINT_MULTIPLICITY_DONE").write_text(
        datetime.now(timezone.utc).isoformat() + "\n", encoding="ascii"
    )
    print(route_results.to_string(index=False))
    print(global_table.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
