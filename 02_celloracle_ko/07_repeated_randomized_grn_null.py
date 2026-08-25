#!/usr/bin/env python3
"""Repeated randomized-GRN audit for prespecified CellOracle perturbations.

One process handles one gene and one GRN arm. The server launcher runs at most
two processes concurrently so the two loaded Oracle objects remain within the
128-GB memory budget. The observed KO and one independently randomized field
are recomputed for every recorded seed. This produces a genuine repeated-null
table; it does not replace the full all-active screen or select parameters from
the HES1 result.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta


DRAW_COLUMNS = [
    "arm",
    "gene",
    "repeat_index",
    "seed",
    "observed_ps",
    "randomized_ps",
    "null_at_least_as_extreme",
    "paired_grid_wilcoxon_p",
    "n_retained_grid_values",
]
FAILURE_COLUMNS = [
    "arm",
    "gene",
    "repeat_index",
    "seed",
    "error_type",
    "error_message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, required=True)
    parser.add_argument("--oracle-pkl", type=Path, required=True)
    parser.add_argument("--arm", choices=("skinatac", "promoter"), required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--repeats", type=int, default=999)
    parser.add_argument("--seed-base", type=int, default=2026072600)
    parser.add_argument("--checkpoint-every", type=int, default=5)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=14)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def load_engine(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location("ssc_celloracle_engine", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import CellOracle engine from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".part")
    frame.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def load_checkpoint(path: Path, columns: list[str], args: argparse.Namespace) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=columns)
    frame = pd.read_table(path)
    missing = set(columns).difference(frame.columns)
    if missing:
        raise RuntimeError(f"{path} lacks checkpoint columns: {sorted(missing)}")
    if not frame.empty:
        if set(frame["arm"].astype(str)) != {args.arm}:
            raise RuntimeError(f"{path} contains a different GRN arm")
        if set(frame["gene"].astype(str)) != {args.gene}:
            raise RuntimeError(f"{path} contains a different gene")
        expected_seeds = set(range(args.seed_base, args.seed_base + args.repeats))
        unexpected = set(frame["seed"].astype(int)).difference(expected_seeds)
        if unexpected:
            raise RuntimeError(
                f"{path} contains seeds outside the requested frozen sequence: "
                f"{sorted(unexpected)[:5]}"
            )
    return frame.loc[:, columns].copy()


def clopper_pearson(exceedances: int, total: int) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    low = 0.0 if exceedances == 0 else float(beta.ppf(0.025, exceedances, total - exceedances + 1))
    high = 1.0 if exceedances == total else float(
        beta.ppf(0.975, exceedances + 1, total - exceedances)
    )
    return low, high


def main() -> None:
    args = parse_args()
    if args.repeats < 99:
        raise ValueError("At least 99 repeated randomized fields are required")
    if args.checkpoint_every < 1:
        raise ValueError("checkpoint-every must be positive")
    if args.n_jobs < 1 or args.threads < 1:
        raise ValueError("--n-jobs and --threads must be positive")
    if not args.oracle_pkl.is_file():
        raise FileNotFoundError(args.oracle_pkl)

    args.outdir.mkdir(parents=True, exist_ok=True)
    draws_path = args.outdir / "repeated_randomized_grn_draws.tsv"
    failure_path = args.outdir / "repeated_randomized_grn_failures.tsv"
    summary_path = args.outdir / "repeated_randomized_grn_summary.json"

    engine = load_engine(args.engine)
    cfg = dict(engine.CFG)
    cfg["seed"] = args.seed_base
    cfg["run_closed_downstream_stages"] = False
    np.random.seed(args.seed_base)

    with args.oracle_pkl.open("rb") as handle:
        oracle = pickle.load(handle)
    if args.gene not in oracle.adata.var_names:
        raise KeyError(f"{args.gene} is absent from {args.oracle_pkl}")

    myo_mask = engine.myo_mask_of(oracle.adata, cfg)
    myo_cell_idx = np.flatnonzero(myo_mask)
    if len(myo_cell_idx) == 0:
        raise RuntimeError("No myofibroblast cells passed the frozen mask")

    gradient = engine.build_dev_gradient(oracle, cfg)
    oracle.simulate_shift(
        perturb_condition={args.gene: 0.0},
        n_propagation=cfg["n_propagation"],
    )
    oracle.calculate_p_mass(
        smooth=0.8,
        n_grid=cfg["vf_n_grid"],
        n_neighbors=cfg["vf_n_neighbors"],
    )
    mass = np.asarray(oracle.total_p_mass)
    min_mass = float(np.quantile(mass[mass > 0], 0.05)) if np.any(mass > 0) else 0.0
    oracle.calculate_mass_filter(min_mass=min_mass, plot=False)

    completed = load_checkpoint(draws_path, DRAW_COLUMNS, args)
    rows = completed.to_dict(orient="records")
    complete_seeds = set(completed.get("seed", pd.Series(dtype=int)).astype(int))

    failures = load_checkpoint(failure_path, FAILURE_COLUMNS, args)
    failure_rows = failures.to_dict(orient="records")
    failed_seeds = set(failures.get("seed", pd.Series(dtype=int)).astype(int))

    from celloracle.applications import Oracle_development_module

    for repeat_index in range(args.repeats):
        seed = args.seed_base + repeat_index
        if seed in complete_seeds or seed in failed_seeds:
            continue
        try:
            np.random.seed(seed)
            oracle.estimate_transition_prob(
                n_neighbors=cfg["vf_n_neighbors"],
                knn_random=True,
                sampled_fraction=cfg["vf_sampled_frac"],
                calculate_randomized=True,
                n_jobs=args.n_jobs,
                threads=args.threads,
                random_seed=seed,
            )
            oracle.calculate_embedding_shift(sigma_corr=cfg["vf_sigma_corr"])
            dev = Oracle_development_module()
            dev.load_differentiation_reference_data(gradient_object=gradient)
            dev.load_perturb_simulation_data(
                oracle_object=oracle,
                cell_idx_use=myo_cell_idx,
                name=f"KO_{args.gene}_seed_{seed}",
            )
            dev.calculate_inner_product()
            dev.calculate_digitized_ip(n_bins=10)
            paired_p, observed_magnitude, null_magnitude = (
                dev.get_negative_PS_p_value(return_ps_sum=True, plot=False)
            )
            observed_ps = -float(observed_magnitude)
            null_ps = -float(null_magnitude)
            if not np.isfinite(observed_ps) or not np.isfinite(null_ps):
                raise ValueError("Non-finite observed or randomized perturbation score")
            rows.append(
                {
                    "arm": args.arm,
                    "gene": args.gene,
                    "repeat_index": repeat_index,
                    "seed": seed,
                    "observed_ps": observed_ps,
                    "randomized_ps": null_ps,
                    "null_at_least_as_extreme": bool(null_ps <= observed_ps),
                    "paired_grid_wilcoxon_p": (
                        float(paired_p) if paired_p is not None else float("nan")
                    ),
                    "n_retained_grid_values": int(len(dev.inner_product_df)),
                }
            )
        except Exception as error:
            failure_rows.append(
                {
                    "arm": args.arm,
                    "gene": args.gene,
                    "repeat_index": repeat_index,
                    "seed": seed,
                    "error_type": type(error).__name__,
                    "error_message": str(error),
                }
            )

        attempted = len(rows) + len(failure_rows)
        if attempted % args.checkpoint_every == 0:
            atomic_tsv(pd.DataFrame(rows, columns=DRAW_COLUMNS), draws_path)
            atomic_tsv(
                pd.DataFrame(failure_rows, columns=FAILURE_COLUMNS), failure_path
            )

    draw_frame = (
        pd.DataFrame(rows, columns=DRAW_COLUMNS)
        .sort_values("seed")
        .drop_duplicates("seed")
    )
    failure_frame = (
        pd.DataFrame(failure_rows, columns=FAILURE_COLUMNS)
        .sort_values("seed")
        .drop_duplicates("seed")
    )
    atomic_tsv(draw_frame, draws_path)
    atomic_tsv(failure_frame, failure_path)

    if len(draw_frame) < max(50, int(0.9 * args.repeats)):
        raise RuntimeError(
            f"Only {len(draw_frame)}/{args.repeats} repeated-null draws succeeded"
        )
    exceedances = int(draw_frame["null_at_least_as_extreme"].sum())
    n_unique_randomized_ps = int(draw_frame["randomized_ps"].nunique(dropna=True))
    if n_unique_randomized_ps < 2:
        raise RuntimeError(
            "Repeated randomized fields collapsed to a single perturbation score; "
            "the CellOracle random seed was not varied effectively"
        )
    empirical_p = float((exceedances + 1) / (len(draw_frame) + 1))
    ci_low, ci_high = clopper_pearson(exceedances, len(draw_frame))
    observed = draw_frame["observed_ps"].to_numpy(dtype=float)
    randomized = draw_frame["randomized_ps"].to_numpy(dtype=float)
    summary = {
        "arm": args.arm,
        "gene": args.gene,
        "oracle_path": str(args.oracle_pkl.resolve()),
        "oracle_size_bytes": args.oracle_pkl.stat().st_size,
        "n_requested": args.repeats,
        "n_successful": int(len(draw_frame)),
        "n_failed": int(len(failure_frame)),
        "n_unique_randomized_ps": n_unique_randomized_ps,
        "seed_base": args.seed_base,
        "celloracle_n_jobs": args.n_jobs,
        "celloracle_threads": args.threads,
        "celloracle_random_seed_explicit": True,
        "alternative": "randomized_ps_less_than_or_equal_to_observed_ps",
        "exceedances": exceedances,
        "empirical_plus_one_p": empirical_p,
        "binomial_95ci_low": ci_low,
        "binomial_95ci_high": ci_high,
        "observed_ps_median": float(np.median(observed)),
        "observed_ps_min": float(np.min(observed)),
        "observed_ps_max": float(np.max(observed)),
        "observed_ps_sd": float(np.std(observed, ddof=1)),
        "randomized_ps_median": float(np.median(randomized)),
        "randomized_ps_q025": float(np.quantile(randomized, 0.025)),
        "randomized_ps_q975": float(np.quantile(randomized, 0.975)),
        "parameter_selection": "outcome_blind_prespecified_gene_and_seed_sequence",
    }
    temporary = summary_path.with_suffix(".json.part")
    temporary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    os.replace(temporary, summary_path)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
