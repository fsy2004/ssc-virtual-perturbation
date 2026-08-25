#!/usr/bin/env python3
"""Donor-resolved LIANA+ analysis for prespecified perivascular Notch routes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from datetime import datetime, timezone
from importlib.metadata import version as package_version
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ATLAS = ROOT / "server_archive" / "fig_atlas" / "atlas_fig.h5ad"
DEFAULT_DONOR_MAP = (
    ROOT
    / "04_manuscript"
    / "revision_20260722"
    / "donor_metadata"
    / "outputs"
    / "sample_to_donor.csv"
)
DEFAULT_OUT = (
    ROOT
    / "04_manuscript"
    / "revision_20260722"
    / "communication_spatial_extension_20260729"
    / "liana_multisample"
)

PRIMARY_ROUTES = (
    ("Pericyte", "JAG1", "NOTCH2"),
    ("Pericyte", "JAG1", "NOTCH3"),
    ("Endothelial", "DLL4", "NOTCH2"),
    ("Endothelial", "DLL4", "NOTCH3"),
)
RECEIVER = "Myofibroblast"
TARGET_GENES = ("JAG1", "DLL4", "NOTCH2", "NOTCH3")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--donor-map", type=Path, default=DEFAULT_DONOR_MAP)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--expr-prop", type=float, default=0.10)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--liana-permutations", type=int, default=1000)
    parser.add_argument("--condition-permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument(
        "--max-donors",
        type=int,
        default=0,
        help="Outcome-blind donor limit for smoke tests; 0 uses all eligible donors.",
    )
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    return restored


def holm_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.maximum.accumulate(
        ranked * (len(ranked) - np.arange(len(ranked)))
    )
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    return restored


def design_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    cohort = pd.get_dummies(frame["cohort"], prefix="cohort", drop_first=True)
    matrix = pd.concat(
        [
            pd.Series(1.0, index=frame.index, name="const"),
            frame["condition"].eq("SSc").astype(float).rename("condition_SSc"),
            cohort.astype(float),
        ],
        axis=1,
    )
    return matrix


def condition_effect(frame: pd.DataFrame, score_col: str) -> tuple[float, object]:
    matrix = design_matrix(frame)
    fit = sm.OLS(frame[score_col].astype(float), matrix).fit(cov_type="HC3")
    return float(fit.params["condition_SSc"]), fit


def permute_within_cohort(
    frame: pd.DataFrame,
    score_col: str,
    observed: float,
    n_permutations: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    permuted = frame.copy()
    betas = np.empty(n_permutations, dtype=float)
    cohort_positions = [
        np.flatnonzero(frame["cohort"].to_numpy() == cohort)
        for cohort in frame["cohort"].unique()
    ]
    original = frame["condition"].to_numpy(copy=True)
    for index in range(n_permutations):
        labels = original.copy()
        for positions in cohort_positions:
            labels[positions] = rng.permutation(labels[positions])
        permuted["condition"] = labels
        betas[index] = condition_effect(permuted, score_col)[0]
    return float((1 + np.sum(np.abs(betas) >= abs(observed))) / (n_permutations + 1))


def prepare_adata(
    atlas_path: Path,
    donor_map_path: Path,
    min_cells: int,
    max_donors: int,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    donor_map = pd.read_csv(donor_map_path)
    required = {
        "sample",
        "cohort",
        "condition",
        "donor_id",
        "analysis_keep",
        "n_cells_atlas",
    }
    missing = required.difference(donor_map.columns)
    if missing:
        raise ValueError(f"Donor map missing columns: {sorted(missing)}")
    kept = donor_map[donor_map["analysis_keep"].eq(True)].copy()
    if (
        kept["donor_id"].nunique() != 230
        or int(kept["n_cells_atlas"].sum()) != 325636
    ):
        raise ValueError(
            "Retained donor contract failed: expected 230 donors and 325,636 cells"
        )

    backed = ad.read_h5ad(atlas_path, backed="r")
    if backed.shape[0] != 423611:
        raise ValueError(f"Unexpected unfiltered atlas shape: {backed.shape}")
    absent_genes = sorted(set(TARGET_GENES).difference(backed.var_names))
    if absent_genes:
        raise ValueError(f"Atlas is missing primary genes: {absent_genes}")

    metadata = kept.set_index("sample")
    obs = backed.obs.copy()
    retained = obs["sample"].astype(str).isin(metadata.index.astype(str))
    obs = obs.loc[retained].copy()
    if len(obs) != 325636:
        raise ValueError(f"Expected 325,636 retained cells, found {len(obs)}")

    obs["donor_id"] = obs["sample"].astype(str).map(metadata["donor_id"])
    obs["cohort"] = obs["sample"].astype(str).map(metadata["cohort"])
    obs["condition"] = obs["sample"].astype(str).map(metadata["condition"])
    obs["ct_work"] = obs["celltype"].astype(str)
    fib = obs["fib_subtype"].astype(str)
    assigned = obs["celltype"].astype(str).eq("Fibroblast") & ~fib.isin(
        ["NA", "nan", "None", ""]
    )
    obs.loc[assigned, "ct_work"] = fib.loc[assigned]

    counts = (
        obs.groupby(["donor_id", "cohort", "condition", "ct_work"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for cell_type in ("Pericyte", "Endothelial", RECEIVER):
        if cell_type not in counts:
            counts[cell_type] = 0
    eligible_any = (
        (counts[RECEIVER] >= min_cells)
        & (
            (counts["Pericyte"] >= min_cells)
            | (counts["Endothelial"] >= min_cells)
        )
    )
    eligible_donors = counts.index.get_level_values("donor_id")[eligible_any]
    if max_donors:
        eligible_donors = np.sort(np.asarray(eligible_donors))[:max_donors]
    relevant = (
        retained
        & backed.obs_names.isin(
            obs.index[
                obs["donor_id"].isin(set(eligible_donors))
                & obs["ct_work"].isin(("Pericyte", "Endothelial", RECEIVER))
            ]
        )
    )
    adata = backed[relevant, list(TARGET_GENES)].to_memory()
    selected_obs = obs.loc[adata.obs_names]
    adata.obs["donor_id"] = selected_obs["donor_id"].astype(str).to_numpy()
    adata.obs["cohort"] = selected_obs["cohort"].astype(str).to_numpy()
    adata.obs["condition"] = selected_obs["condition"].astype(str).to_numpy()
    adata.obs["ct_work"] = pd.Categorical(selected_obs["ct_work"].astype(str))

    route_eligibility = counts.reset_index()
    route_rows = []
    for source, ligand, receptor in PRIMARY_ROUTES:
        for row in route_eligibility.itertuples(index=False):
            route_rows.append(
                {
                    "donor_id": row.donor_id,
                    "cohort": row.cohort,
                    "condition": row.condition,
                    "source": source,
                    "target": RECEIVER,
                    "ligand": ligand,
                    "receptor": receptor,
                    "n_sender": int(getattr(row, source)),
                    "n_receiver": int(getattr(row, RECEIVER)),
                    "eligible": bool(
                        getattr(row, source) >= min_cells
                        and getattr(row, RECEIVER) >= min_cells
                    ),
                }
            )
    route_eligibility = pd.DataFrame(route_rows)
    if max_donors:
        route_eligibility = route_eligibility[
            route_eligibility["donor_id"].isin(set(eligible_donors))
        ].copy()
    return adata, kept, route_eligibility


def callable_donor_audit(
    adata: ad.AnnData, expr_prop: float
) -> pd.DataFrame:
    gene_position = {gene: adata.var_names.get_loc(gene) for gene in TARGET_GENES}
    rows = []
    donor_values = adata.obs["donor_id"].astype(str).to_numpy()
    cell_type_values = adata.obs["ct_work"].astype(str).to_numpy()
    for donor in sorted(np.unique(donor_values)):
        fractions: dict[tuple[str, str], float] = {}
        for cell_type in ("Pericyte", "Endothelial", RECEIVER):
            positions = np.flatnonzero(
                (donor_values == donor) & (cell_type_values == cell_type)
            )
            for gene in TARGET_GENES:
                if len(positions):
                    column = adata.X[positions, gene_position[gene]]
                    nonzero = column.nnz if sp.issparse(column) else np.count_nonzero(column)
                    value = nonzero / len(positions)
                else:
                    value = 0.0
                fractions[(cell_type, gene)] = float(value)
        route_callable = []
        for source, ligand, receptor in PRIMARY_ROUTES:
            route_callable.append(
                fractions[(source, ligand)] >= expr_prop
                and fractions[(RECEIVER, receptor)] >= expr_prop
            )
        row = {
            "donor_id": donor,
            "liana_callable": bool(any(route_callable)),
            "n_callable_primary_routes": int(sum(route_callable)),
        }
        for (cell_type, gene), value in fractions.items():
            row[f"fraction_{cell_type}_{gene}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


def run_liana(
    adata: ad.AnnData, args: argparse.Namespace
) -> tuple[pd.DataFrame, pd.DataFrame]:
    interactions = [(ligand, receptor) for _, ligand, receptor in PRIMARY_ROUTES]
    callability = callable_donor_audit(adata, args.expr_prop)
    callable_donors = set(
        callability.loc[callability["liana_callable"], "donor_id"].astype(str)
    )
    if not callable_donors:
        raise RuntimeError("No donor has a callable primary LIANA+ interaction")
    donor_values = adata.obs["donor_id"].astype(str)
    result_frames = []
    statuses = []
    empty_error = "cannot set a frame with no defined index and a scalar"
    for donor_index, donor in enumerate(sorted(callable_donors), start=1):
        donor_adata = adata[donor_values.eq(donor)].copy()
        try:
            donor_result = li.mt.rank_aggregate(
                donor_adata,
                groupby="ct_work",
                resource_name="consensus",
                interactions=interactions,
                expr_prop=args.expr_prop,
                min_cells=args.min_cells,
                use_raw=False,
                n_perms=args.liana_permutations,
                seed=args.seed,
                n_jobs=args.n_jobs,
                inplace=False,
                verbose=False,
            )
        except ValueError as error:
            if empty_error not in str(error):
                raise
            donor_result = None
            status = "empty_after_liana_internal_filter"
        else:
            status = (
                "returned"
                if donor_result is not None and not donor_result.empty
                else "empty"
            )
        n_rows = 0 if donor_result is None else len(donor_result)
        statuses.append(
            {"donor_id": donor, "liana_status": status, "n_liana_rows": n_rows}
        )
        if donor_result is not None and not donor_result.empty:
            donor_result = donor_result.copy()
            donor_result["donor_id"] = donor
            result_frames.append(donor_result)
        print(
            f"[LIANA+] donor {donor_index}/{len(callable_donors)} "
            f"{donor}: {status} ({n_rows} rows)",
            flush=True,
        )
    callability = callability.merge(
        pd.DataFrame(statuses), on="donor_id", how="left", validate="one_to_one"
    )
    callability["liana_status"] = callability["liana_status"].fillna(
        "not_callable_by_expression"
    )
    callability["n_liana_rows"] = callability["n_liana_rows"].fillna(0).astype(int)
    if not result_frames:
        raise RuntimeError("LIANA+ returned no donor-level interactions")
    result = pd.concat(result_frames, ignore_index=True)
    return result, callability


def assemble_scores(
    liana_result: pd.DataFrame, eligibility: pd.DataFrame
) -> pd.DataFrame:
    donor_column = "donor_id" if "donor_id" in liana_result else "sample"
    required = {
        donor_column,
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "magnitude_rank",
    }
    missing = required.difference(liana_result.columns)
    if missing:
        raise ValueError(f"LIANA+ result missing columns: {sorted(missing)}")
    available = liana_result.rename(
        columns={
            donor_column: "donor_id",
            "ligand_complex": "ligand",
            "receptor_complex": "receptor",
        }
    )
    available = available[
        [
            "donor_id",
            "source",
            "target",
            "ligand",
            "receptor",
            "magnitude_rank",
            "specificity_rank",
        ]
    ].copy()
    available["donor_id"] = available["donor_id"].astype(str)
    scores = eligibility[eligibility["eligible"]].merge(
        available,
        how="left",
        on=["donor_id", "source", "target", "ligand", "receptor"],
        validate="one_to_one",
    )
    scores["returned_by_liana"] = scores["magnitude_rank"].notna()
    scores["magnitude_rank"] = scores["magnitude_rank"].fillna(1.0)
    scores["specificity_rank"] = scores["specificity_rank"].fillna(1.0)
    scores["magnitude_score"] = -np.log10(
        np.maximum(scores["magnitude_rank"].astype(float), 1e-12)
    )
    return scores


def fit_primary(
    scores: pd.DataFrame, n_permutations: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    results = []
    loo_rows = []
    for route_index, (source, ligand, receptor) in enumerate(PRIMARY_ROUTES):
        route = scores[
            scores["source"].eq(source)
            & scores["ligand"].eq(ligand)
            & scores["receptor"].eq(receptor)
        ].copy()
        observed, fit = condition_effect(route, "magnitude_score")
        ci = fit.conf_int().loc["condition_SSc"]
        permutation_p = permute_within_cohort(
            route,
            "magnitude_score",
            observed,
            n_permutations,
            seed + route_index,
        )
        results.append(
            {
                "source": source,
                "target": RECEIVER,
                "ligand": ligand,
                "receptor": receptor,
                "n_donors": len(route),
                "n_SSc": int(route["condition"].eq("SSc").sum()),
                "n_HC": int(route["condition"].eq("HC").sum()),
                "median_score_SSc": float(
                    route.loc[route["condition"].eq("SSc"), "magnitude_score"].median()
                ),
                "median_score_HC": float(
                    route.loc[route["condition"].eq("HC"), "magnitude_score"].median()
                ),
                "adjusted_difference": observed,
                "ci_low": float(ci.iloc[0]),
                "ci_high": float(ci.iloc[1]),
                "hc3_p": float(fit.pvalues["condition_SSc"]),
                "permutation_p": permutation_p,
                "returned_fraction": float(route["returned_by_liana"].mean()),
            }
        )
        for cohort in sorted(route["cohort"].unique()):
            subset = route[~route["cohort"].eq(cohort)].copy()
            beta, _ = condition_effect(subset, "magnitude_score")
            loo_rows.append(
                {
                    "source": source,
                    "ligand": ligand,
                    "receptor": receptor,
                    "left_out_cohort": cohort,
                    "n_donors": len(subset),
                    "adjusted_difference": beta,
                }
            )

    result = pd.DataFrame(results)
    result["BH_q"] = bh_adjust(result["permutation_p"])
    result["Holm_p"] = holm_adjust(result["permutation_p"])
    return result, pd.DataFrame(loo_rows)


def main() -> int:
    args = parse_args()
    started = time.time()
    outdir = args.outdir
    if args.max_donors:
        outdir = outdir / f"smoke_{args.max_donors}_donors"
    outdir.mkdir(parents=True, exist_ok=True)

    adata, kept, eligibility = prepare_adata(
        args.atlas, args.donor_map, args.min_cells, args.max_donors
    )
    eligibility.to_csv(outdir / "liana_route_eligibility.csv", index=False)
    liana_result, callability = run_liana(adata, args)
    callability.to_csv(outdir / "liana_donor_callability.csv", index=False)
    liana_result.to_csv(outdir / "liana_by_donor_raw.csv", index=False)
    scores = assemble_scores(liana_result, eligibility)
    scores.to_csv(outdir / "liana_primary_donor_scores.csv", index=False)
    primary, loo = fit_primary(
        scores, args.condition_permutations, args.seed
    )
    primary.to_csv(outdir / "liana_primary_condition_results.csv", index=False)
    loo.to_csv(outdir / "liana_primary_leave_one_cohort_out.csv", index=False)

    audit = {
        "atlas_cells_unfiltered": 423611,
        "atlas_cells_retained": 325636,
        "retained_donors": 230,
        "retained_SSc_donors": int(
            kept.loc[kept["condition"].eq("SSc"), "donor_id"].nunique()
        ),
        "retained_HC_donors": int(
            kept.loc[kept["condition"].eq("HC"), "donor_id"].nunique()
        ),
        "analysis_cells": int(adata.n_obs),
        "analysis_donors": int(adata.obs["donor_id"].nunique()),
        "liana_callable_donors": int(callability["liana_callable"].sum()),
        "eligible_donor_routes": int(eligibility["eligible"].sum()),
    }
    (outdir / "liana_input_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command_parameters": vars(args) | {"outdir": str(outdir)},
        "inputs": {
            "atlas": {
                "path": str(args.atlas.resolve()),
                "bytes": args.atlas.stat().st_size,
                "sha256": sha256(args.atlas),
            },
            "donor_map": {
                "path": str(args.donor_map.resolve()),
                "bytes": args.donor_map.stat().st_size,
                "sha256": sha256(args.donor_map),
            },
        },
        "primary_routes": [
            {
                "source": source,
                "ligand": ligand,
                "receptor": receptor,
                "target": RECEIVER,
            }
            for source, ligand, receptor in PRIMARY_ROUTES
        ],
        "audit": audit,
        "runtime_seconds": time.time() - started,
        "versions": {
            "python": platform.python_version(),
            "anndata": package_version("anndata"),
            "liana": package_version("liana"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": package_version("statsmodels"),
        },
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )
    if not args.max_donors:
        (outdir / "LIANA_MULTISAMPLE_DONE").write_text(
            datetime.now(timezone.utc).isoformat() + "\n", encoding="ascii"
        )
    print(primary.to_string(index=False))
    print(f"LIANA+ donor analysis complete: {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
