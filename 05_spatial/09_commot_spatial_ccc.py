#!/usr/bin/env python3
"""Sample-resolved COMMOT analysis of perivascular Notch communication."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import time
from datetime import datetime, timezone
from importlib.metadata import version as package_version
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# COMMOT 0.0.3 still uses the removed NumPy 1.x alias.
np.Inf = np.inf

import commot as ct
import pandas as pd
import scanpy as sc
import scipy
import scipy.sparse as sp
import seaborn as sns
import statsmodels.api as sm
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VISIUM = (
    ROOT
    / "server_archive"
    / "tier2_large"
    / "spatial_extract"
    / "raw"
    / "spatial"
    / "zenodo_14577696"
    / "visium_all.h5ad"
)
DEFAULT_XENIUM = (
    ROOT
    / ".codex_tmp"
    / "xenium_raw"
    / "raw"
    / "spatial"
    / "GSE312932"
    / "xenium_all.h5ad"
)
DEFAULT_RCTD = (
    ROOT / "server_archive" / "tier1" / "powered" / "spatial" / "rctd"
)
DEFAULT_OUT = (
    ROOT
    / "04_manuscript"
    / "revision_20260722"
    / "communication_spatial_extension_20260729"
    / "commot"
)
DATABASE_NAME = "notch4"
PRIMARY_ROUTES = (
    ("JAG1", "NOTCH2"),
    ("JAG1", "NOTCH3"),
    ("DLL4", "NOTCH2"),
    ("DLL4", "NOTCH3"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visium-h5ad", type=Path, default=DEFAULT_VISIUM)
    parser.add_argument("--xenium-h5ad", type=Path, default=DEFAULT_XENIUM)
    parser.add_argument("--rctd-dir", type=Path, default=DEFAULT_RCTD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--mode", choices=("visium", "xenium", "both"), default="both"
    )
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--cot-nitermax", type=int, default=10000)
    parser.add_argument(
        "--max-samples",
        type=int,
        default=0,
        help="Outcome-blind sample limit for smoke tests; 0 uses all samples.",
    )
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: chunk_size and handle.read(chunk_size), b""):
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


def lr_resource() -> pd.DataFrame:
    cellchat = ct.pp.ligand_receptor_database(
        database="CellChat",
        species="human",
        signaling_type="Cell-Cell Contact",
    )
    cellchat = cellchat.iloc[:, :3].copy()
    cellchat.columns = ["ligand", "receptor", "pathway"]
    cellchat["source_database"] = "CellChat"
    cpdb = ct.pp.ligand_receptor_database(
        database="CellPhoneDB_v4.0",
        species="human",
        signaling_type="Cell-Cell Contact",
    )
    cpdb = cpdb.iloc[:, :3].copy()
    cpdb.columns = ["ligand", "receptor", "pathway"]
    cpdb["source_database"] = "CellPhoneDB_v4.0"
    combined = pd.concat([cellchat, cpdb], ignore_index=True)
    target = pd.MultiIndex.from_tuples(PRIMARY_ROUTES)
    combined_index = pd.MultiIndex.from_frame(combined[["ligand", "receptor"]])
    combined = combined[combined_index.isin(target)].copy()
    sources = (
        combined.groupby(["ligand", "receptor"], as_index=False)["source_database"]
        .agg(lambda values: ";".join(sorted(set(values))))
    )
    resource = pd.DataFrame(PRIMARY_ROUTES, columns=["ligand", "receptor"])
    resource["pathway"] = "NOTCH"
    resource = resource.merge(
        sources, on=["ligand", "receptor"], how="left", validate="one_to_one"
    )
    if resource["source_database"].isna().any():
        raise ValueError("A primary route is absent from both COMMOT source databases")
    return resource


def load_rctd_weights(
    adata: ad.AnnData, path: Path
) -> tuple[pd.DataFrame, float, int, int]:
    weights = pd.read_csv(path, index_col=0)
    weights.index = weights.index.astype(str)
    common = adata.obs_names.intersection(weights.index, sort=False)
    coverage = len(common) / adata.n_obs
    aligned_weights = weights.reindex(adata.obs_names).fillna(0.0)
    return aligned_weights, coverage, len(weights), len(common)


def align_rctd(
    adata: ad.AnnData, path: Path, minimum_coverage: float = 0.50
) -> tuple[ad.AnnData, pd.DataFrame, float]:
    aligned_weights, coverage, _, n_common = load_rctd_weights(adata, path)
    if coverage < minimum_coverage:
        raise ValueError(
            f"{path.name} covers only {n_common}/{adata.n_obs} "
            f"spatial units ({coverage:.1%})"
        )
    return adata, aligned_weights, coverage


def geometry_threshold(xy: np.ndarray, k: int) -> float:
    if len(xy) <= k:
        raise ValueError(f"Too few spatial units for k={k}: {len(xy)}")
    distances = cKDTree(xy).query(xy, k=k + 1)[0][:, -1]
    return float(np.median(distances) * 1.05)


def normalise_primary_genes(adata: ad.AnnData) -> ad.AnnData:
    genes = sorted({gene for pair in PRIMARY_ROUTES for gene in pair})
    missing = sorted(set(genes).difference(adata.var_names))
    if missing:
        raise ValueError(f"Spatial sample is missing genes: {missing}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata[:, genes].copy()


def commot_matrices(
    adata: ad.AnnData,
    resource: pd.DataFrame,
    distance_threshold: float,
    cot_nitermax: int,
) -> dict[tuple[str, str], sp.spmatrix]:
    ct.tl.spatial_communication(
        adata,
        database_name=DATABASE_NAME,
        df_ligrec=resource[["ligand", "receptor", "pathway"]],
        pathway_sum=True,
        dis_thr=distance_threshold,
        cot_nitermax=cot_nitermax,
    )
    matrices = {}
    for ligand, receptor in PRIMARY_ROUTES:
        key = f"commot-{DATABASE_NAME}-{ligand}-{receptor}"
        if key not in adata.obsp:
            raise KeyError(f"COMMOT did not create {key}")
        matrices[(ligand, receptor)] = adata.obsp[key]
    return matrices


def weighted_mass(
    matrix: sp.spmatrix, sender_weight: np.ndarray, receiver_weight: np.ndarray
) -> float:
    return float(sender_weight @ (matrix @ receiver_weight))


def process_visium(
    path: Path,
    rctd_dir: Path,
    resource: pd.DataFrame,
    outdir: Path,
    max_samples: int,
    cot_nitermax: int,
) -> pd.DataFrame:
    backed = ad.read_h5ad(path, backed="r")
    required_obs = {"sample", "condition"}
    missing_obs = required_obs.difference(backed.obs.columns)
    if missing_obs or "spatial" not in backed.obsm:
        raise ValueError(
            f"Visium object contract failed: missing obs={sorted(missing_obs)}"
        )
    samples = sorted(backed.obs["sample"].astype(str).unique())
    if len(samples) != 14:
        raise ValueError(f"Expected 14 Visium samples, found {len(samples)}")
    if max_samples:
        samples = samples[:max_samples]

    rows = []
    for sample in samples:
        mask = backed.obs["sample"].astype(str).eq(sample).to_numpy()
        adata = backed[mask].to_memory()
        condition_values = adata.obs["condition"].astype(str).unique()
        if len(condition_values) != 1:
            raise ValueError(f"Non-unique condition for {sample}: {condition_values}")
        condition = condition_values[0]
        weights_path = rctd_dir / f"Zenodo_visium__{sample}_weights.csv"
        n_input_spots = adata.n_obs
        adata, weights, rctd_coverage = align_rctd(adata, weights_path)
        required_weights = {"Endothelial", "Pericyte", "Fib_proFibrotic"}
        absent = required_weights.difference(weights.columns)
        if absent:
            raise ValueError(f"{weights_path.name} missing RCTD columns: {sorted(absent)}")
        xy = np.asarray(adata.obsm["spatial"], dtype=float)
        distance_threshold = geometry_threshold(xy, k=6)
        n_edges = int(
            cKDTree(xy).query_pairs(distance_threshold, output_type="ndarray").shape[0]
        )
        adata = normalise_primary_genes(adata)
        matrices = commot_matrices(
            adata, resource, distance_threshold, cot_nitermax
        )
        for ligand, receptor in PRIMARY_ROUTES:
            sender = (
                weights["Endothelial"] + weights["Pericyte"]
                if ligand == "JAG1"
                else weights["Endothelial"]
            ).to_numpy(dtype=float)
            receiver = weights["Fib_proFibrotic"].to_numpy(dtype=float)
            mass = weighted_mass(matrices[(ligand, receptor)], sender, receiver)
            rows.append(
                {
                    "platform": "Visium",
                    "sample": sample,
                    "condition": condition,
                    "ligand": ligand,
                    "receptor": receptor,
                    "n_input_spots": n_input_spots,
                    "n_spots": adata.n_obs,
                    "rctd_coverage": rctd_coverage,
                    "n_spatial_edges": n_edges,
                    "distance_threshold": distance_threshold,
                    "weighted_transport_mass": mass,
                    "mass_per_spot": mass / adata.n_obs,
                    "mass_per_edge": mass / max(n_edges, 1),
                    "primary_score": float(np.log1p(mass / adata.n_obs)),
                    "edge_score": float(np.log1p(mass / max(n_edges, 1))),
                }
            )
        print(f"[COMMOT] Visium {sample}: {adata.n_obs} spots", flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(outdir / "commot_visium_sample_scores.csv", index=False)
    return result


def exact_condition_test(
    frame: pd.DataFrame, observed: float
) -> tuple[float, int]:
    samples = frame["sample"].to_numpy()
    scores = frame["primary_score"].to_numpy(dtype=float)
    n_hc = int(frame["condition"].eq("HC").sum())
    betas = []
    for hc_positions in itertools.combinations(range(len(frame)), n_hc):
        is_ssc = np.ones(len(frame), dtype=float)
        is_ssc[list(hc_positions)] = 0.0
        matrix = sm.add_constant(is_ssc)
        betas.append(float(sm.OLS(scores, matrix).fit().params[1]))
    betas = np.asarray(betas)
    return float(np.mean(np.abs(betas) >= abs(observed) - 1e-15)), len(betas)


def fit_visium(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ligand, receptor in PRIMARY_ROUTES:
        route = frame[
            frame["ligand"].eq(ligand) & frame["receptor"].eq(receptor)
        ].copy()
        is_ssc = route["condition"].eq("SSc").astype(float)
        fit = sm.OLS(route["primary_score"], sm.add_constant(is_ssc)).fit(
            cov_type="HC3"
        )
        beta = float(fit.params["condition"])
        exact_p, assignments = exact_condition_test(route, beta)
        ci = fit.conf_int().loc["condition"]
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "n_samples": len(route),
                "n_SSc": int(route["condition"].eq("SSc").sum()),
                "n_HC": int(route["condition"].eq("HC").sum()),
                "median_score_SSc": float(
                    route.loc[route["condition"].eq("SSc"), "primary_score"].median()
                ),
                "median_score_HC": float(
                    route.loc[route["condition"].eq("HC"), "primary_score"].median()
                ),
                "difference_SSc_minus_HC": beta,
                "ci_low": float(ci.iloc[0]),
                "ci_high": float(ci.iloc[1]),
                "hc3_p": float(fit.pvalues["condition"]),
                "exact_permutation_p": exact_p,
                "permutation_assignments": assignments,
            }
        )
    result = pd.DataFrame(rows)
    result["BH_q"] = bh_adjust(result["exact_permutation_p"])
    result["Holm_p"] = holm_adjust(result["exact_permutation_p"])
    return result


def process_xenium(
    path: Path,
    rctd_dir: Path,
    resource: pd.DataFrame,
    outdir: Path,
    max_samples: int,
    cot_nitermax: int,
) -> pd.DataFrame:
    backed = ad.read_h5ad(path, backed="r")
    samples = sorted(backed.obs["sample"].astype(str).unique())
    if len(samples) != 10:
        raise ValueError(f"Expected 10 Xenium sections, found {len(samples)}")
    if max_samples:
        samples = samples[:max_samples]
    partial_path = outdir / "commot_xenium_section_scores.partial.csv"
    if partial_path.exists():
        partial = pd.read_csv(partial_path)
        rows = partial.to_dict(orient="records")
        completed_samples = set(partial["sample"].astype(str))
    else:
        rows = []
        completed_samples = set()
    coverage_rows = []
    for sample in samples:
        mask = backed.obs["sample"].astype(str).eq(sample).to_numpy()
        adata = backed[mask].to_memory()
        weights_path = rctd_dir / f"GSE312932_Xenium__{sample}_weights.csv"
        n_input_cells = adata.n_obs
        weights, rctd_coverage, n_weight_rows, n_common = load_rctd_weights(
            adata, weights_path
        )
        eligible = rctd_coverage >= 0.50
        coverage_rows.append(
            {
                "sample": sample,
                "n_input_cells": n_input_cells,
                "n_weight_rows": n_weight_rows,
                "n_common_cells": n_common,
                "rctd_coverage": rctd_coverage,
                "eligible_at_50pct": eligible,
                "reason": (
                    "eligible"
                    if eligible
                    else "RCTD coverage below prespecified 50% threshold"
                ),
            }
        )
        if not eligible:
            print(
                f"[COMMOT] Xenium {sample}: excluded at "
                f"{n_common}/{n_input_cells} RCTD coverage "
                f"({rctd_coverage:.1%})",
                flush=True,
            )
            continue
        if sample in completed_samples:
            print(f"[COMMOT] Xenium {sample}: resumed from checkpoint", flush=True)
            continue
        label = weights.idxmax(axis=1)
        confidence = weights.max(axis=1)
        xy = np.asarray(adata.obsm["spatial"], dtype=float)
        distance_threshold = geometry_threshold(xy, k=6)
        adata = normalise_primary_genes(adata)
        matrices = commot_matrices(
            adata, resource, distance_threshold, cot_nitermax
        )
        for ligand, receptor in PRIMARY_ROUTES:
            sender_types = (
                ("Endothelial", "Pericyte")
                if ligand == "JAG1"
                else ("Endothelial",)
            )
            sender = (
                label.isin(sender_types) & confidence.ge(0.50)
            ).astype(float).to_numpy()
            receiver = (
                label.eq("Fib_proFibrotic") & confidence.ge(0.50)
            ).astype(float).to_numpy()
            mass = weighted_mass(matrices[(ligand, receptor)], sender, receiver)
            rows.append(
                {
                    "platform": "Xenium",
                    "sample": sample,
                    "condition": "SSc",
                    "ligand": ligand,
                    "receptor": receptor,
                    "n_input_cells": n_input_cells,
                    "n_cells": adata.n_obs,
                    "rctd_coverage": rctd_coverage,
                    "n_sender_cells": int(sender.sum()),
                    "n_receiver_cells": int(receiver.sum()),
                    "distance_threshold": distance_threshold,
                    "restricted_transport_mass": mass,
                    "mass_per_cell": mass / adata.n_obs,
                    "localisation_score": float(np.log1p(mass / adata.n_obs)),
                }
            )
        print(f"[COMMOT] Xenium {sample}: {adata.n_obs} cells", flush=True)
        pd.DataFrame(rows).to_csv(partial_path, index=False)
    result = pd.DataFrame(rows)
    result.to_csv(outdir / "commot_xenium_section_scores.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(
        outdir / "commot_xenium_coverage_audit.csv", index=False
    )
    return result


def plot_visium(
    scores: pd.DataFrame, results: pd.DataFrame, outdir: Path
) -> None:
    matplotlib.use("Agg")
    sns.set_theme(style="whitegrid", context="paper")
    plot = scores.copy()
    plot["route"] = plot["ligand"] + "-" + plot["receptor"]
    order = [f"{ligand}-{receptor}" for ligand, receptor in PRIMARY_ROUTES]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    sns.stripplot(
        data=plot,
        x="route",
        y="primary_score",
        hue="condition",
        order=order,
        palette={"HC": "#4477AA", "SSc": "#CC6677"},
        dodge=True,
        size=5,
        ax=axes[0],
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel("log1p RCTD-weighted COMMOT mass per spot")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].legend(title="", frameon=False)

    forest = results.copy()
    forest["route"] = forest["ligand"] + "-" + forest["receptor"]
    forest = forest.set_index("route").loc[order].reset_index()
    y = np.arange(len(forest))
    axes[1].axvline(0, color="#777777", linewidth=0.8)
    axes[1].errorbar(
        forest["difference_SSc_minus_HC"],
        y,
        xerr=[
            forest["difference_SSc_minus_HC"] - forest["ci_low"],
            forest["ci_high"] - forest["difference_SSc_minus_HC"],
        ],
        fmt="o",
        color="#228833",
        capsize=3,
    )
    axes[1].set_yticks(y, forest["route"])
    axes[1].set_xlabel("SSc-HC difference (HC3 95% CI)")
    axes[1].set_ylabel("")
    for index, row in forest.iterrows():
        axes[1].text(
            axes[1].get_xlim()[1],
            index,
            f"q={row.BH_q:.3f}",
            ha="right",
            va="bottom",
            fontsize=7,
        )
    fig.tight_layout()
    fig.savefig(outdir / "COMMOT_Visium_primary.pdf", bbox_inches="tight")
    fig.savefig(
        outdir / "COMMOT_Visium_primary.png", dpi=400, bbox_inches="tight"
    )
    plt.close(fig)


def main() -> int:
    args = parse_args()
    started = time.time()
    outdir = args.outdir
    if args.max_samples:
        outdir = outdir / f"smoke_{args.max_samples}_samples_{args.mode}"
    outdir.mkdir(parents=True, exist_ok=True)
    resource = lr_resource()
    resource.to_csv(outdir / "commot_notch_resource.tsv", sep="\t", index=False)

    visium_scores = None
    visium_results = None
    xenium_scores = None
    if args.mode in ("visium", "both"):
        visium_scores = process_visium(
            args.visium_h5ad,
            args.rctd_dir,
            resource,
            outdir,
            args.max_samples,
            args.cot_nitermax,
        )
        if not args.max_samples:
            visium_results = fit_visium(visium_scores)
            visium_results.to_csv(
                outdir / "commot_visium_condition_results.csv", index=False
            )
            plot_visium(visium_scores, visium_results, outdir)
    if args.mode in ("xenium", "both"):
        xenium_scores = process_xenium(
            args.xenium_h5ad,
            args.rctd_dir,
            resource,
            outdir,
            args.max_samples,
            args.cot_nitermax,
        )

    inputs = {}
    for name, path in (
        ("visium_h5ad", args.visium_h5ad),
        ("xenium_h5ad", args.xenium_h5ad),
    ):
        if path.exists() and (
            args.mode == "both" or name.startswith(args.mode)
        ):
            inputs[name] = {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": vars(args) | {"outdir": str(outdir)},
        "inputs": inputs,
        "resource_rows": resource.to_dict(orient="records"),
        "visium_samples": (
            0 if visium_scores is None else int(visium_scores["sample"].nunique())
        ),
        "xenium_sections": (
            0 if xenium_scores is None else int(xenium_scores["sample"].nunique())
        ),
        "runtime_seconds": time.time() - started,
        "versions": {
            "python": platform.python_version(),
            "anndata": package_version("anndata"),
            "commot": package_version("commot"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scanpy": package_version("scanpy"),
            "scipy": scipy.__version__,
            "statsmodels": package_version("statsmodels"),
        },
    }
    coverage_path = outdir / "commot_xenium_coverage_audit.csv"
    if coverage_path.exists():
        coverage = pd.read_csv(coverage_path)
        manifest["xenium_screened_sections"] = int(len(coverage))
        manifest["xenium_excluded_sections"] = coverage.loc[
            ~coverage["eligible_at_50pct"].astype(bool), "sample"
        ].astype(str).tolist()
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str), encoding="utf-8"
    )
    if not args.max_samples:
        (outdir / "COMMOT_SPATIAL_DONE").write_text(
            datetime.now(timezone.utc).isoformat() + "\n", encoding="ascii"
        )
    if visium_results is not None:
        print(visium_results.to_string(index=False))
    print(f"COMMOT spatial analysis complete: {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
