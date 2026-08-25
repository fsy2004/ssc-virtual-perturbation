#!/usr/bin/env python
"""Donor-level cell composition analysis for revised Figure 1c.

Primary endpoints are the total fibroblast fraction and each of the seven
atlas-defined fibroblast-state contributions, all divided by the number of
profiled skin cells per resolved donor. Beta-binomial models adjust for cohort.
The complete 14-cell-type family and within-fibroblast denominators are exported
as pre-specified supporting analyses.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "server_archive" / "fig_atlas" / "atlas_fig.h5ad"
DONOR_MAP = ROOT / "04_manuscript" / "revision_20260722" / "donor_metadata" / "outputs" / "sample_to_donor.csv"
OUT = ROOT / "04_manuscript" / "revision_20260722" / "figure1" / "outputs"
PLOT_DATA = ROOT / "04_manuscript" / "plot_data_local" / "figure1_revision"
R_SCRIPT = ROOT / "06_code_reproducibility" / "01_scrna_atlas" / "04_donor_composition_glmmTMB.R"
SCCOMP_R_SCRIPT = ROOT / "06_code_reproducibility" / "01_scrna_atlas" / "04_donor_composition_sccomp.R"
RSCRIPT = os.environ.get("RSCRIPT", r"C:\Program Files\R\R-4.4.3\bin\Rscript.exe")

SUBTYPES = [
    "Myofibroblast",
    "SFRP4_proFib",
    "SFRP2_DPP4",
    "Adipogenic",
    "FMO1_LSP1",
    "LGR5_Gur",
    "Inflammatory",
]
DISPLAY = {
    "Fibroblast": "Fibroblast total",
    "Myofibroblast": "Myofibroblast",
    "SFRP4_proFib": "SFRP4+ pro-fibrotic",
    "SFRP2_DPP4": "SFRP2+/DPP4+",
    "Adipogenic": "Adipogenic",
    "FMO1_LSP1": "FMO1+/LSP1+",
    "LGR5_Gur": "LGR5+",
    "Inflammatory": "Inflammatory-like*",
}
COHORT_COLORS = {
    "Gur_GSE195452": "#0072B2",
    "GSE249279": "#E69F00",
    "Tabib_GSE138669": "#009E73",
    "GSE292979": "#CC79A7",
    "GSE236111": "#56B4E9",
}
CONDITION_COLORS = {"HC": "#2166AC", "SSc": "#D55E00"}


def log(message: str) -> None:
    print(f"[figure1c] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_kept_cells() -> tuple[pd.DataFrame, dict]:
    donor_map = pd.read_csv(DONOR_MAP, dtype={"sample": str, "donor_id": str})
    donor_map["analysis_keep"] = donor_map["analysis_keep"].astype(str).str.lower().eq("true")
    kept_map = donor_map[donor_map["analysis_keep"]][
        ["sample", "donor_id", "cohort", "condition"]
    ].copy()
    if kept_map["sample"].duplicated().any():
        raise ValueError("Retained sample-to-donor mapping is not unique")

    atlas = ad.read_h5ad(ATLAS, backed="r")
    obs = atlas.obs[["sample", "celltype", "fib_subtype"]].copy().reset_index(names="cell_id")
    atlas_shape = [int(atlas.n_obs), int(atlas.n_vars)]
    atlas.file.close()
    obs["sample"] = obs["sample"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    obs["fib_subtype"] = obs["fib_subtype"].astype(str)
    obs = obs.merge(kept_map, on="sample", how="inner", validate="many_to_one")
    if obs.empty:
        raise ValueError("No atlas cells remain after donor/tissue/deduplication filters")
    donor_conditions = obs.groupby("donor_id")["condition"].nunique()
    if (donor_conditions > 1).any():
        raise ValueError("Resolved donor has conflicting conditions")
    audit = {
        "atlas_shape": atlas_shape,
        "retained_cells": int(len(obs)),
        "retained_donors": int(obs["donor_id"].nunique()),
        "retained_libraries": int(obs["sample"].nunique()),
        "retained_by_condition": obs.drop_duplicates("donor_id")["condition"].value_counts().to_dict(),
        "retained_by_cohort": obs.drop_duplicates("donor_id")["cohort"].value_counts().to_dict(),
    }
    return obs, audit


def build_count_table(obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    donor_meta = (
        obs.groupby("donor_id", observed=True)
        .agg(
            cohort=("cohort", "first"),
            condition=("condition", "first"),
            n_libraries=("sample", "nunique"),
            total_cells=("cell_id", "size"),
        )
        .reset_index()
    )
    donor_index = donor_meta.set_index("donor_id")

    cell_counts = (
        obs.groupby(["donor_id", "celltype"], observed=True).size().unstack(fill_value=0)
    )
    subtype_obs = obs[obs["celltype"].eq("Fibroblast") & obs["fib_subtype"].isin(SUBTYPES)]
    subtype_counts = (
        subtype_obs.groupby(["donor_id", "fib_subtype"], observed=True).size().unstack(fill_value=0)
    )
    assigned_fibro = subtype_counts.reindex(columns=SUBTYPES, fill_value=0).sum(axis=1)

    rows: list[dict] = []

    def add_family(family: str, denominator: str, endpoints: list[str], totals: pd.Series, source: pd.DataFrame) -> None:
        for donor_id, meta in donor_index.iterrows():
            total = int(totals.get(donor_id, 0))
            if total <= 0:
                continue
            for endpoint in endpoints:
                count = int(source.get(endpoint, pd.Series(dtype=int)).get(donor_id, 0))
                rows.append(
                    {
                        "family": family,
                        "denominator": denominator,
                        "endpoint": endpoint,
                        "display": DISPLAY.get(endpoint, endpoint),
                        "donor_id": donor_id,
                        "cohort": meta["cohort"],
                        "condition": meta["condition"],
                        "n_libraries": int(meta["n_libraries"]),
                        "count": count,
                        "total": total,
                        "fraction": count / total,
                    }
                )

    total_cells = donor_index["total_cells"]
    primary_source = subtype_counts.reindex(index=donor_index.index, columns=SUBTYPES, fill_value=0).copy()
    primary_source["Fibroblast"] = cell_counts.reindex(donor_index.index, fill_value=0).get(
        "Fibroblast", pd.Series(0, index=donor_index.index)
    )
    add_family(
        "fibro_primary", "all_skin_cells", ["Fibroblast"] + SUBTYPES,
        total_cells, primary_source,
    )
    major_endpoints = sorted(obs["celltype"].unique())
    add_family(
        "major_celltype", "all_skin_cells", major_endpoints,
        total_cells, cell_counts.reindex(donor_index.index, fill_value=0),
    )
    add_family(
        "fibro_within", "assigned_fibroblasts", SUBTYPES,
        assigned_fibro, subtype_counts.reindex(donor_index.index, columns=SUBTYPES, fill_value=0),
    )
    counts = pd.DataFrame(rows).sort_values(
        ["family", "endpoint", "condition", "cohort", "donor_id"]
    )
    return counts, donor_meta


def run_models(counts_path: Path) -> None:
    log("running glmmTMB beta-binomial models")
    commands = [
        ("glmmTMB", [RSCRIPT, str(R_SCRIPT), str(counts_path), str(OUT)]),
        ("sccomp", [RSCRIPT, str(SCCOMP_R_SCRIPT), str(counts_path), str(OUT)]),
    ]
    for name, command in commands:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="")
        if completed.returncode != 0:
            raise RuntimeError(f"{name} R model script failed with exit code {completed.returncode}")


def stat_text(row: pd.Series) -> str:
    if not np.isfinite(row.get("OR", np.nan)):
        return "not estimable"
    q_value = float(row.q_BH)
    q_label = "<0.001" if q_value < 0.001 else f"={q_value:.3f}"
    return (
        f"OR {row.OR:.2f} [{row.CI_low:.2f}, {row.CI_high:.2f}]\n"
        f"BH q{q_label}"
    )


def draw_box(ax: plt.Axes, values: np.ndarray, position: float, color: str) -> None:
    if len(values) == 0:
        return
    bp = ax.boxplot(
        values, positions=[position], widths=0.36, patch_artist=True, showfliers=False,
        medianprops={"color": color, "linewidth": 1.2},
        boxprops={"facecolor": "white", "edgecolor": color, "linewidth": 1.0, "alpha": 0.9},
        whiskerprops={"color": color, "linewidth": 0.9},
        capprops={"color": color, "linewidth": 0.9},
    )
    for patch in bp["boxes"]:
        patch.set_zorder(2)


def plot_primary(counts: pd.DataFrame, results: pd.DataFrame) -> None:
    primary = counts[(counts["family"] == "fibro_primary") & (counts["denominator"] == "all_skin_cells")]
    primary_results = results[
        (results["family"] == "fibro_primary")
        & (results["denominator"] == "all_skin_cells")
        & (results["omitted_cohort"] == "none")
    ].set_index("endpoint")
    order = ["Fibroblast"] + SUBTYPES
    fig, axes = plt.subplots(2, 4, figsize=(11.8, 6.4), constrained_layout=True)
    rng = np.random.default_rng(20260723)
    for panel_i, (ax, endpoint) in enumerate(zip(axes.flat, order)):
        d = primary[primary["endpoint"] == endpoint].copy()
        for x, condition in enumerate(["HC", "SSc"]):
            group = d[d["condition"] == condition]
            values = group["fraction"].to_numpy(float) * 100
            draw_box(ax, values, x, CONDITION_COLORS[condition])
            jitter = rng.normal(0, 0.075, len(group))
            colors = [COHORT_COLORS.get(c, "#777777") for c in group["cohort"]]
            ax.scatter(
                x + jitter, values, s=11, c=colors, alpha=0.58,
                edgecolors="white", linewidths=0.25, rasterized=True, zorder=3,
            )
        max_value = float(d["fraction"].max() * 100) if len(d) else 0.0
        upper = max(1.0, max_value * 1.28)
        if endpoint == "Fibroblast":
            upper = min(102.0, upper)
        ax.set_ylim(0, upper)
        ax.set_xlim(-0.55, 1.55)
        ax.set_xticks([0, 1], ["HC", "SSc"])
        ax.set_title(DISPLAY[endpoint], fontsize=9, fontweight="bold", pad=4)
        if endpoint in primary_results.index:
            ax.text(
                0.98, 0.97, stat_text(primary_results.loc[endpoint]),
                transform=ax.transAxes, ha="right", va="top", fontsize=6.7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
            )
        if endpoint == "Inflammatory":
            n_cells = int(d["count"].sum())
            n_positive = int((d["count"] > 0).sum())
            ax.text(
                0.02, 0.97, f"{n_cells} cells; {n_positive} donor",
                transform=ax.transAxes, ha="left", va="top", fontsize=6.7,
            )
        if panel_i % 4 == 0:
            ax.set_ylabel("Fraction of all skin cells (%)", fontsize=8)
        else:
            ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=7)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.6, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color,
               markeredgecolor="white", markersize=5.5, label=cohort)
        for cohort, color in COHORT_COLORS.items()
    ]
    fig.legend(
        handles=handles, loc="outside lower center", ncol=len(handles),
        frameon=False, fontsize=7, title="Cohort", title_fontsize=7.5,
    )
    fig.suptitle(
        "Donor-level fibroblast composition in SSc and healthy skin",
        fontsize=11, fontweight="bold",
    )
    png = OUT / "Figure1c_donor_fibroblast_composition.png"
    pdf = OUT / "Figure1c_donor_fibroblast_composition.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    Image.open(png).convert("L").save(OUT / "Figure1c_donor_fibroblast_composition_grayscale.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT_DATA.mkdir(parents=True, exist_ok=True)
    obs, audit = load_kept_cells()
    counts, donor_meta = build_count_table(obs)
    counts_path = OUT / "Figure1c_donor_counts.csv"
    counts.to_csv(counts_path, index=False)
    counts.to_csv(PLOT_DATA / "Figure1c_donor_counts.csv", index=False)
    donor_meta.to_csv(OUT / "Figure1c_donor_audit.csv", index=False)
    run_models(counts_path)

    results = pd.read_csv(OUT / "Figure1c_glmmTMB_results.csv")
    primary_results = results[
        (results["family"] == "fibro_primary")
        & (results["denominator"] == "all_skin_cells")
        & (results["omitted_cohort"] == "none")
    ]
    estimable = primary_results[primary_results["OR"].notna()]
    if not estimable["converged"].astype(bool).all():
        raise RuntimeError("At least one estimable primary composition model did not converge")
    plot_primary(counts, results)

    source = counts[
        (counts["family"] == "fibro_primary") & (counts["denominator"] == "all_skin_cells")
    ].merge(
        primary_results[
            ["endpoint", "OR", "CI_low", "CI_high", "p_value", "q_BH", "p_Holm", "model_message"]
        ],
        on="endpoint", how="left", validate="many_to_one",
    )
    source.to_csv(OUT / "Figure1c_source_data.csv", index=False)
    manifest = {
        **audit,
        "independent_unit": "resolved donor",
        "primary_denominator": "all retained SSc/HC skin cells per donor",
        "primary_endpoints": ["Fibroblast"] + SUBTYPES,
        "model": "glmmTMB beta-binomial: cbind(count, total-count) ~ condition + cohort",
        "composition_sensitivity": (
            "sccomp sum-constrained beta-binomial: 14 mutually exclusive major cell types; "
            "seven fibroblast states plus other skin cells; composition ~ condition + cohort"
        ),
        "multiple_testing": "BH and Holm within each endpoint family and denominator",
        "sensitivity_analyses": [
            "all 14 major cell types with all-skin-cell denominator",
            "seven fibroblast states with assigned-fibroblast denominator",
            "leave-one-cohort-out refits",
        ],
        "atlas_sha256": sha256(ATLAS),
        "donor_map_sha256": sha256(DONOR_MAP),
        "r_script_sha256": sha256(R_SCRIPT),
        "sccomp_r_script_sha256": sha256(SCCOMP_R_SCRIPT),
    }
    (OUT / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    log(f"completed donor composition analysis: {OUT}")


if __name__ == "__main__":
    main()
