#!/usr/bin/env python
"""True-donor Figure 6d analysis of myofibroblast receiver-gene detection.

Retained skin libraries are mapped to deposited donor identifiers before any
counting. Multiple sorted or technical libraries from one participant are
summed within donor. The primary beta-binomial model includes every donor with
at least one retained myofibroblast; denominator sensitivities use outcome-blind
quartiles of the donor myofibroblast-count distribution.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "server_archive" / "fig_atlas" / "atlas_fig.h5ad"
DONOR_MAP = (
    ROOT
    / "04_manuscript"
    / "revision_20260722"
    / "donor_metadata"
    / "outputs"
    / "sample_to_donor.csv"
)
OUT = ROOT / "04_manuscript" / "revision_20260722" / "figure6d" / "outputs"
PLOT_DATA = ROOT / "04_manuscript" / "plot_data_local" / "ccc_nichenet" / "donor_notch_receiver"
R_SCRIPT = ROOT / "06_code_reproducibility" / "06_ccc_nichenet" / "04_donor_notch_receiver_glmmTMB.R"
RSCRIPT = os.environ.get("RSCRIPT", r"C:\Program Files\R\R-4.4.3\bin\Rscript.exe")

GENES = ["NOTCH1", "NOTCH2", "NOTCH3", "HES1"]
MYO_LABEL = "Myofibroblast"
CONDITION_COLORS = {"HC": "#2166AC", "SSc": "#D55E00"}
COHORT_COLORS = {
    "Gur_GSE195452": "#0072B2",
    "GSE249279": "#E69F00",
    "Tabib_GSE138669": "#009E73",
    "GSE292979": "#CC79A7",
    "GSE236111": "#56B4E9",
}


def log(message: str) -> None:
    print(f"[fig6d] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_retained_map() -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = pd.read_csv(DONOR_MAP, dtype={"sample": str, "donor_id": str})
    mapping["analysis_keep"] = mapping["analysis_keep"].astype(str).str.lower().eq("true")
    kept = mapping.loc[
        mapping["analysis_keep"], ["sample", "donor_id", "cohort", "condition"]
    ].copy()
    if kept["sample"].duplicated().any():
        raise ValueError("Retained sample-to-donor map is not unique")
    consistency = kept.groupby("donor_id")[["cohort", "condition"]].nunique()
    if (consistency > 1).any().any():
        raise ValueError("A retained donor spans multiple cohorts or conditions")
    roster = (
        kept.groupby("donor_id", observed=True)
        .agg(
            cohort=("cohort", "first"),
            condition=("condition", "first"),
            n_libraries=("sample", "nunique"),
            sample_ids=("sample", lambda x: ";".join(sorted(x))),
        )
        .reset_index()
    )
    if len(roster) != 230:
        raise ValueError(f"Expected 230 retained skin donors, found {len(roster)}")
    return kept, roster


def extract_counts() -> tuple[pd.DataFrame, dict, list[int]]:
    kept, roster = load_retained_map()
    atlas = ad.read_h5ad(ATLAS, backed="r")
    missing = [gene for gene in GENES if gene not in atlas.var_names]
    if missing:
        atlas.file.close()
        raise ValueError(f"Target genes absent from atlas: {missing}")
    required = {"sample", "fib_subtype"}
    if not required.issubset(atlas.obs.columns):
        atlas.file.close()
        raise ValueError(f"Atlas obs is missing: {sorted(required - set(atlas.obs.columns))}")

    sample = atlas.obs["sample"].astype(str)
    fib_subtype = atlas.obs["fib_subtype"].astype(str)
    retained_mask = sample.isin(kept["sample"]).to_numpy()
    myo_mask = retained_mask & fib_subtype.eq(MYO_LABEL).to_numpy()
    myo_meta = pd.DataFrame({"sample": sample.to_numpy()[myo_mask]})
    myo_meta = myo_meta.merge(kept, on="sample", how="left", validate="many_to_one")
    if myo_meta["donor_id"].isna().any():
        atlas.file.close()
        raise ValueError("At least one retained myofibroblast lacks donor metadata")

    gene_indices = [int(atlas.var_names.get_loc(gene)) for gene in GENES]
    matrix = atlas.X[myo_mask, :][:, gene_indices]
    positive = matrix.toarray() > 0 if sp.issparse(matrix) else np.asarray(matrix) > 0
    atlas_shape = [int(atlas.n_obs), int(atlas.n_vars)]
    atlas.file.close()
    log(f"retained myofibroblasts: {len(myo_meta):,}")

    rows: list[pd.DataFrame] = []
    for gene_index, gene in enumerate(GENES):
        gene_meta = myo_meta[["donor_id", "cohort", "condition"]].copy()
        gene_meta["positive_cell"] = positive[:, gene_index].astype(int)
        aggregated = (
            gene_meta.groupby(["donor_id", "cohort", "condition"], observed=True)
            .agg(positive=("positive_cell", "sum"), total=("positive_cell", "size"))
            .reset_index()
        )
        aggregated["gene"] = gene
        rows.append(aggregated)
    observed = pd.concat(rows, ignore_index=True)

    complete = pd.MultiIndex.from_product(
        [roster["donor_id"], GENES], names=["donor_id", "gene"]
    ).to_frame(index=False)
    counts = complete.merge(
        observed[["donor_id", "gene", "positive", "total"]],
        on=["donor_id", "gene"],
        how="left",
        validate="one_to_one",
    ).merge(roster, on="donor_id", how="left", validate="many_to_one")
    counts[["positive", "total"]] = counts[["positive", "total"]].fillna(0).astype(int)
    counts["negative"] = counts["total"] - counts["positive"]
    counts["fraction"] = np.where(counts["total"] > 0, counts["positive"] / counts["total"], np.nan)
    counts = counts[
        [
            "gene", "donor_id", "cohort", "condition", "positive", "negative",
            "total", "fraction", "n_libraries", "sample_ids",
        ]
    ].sort_values(["gene", "condition", "cohort", "donor_id"])

    denominator = counts.loc[counts["gene"].eq(GENES[0]) & counts["total"].gt(0), "total"]
    if denominator.empty:
        raise ValueError("No retained donor has a myofibroblast denominator")
    quartiles = np.quantile(denominator.to_numpy(), [0.25, 0.50, 0.75], method="nearest")
    thresholds = sorted({1, *(int(value) for value in quartiles)})

    audit = {
        "atlas_shape": atlas_shape,
        "retained_skin_cells": int(kept["sample"].map(
            pd.read_csv(DONOR_MAP, dtype={"sample": str}).set_index("sample")["n_cells_atlas"]
        ).astype(int).sum()),
        "retained_donors": int(roster["donor_id"].nunique()),
        "retained_libraries": int(kept["sample"].nunique()),
        "retained_donors_by_condition": roster["condition"].value_counts().to_dict(),
        "retained_donors_by_cohort": roster["cohort"].value_counts().to_dict(),
        "myofibroblast_cells": int(len(myo_meta)),
        "donors_with_myofibroblasts": int(denominator.shape[0]),
        "donors_without_myofibroblasts": int(len(roster) - denominator.shape[0]),
        "outcome_blind_threshold_rule": "1 plus nearest 25th, 50th and 75th percentiles of donor myofibroblast totals",
        "thresholds_min_total": thresholds,
        "myofibroblast_denominator_quantiles": {
            "q25": int(quartiles[0]), "median": int(quartiles[1]), "q75": int(quartiles[2])
        },
    }
    return counts, audit, thresholds


def build_audit_tables(counts: pd.DataFrame, thresholds: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    donor_audit = (
        counts.groupby(["gene", "condition"], observed=True)
        .agg(
            donors=("donor_id", "nunique"),
            donors_with_myo=("total", lambda x: int((x > 0).sum())),
            median_total=("total", "median"),
            min_total=("total", "min"),
            max_total=("total", "max"),
        )
        .reset_index()
    )
    rows: list[dict] = []
    for threshold in thresholds:
        eligible = counts[counts["total"].ge(threshold)]
        for gene, group in eligible.groupby("gene", observed=True):
            condition_counts = group["condition"].value_counts()
            rows.append(
                {
                    "threshold_min_total": threshold,
                    "gene": gene,
                    "n_donors": int(group["donor_id"].nunique()),
                    "n_SSc": int(condition_counts.get("SSc", 0)),
                    "n_HC": int(condition_counts.get("HC", 0)),
                }
            )
    return donor_audit, pd.DataFrame(rows)


def run_models(counts_path: Path, thresholds: list[int]) -> pd.DataFrame:
    command = [RSCRIPT, str(R_SCRIPT), str(counts_path), str(OUT), ",".join(map(str, thresholds))]
    log(f"running donor beta-binomial models at thresholds {thresholds}")
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
        raise RuntimeError(f"Figure 6d R model failed with exit code {completed.returncode}")
    return pd.read_csv(OUT / "Figure6d_glmmTMB_results.csv")


def draw_box(ax: plt.Axes, values: np.ndarray, x: float, color: str) -> None:
    if not len(values):
        return
    ax.boxplot(
        values,
        positions=[x],
        widths=0.38,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": color, "linewidth": 1.2},
        boxprops={"facecolor": "white", "edgecolor": color, "linewidth": 1.0},
        whiskerprops={"color": color, "linewidth": 0.9},
        capprops={"color": color, "linewidth": 0.9},
    )


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"={value:.3f}"


def draw_panel(counts: pd.DataFrame, results: pd.DataFrame) -> None:
    plot_counts = counts[counts["total"].ge(1)].copy()
    main = results[results["threshold_min_total"].eq(1)].set_index("gene")
    rng = np.random.default_rng(20260723)
    fig, axes = plt.subplots(1, len(GENES), figsize=(9.2, 3.15), sharey=True)
    root_denominator = np.sqrt(plot_counts["total"].to_numpy(float))
    denominator_min = float(root_denominator.min())
    denominator_max = float(root_denominator.max())

    for gene_index, (ax, gene) in enumerate(zip(axes, GENES)):
        gene_data = plot_counts[plot_counts["gene"].eq(gene)]
        for x, condition in enumerate(["HC", "SSc"]):
            group = gene_data[gene_data["condition"].eq(condition)]
            values = group["fraction"].to_numpy(float)
            draw_box(ax, values, x, CONDITION_COLORS[condition])
            root_size = np.sqrt(group["total"].to_numpy(float))
            if denominator_max > denominator_min:
                sizes = 12 + (root_size - denominator_min) / (denominator_max - denominator_min) * 58
            else:
                sizes = np.repeat(30.0, len(group))
            ax.scatter(
                x + rng.normal(0, 0.065, len(group)),
                values,
                s=sizes,
                c=[COHORT_COLORS.get(cohort, "#777777") for cohort in group["cohort"]],
                alpha=0.62,
                edgecolors="white",
                linewidths=0.3,
                rasterized=True,
                zorder=3,
            )
        row = main.loc[gene]
        if np.isfinite(row.OR):
            title = (
                f"{gene}\nOR {row.OR:.2f} [{row.CI_low:.2f}, {row.CI_high:.2f}]\n"
                f"BH q{format_p(float(row.q_BH))}; Holm p{format_p(float(row.p_Holm))}"
            )
        else:
            title = f"{gene}\n{row.model_message}"
        ax.set_title(title, fontsize=7.5, fontweight="bold", linespacing=1.15)
        ax.set_xticks([0, 1], ["HC", "SSc"])
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(axis="y", color="#E1E1E1", linewidth=0.55, linestyle=":")
        ax.spines[["top", "right"]].set_visible(False)
        if gene_index:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", left=False)
    axes[0].set_ylabel("Positive myofibroblast fraction", fontsize=8)

    cohort_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor=color,
            markeredgecolor="white", markersize=5.2, label=cohort,
        )
        for cohort, color in COHORT_COLORS.items()
    ]
    fig.legend(
        handles=cohort_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=len(cohort_handles),
        frameon=False,
        fontsize=6.8,
        title="Cohort (point area scales with donor myofibroblast count)",
        title_fontsize=7,
    )
    fig.suptitle("True-donor receiver-gene detection in myofibroblasts", fontsize=10, fontweight="bold", y=1.03)
    fig.tight_layout()
    fig.savefig(OUT / "Figure6d_donor_notch_receiver.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "Figure6d_donor_notch_receiver.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    source = plot_counts.merge(
        main.reset_index()[
            ["gene", "n_donors", "n_SSc", "n_HC", "OR", "CI_low", "CI_high", "p_value", "q_BH", "p_Holm", "converged", "model_message"]
        ],
        on="gene",
        how="left",
        validate="many_to_one",
    )
    source.to_csv(OUT / "Figure6d_source_data.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT_DATA.mkdir(parents=True, exist_ok=True)
    counts, audit, thresholds = extract_counts()
    counts_path = OUT / "Figure6d_donor_counts.csv"
    counts.to_csv(counts_path, index=False)
    donor_audit, threshold_audit = build_audit_tables(counts, thresholds)
    donor_audit.to_csv(OUT / "Figure6d_donor_audit.csv", index=False)
    threshold_audit.to_csv(OUT / "Figure6d_threshold_audit.csv", index=False)
    results = run_models(counts_path, thresholds)
    primary = results[results["threshold_min_total"].eq(1) & results["OR"].notna()]
    if len(primary) != len(GENES) or not primary["converged"].astype(bool).all():
        raise RuntimeError("At least one primary Figure 6d model failed")
    draw_panel(counts, results)

    for path in OUT.glob("Figure6d_*"):
        if path.is_file():
            shutil.copy2(path, PLOT_DATA / path.name)
    manifest = {
        **audit,
        "input_atlas": str(ATLAS),
        "donor_map": str(DONOR_MAP),
        "independent_unit": "deposited donor/patient identifier",
        "technical_library_handling": "sum positive and total myofibroblast cells within donor",
        "positive_definition": "atlas expression value > 0 in retained Myofibroblast cells",
        "primary_model": "glmmTMB beta-binomial: cbind(positive, total-positive) ~ condition + cohort",
        "multiple_testing": "BH and Holm across NOTCH1, NOTCH2, NOTCH3 and HES1 within threshold",
        "atlas_sha256": sha256(ATLAS),
        "donor_map_sha256": sha256(DONOR_MAP),
        "r_script_sha256": sha256(R_SCRIPT),
    }
    (OUT / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    log(f"completed: {OUT}")


if __name__ == "__main__":
    main()
