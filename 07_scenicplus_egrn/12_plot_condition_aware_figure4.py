#!/usr/bin/env python3
"""Create a publication candidate for the condition-aware Figure 4 extension.

The script writes a separate candidate figure. It never replaces the current
manuscript Figure 4 automatically.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PRIORITY_TFS = ("HES1", "RBPJ", "SMAD3")
COLORS = {
    "positive": "#B9473F",
    "negative": "#2C6E9B",
    "neutral": "#B8B8B8",
    "direct": "#006D77",
    "extended": "#E08E45",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dar", type=Path, required=True)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--label-top", type=int, default=8)
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return path


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.04,
        label,
        transform=axis.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )


def main() -> None:
    args = parse_args()
    dar = pd.read_csv(require_file(args.dar))
    support = pd.read_csv(require_file(args.support))
    dar_required = {"logFC", "logCPM", "FDR"}
    support_required = {
        "scope",
        "TF",
        "n_DAR_tested_regions",
        "SSc_positive_region_fraction",
        "median_logFC",
        "matched_null_median",
        "matched_null_q025",
        "matched_null_q975",
        "median_logFC_minus_null",
        "empirical_competition_bh_q",
    }
    if missing := dar_required.difference(dar.columns):
        raise ValueError(f"DAR input lacks {sorted(missing)}")
    if missing := support_required.difference(support.columns):
        raise ValueError(f"Support input lacks {sorted(missing)}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(14.2, 10.3), facecolor="white")
    grid = GridSpec(
        2,
        2,
        figure=figure,
        width_ratios=[1.05, 0.95],
        height_ratios=[1.0, 1.05],
        wspace=0.34,
        hspace=0.36,
    )

    axis = figure.add_subplot(grid[0, 0])
    significant = dar["FDR"] < 0.05
    point_colors = np.where(
        significant & (dar["logFC"] > 0),
        COLORS["positive"],
        np.where(
            significant & (dar["logFC"] < 0),
            COLORS["negative"],
            COLORS["neutral"],
        ),
    )
    axis.scatter(
        dar["logCPM"],
        dar["logFC"],
        c=point_colors,
        s=np.where(significant, 9, 4),
        alpha=np.where(significant, 0.72, 0.24),
        linewidths=0,
        rasterized=True,
    )
    axis.axhline(0, color="#555555", linewidth=0.8)
    axis.set_xlabel("Mean accessibility (edgeR logCPM)")
    axis.set_ylabel("SSc vs HC accessibility (logFC)")
    axis.set_title("Donor-pseudobulk fibroblast accessibility")
    axis.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                color=COLORS["positive"],
                label="SSc-positive, FDR < 0.05",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                color=COLORS["negative"],
                label="SSc-negative, FDR < 0.05",
            ),
        ],
        frameon=False,
        fontsize=7,
        loc="upper right",
    )
    axis.spines[["top", "right"]].set_visible(False)
    panel_label(axis, "A")

    axis = figure.add_subplot(grid[0, 1])
    direct = support.loc[
        support["scope"].eq("direct")
        & support["empirical_competition_bh_q"].notna()
    ].copy()
    if direct.empty:
        raise RuntimeError("No direct TF competition result is available to plot")
    direct["minus_log10_q"] = -np.log10(
        direct["empirical_competition_bh_q"].clip(lower=1e-300)
    )
    colors = np.where(
        direct["TF"].isin(PRIORITY_TFS),
        COLORS["positive"],
        np.where(
            direct["median_logFC_minus_null"] > 0,
            COLORS["direct"],
            COLORS["neutral"],
        ),
    )
    sizes = 15 + 120 * np.sqrt(
        direct["n_DAR_tested_regions"]
        / direct["n_DAR_tested_regions"].max()
    )
    axis.scatter(
        direct["median_logFC_minus_null"],
        direct["minus_log10_q"],
        c=colors,
        s=sizes,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.5,
    )
    axis.axvline(0, color="#666666", linewidth=0.8)
    axis.axhline(-np.log10(0.05), color="#777777", linewidth=0.7, linestyle="--")
    labels = pd.concat(
        [
            direct.loc[direct["TF"].isin(PRIORITY_TFS)],
            direct.nlargest(args.label_top, "minus_log10_q"),
        ]
    ).drop_duplicates("TF")
    for row in labels.itertuples():
        axis.annotate(
            row.TF,
            (row.median_logFC_minus_null, row.minus_log10_q),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7,
        )
    axis.set_xlabel("Median logFC minus matched-region null")
    axis.set_ylabel(r"$-\log_{10}$(empirical BH q)")
    axis.set_title("All-TF direct-region competition")
    axis.text(
        0.02,
        0.98,
        "Point area = tested eRegulon regions",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    axis.spines[["top", "right"]].set_visible(False)
    panel_label(axis, "B")

    axis = figure.add_subplot(grid[1, 0])
    priority = support.loc[support["TF"].isin(PRIORITY_TFS)].copy()
    observed_tfs = set(priority["TF"])
    if observed_tfs != set(PRIORITY_TFS):
        raise RuntimeError(
            "Priority TF support rows are incomplete: "
            f"missing={sorted(set(PRIORITY_TFS) - observed_tfs)}, "
            f"extra={sorted(observed_tfs - set(PRIORITY_TFS))}"
        )
    priority["display"] = (
        priority["TF"]
        + " "
        + priority["scope"].map({"direct": "direct", "extended": "extended"})
    )
    priority = priority.sort_values(["TF", "scope"], ascending=[False, True])
    y = np.arange(len(priority))
    competition_ready = priority["matched_null_median"].notna()
    axis.hlines(
        y[competition_ready],
        priority.loc[competition_ready, "matched_null_q025"],
        priority.loc[competition_ready, "matched_null_q975"],
        color="#9A9A9A",
        linewidth=2.0,
        label="Matched-null 95% interval",
    )
    axis.scatter(
        priority.loc[competition_ready, "matched_null_median"],
        y[competition_ready],
        marker="|",
        s=120,
        color="#444444",
        linewidth=1.5,
        zorder=3,
    )
    axis.scatter(
        priority["median_logFC"],
        y,
        s=56,
        color=priority["scope"].map(COLORS),
        edgecolor="white",
        linewidth=0.7,
        zorder=4,
    )
    for row_y, row in zip(y[~competition_ready], priority.loc[~competition_ready].itertuples()):
        axis.annotate(
            f"n={row.n_DAR_tested_regions}; below 10-region gate",
            (row.median_logFC, row_y),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=6.5,
            color="#555555",
        )
    axis.axvline(0, color="#666666", linewidth=0.7)
    axis.set_yticks(y)
    axis.set_yticklabels(priority["display"], fontsize=8)
    axis.set_xlabel("Median SSc-vs-HC peak logFC")
    axis.set_title("Prespecified pathway TF region support")
    axis.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                color=COLORS["direct"],
                label="Observed direct",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                color=COLORS["extended"],
                label="Observed extended",
            ),
            Line2D(
                [0],
                [0],
                color="#9A9A9A",
                linewidth=2,
                label="Matched-null 95% interval (when n >= 10)",
            ),
        ],
        frameon=False,
        fontsize=7,
        loc="best",
    )
    axis.spines[["top", "right"]].set_visible(False)
    panel_label(axis, "C")

    axis = figure.add_subplot(grid[1, 1])
    priority = priority.sort_values(["TF", "scope"]).reset_index(drop=True)
    x = np.arange(len(priority))
    scatter = axis.scatter(
        x,
        priority["SSc_positive_region_fraction"],
        s=35 + 120 * np.sqrt(
            priority["n_DAR_tested_regions"]
            / priority["n_DAR_tested_regions"].max()
        ),
        c=priority["median_logFC"],
        cmap="RdBu_r",
        vmin=-max(abs(priority["median_logFC"]).max(), 1e-6),
        vmax=max(abs(priority["median_logFC"]).max(), 1e-6),
        edgecolor="white",
        linewidth=0.7,
    )
    axis.axhline(0.5, color="#777777", linewidth=0.7, linestyle="--")
    axis.set_xticks(x)
    axis.set_xticklabels(priority["display"], rotation=35, ha="right", fontsize=8)
    axis.set_ylim(0, 1)
    axis.set_ylabel("Fraction of tested regions with logFC > 0")
    axis.set_title("Direction and coverage of accessibility shifts")
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.045, pad=0.03)
    colorbar.set_label("Median region logFC", fontsize=7)
    colorbar.ax.tick_params(labelsize=6)
    axis.text(
        0.02,
        0.98,
        "Point area = tested regions",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    axis.spines[["top", "right"]].set_visible(False)
    panel_label(axis, "D")

    figure.suptitle(
        "Condition-aware accessibility support for the archived skin eRegulons",
        fontsize=13,
        y=0.995,
    )
    figure.text(
        0.5,
        0.018,
        "Shared-source condition-aware analysis; eRegulons and DAR both derive from GSE312129.",
        ha="center",
        fontsize=7,
        color="#4A4A4A",
    )
    figure.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.10)
    pdf_path = args.outdir / "Figure4_condition_aware_candidate.pdf"
    png_path = args.outdir / "Figure4_condition_aware_candidate.png"
    figure.savefig(pdf_path, dpi=600, bbox_inches="tight")
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(pdf_path)


if __name__ == "__main__":
    main()
