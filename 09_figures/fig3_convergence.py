"""Draw revised Figure 3 from donor-level and held-out-cohort audit outputs."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figstyle as fs  # noqa: E402


PROJECT = HERE.parents[1]
SOURCE = PROJECT / "04_manuscript/revision_20260722/figure3/outputs"
Q3 = PROJECT / "04_manuscript/plot_data_local/atlas_queries/Q3b_HES1_regulon_activity_by_celltype.csv"
MANUSCRIPT_FIGURES = PROJECT / "04_manuscript/manuscript_final/figures"

TF_ORDER = ["SMAD3", "HES1", "MEF2C"]
MODEL_ORDER = ["hes1_only", "full_panel", "remove_hes1"]
MODEL_LABEL = {"hes1_only": "HES1 only", "full_panel": "Full panel", "remove_hes1": "Panel minus HES1"}
MODEL_COLOR = {"hes1_only": "#6B6B6B", "full_panel": "#0072B2", "remove_hes1": "#D55E00"}
COHORT_ORDER = [
    "GSE130955",
    "GSE58095",
    "GSE249550",
    "GSE181549",
    "GSE9285",
    "GSE32413",
    "GSE76807",
    "GSE231692",
]
COHORT_COLOR = {
    "GSE236111": "#0072B2",
    "GSE249279": "#E69F00",
    "GSE292979": "#009E73",
    "Gur_GSE195452": "#CC79A7",
    "Tabib_GSE138669": "#56B4E9",
    "GSE130955": "#0072B2",
    "GSE58095": "#D55E00",
    "GSE249550": "#009E73",
}

CELLTYPE_LABEL = {
    "Myofibroblast": "Myofibroblast",
    "SFRP4_proFib": "SFRP4+ pro-fibrotic",
    "SFRP2_DPP4": "SFRP2+/DPP4+",
    "Adipogenic": "Adipogenic",
    "FMO1_LSP1": "FMO1+/LSP1+",
    "LGR5_Gur": "LGR5+",
    "Inflammatory": "Inflammatory-like*",
    "Fibroblast_other": "Other fibroblast (unassigned)",
    "Tcell": "T cell",
    "Bcell": "B cell",
    "SmoothMuscle": "Smooth muscle",
}
FIBROBLAST_STATES = {
    "Myofibroblast",
    "SFRP4_proFib",
    "SFRP2_DPP4",
    "Adipogenic",
    "FMO1_LSP1",
    "LGR5_Gur",
    "Inflammatory",
    "Fibroblast_other",
}


def panel_a(ax: plt.Axes) -> None:
    source = pd.read_csv(SOURCE / "Figure3a_donor_activity_source_data.csv")
    result = pd.read_csv(SOURCE / "Figure3a_donor_activity_results.csv").set_index("TF").loc["HES1"]
    data = source.loc[source["TF"].eq("HES1")].copy()
    conditions = ["HC", "SSc"]
    rng = np.random.default_rng(20260723)
    values = [data.loc[data["condition"].eq(condition), "activity"] for condition in conditions]
    boxes = ax.boxplot(
        values,
        positions=[0, 1],
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.4},
        boxprops={"facecolor": "white", "edgecolor": "#555555", "linewidth": 0.9},
        whiskerprops={"color": "#555555", "linewidth": 0.8},
        capprops={"color": "#555555", "linewidth": 0.8},
    )
    for box in boxes["boxes"]:
        box.set_alpha(0.75)
    for x, condition in enumerate(conditions):
        subset = data.loc[data["condition"].eq(condition)]
        for cohort in sorted(subset["cohort"].unique()):
            cohort_data = subset.loc[subset["cohort"].eq(cohort)]
            jitter = rng.uniform(-0.18, 0.18, len(cohort_data))
            ax.scatter(
                x + jitter,
                cohort_data["activity"],
                s=14,
                color=COHORT_COLOR[cohort],
                alpha=0.68,
                edgecolor="white",
                linewidth=0.25,
                label=cohort if x == 0 else None,
                zorder=3,
            )
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"HC\n(n={int(result.n_HC)})", f"SSc\n(n={int(result.n_SSc)})"])
    ax.set_ylabel("HES1 regulon activity\n(ULM t-value)")
    ax.set_title("Donor-level myofibroblast activity", fontsize=9.5)
    ax.text(
        0.02,
        0.98,
        "Cohort-adjusted difference\n"
        f"{result.adjusted_mean_difference:.2f} "
        f"[{result.adjusted_ci_low:.2f}, {result.adjusted_ci_high:.2f}]\n"
        f"BH q={result.q_bh:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2},
    )
    ax.legend(fontsize=5.8, frameon=False, loc="lower right", handletextpad=0.3, labelspacing=0.25)
    fs.despine(ax)
    fs.panel_label(ax, "a")


def panel_b(ax: plt.Axes) -> None:
    source = pd.read_csv(Q3)
    pivot = source.pivot_table(index="TF", columns="ct_work", values="mean").loc[TF_ORDER]
    columns = pivot.loc["HES1"].sort_values(ascending=False).index.tolist()
    pivot = pivot[columns]
    values = pivot.to_numpy()
    vmax = float(np.nanmax(np.abs(values)))
    image = ax.imshow(values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
    ax.set_yticks(range(len(TF_ORDER)))
    ax.set_yticklabels(["SMAD3\n(positive control)", "HES1", "MEF2C\n(negative control)"], fontsize=7.2)
    ax.set_xticks(range(len(columns)))
    labels = [CELLTYPE_LABEL.get(column, column) for column in columns]
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=5.8)
    for tick, column in zip(ax.get_xticklabels(), columns):
        if column == "Myofibroblast":
            tick.set_color("#D55E00")
            tick.set_fontweight("bold")
        elif column in FIBROBLAST_STATES:
            tick.set_color("#0072B2")
    ax.set_title("Regulon activity across annotated cell states", fontsize=9.5)
    colorbar = plt.colorbar(image, ax=ax, fraction=0.020, pad=0.012)
    colorbar.set_label("ULM t-value", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6)
    ax.text(
        0.995,
        0.025,
        "*54 cells from one HC donor; descriptive only",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.6,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    fs.panel_label(ax, "b")


def panel_c(ax: plt.Axes) -> None:
    source = pd.read_csv(SOURCE / "Figure3b_donor_trajectory_source_data.csv")
    source = source.loc[source["n_cells"].ge(30)]
    result = pd.read_csv(SOURCE / "Figure3b_trajectory_results.csv").set_index("TF")
    order = ["SMAD3", "HES1", "MEF2C"]
    colors = {"SMAD3": "#D55E00", "HES1": "#0072B2", "MEF2C": "#999999"}
    rng = np.random.default_rng(20260723)
    positions = np.arange(len(order))[::-1]
    for position, tf in zip(positions, order):
        values = source.loc[source["TF"].eq(tf), "rho"].dropna().to_numpy()
        violin = ax.violinplot(
            [values], positions=[position], orientation="horizontal", widths=0.72, showextrema=False
        )
        for body in violin["bodies"]:
            body.set_facecolor(colors[tf])
            body.set_edgecolor(colors[tf])
            body.set_alpha(0.24)
        jitter = rng.uniform(-0.14, 0.14, len(values))
        ax.scatter(values, position + jitter, s=5, color=colors[tf], alpha=0.36, linewidth=0)
        row = result.loc[tf]
        ax.errorbar(
            row.median_rho,
            position,
            xerr=[[row.median_rho - row.median_rho_ci_low], [row.median_rho_ci_high - row.median_rho]],
            fmt="D",
            markersize=5.5,
            color=colors[tf],
            markeredgecolor="white",
            markeredgewidth=0.7,
            capsize=2.5,
            linewidth=1.2,
            zorder=5,
        )
        ax.text(
            0.72,
            position,
            f"median rho={row.median_rho:+.2f}\n{100 * row.positive_fraction:.0f}% positive",
            ha="left",
            va="center",
            fontsize=6.7,
        )
    ax.axvline(0, color="#555555", linewidth=0.7)
    ax.set_yticks(positions)
    ax.set_yticklabels(["SMAD3\n(positive control)", "HES1", "MEF2C\n(negative control)"], fontsize=7.2)
    ax.set_xlabel("Within-donor Spearman rho vs myofibroblast fate")
    ax.set_title("Trajectory coupling (206 donors; >=30 cells each)", fontsize=9.5)
    ax.set_xlim(-0.55, 1.08)
    fs.despine(ax)
    fs.panel_label(ax, "c")


def panel_d(subgrid) -> None:
    metrics = pd.read_csv(SOURCE / "Figure3d_bulk_ablation_metrics.csv")
    metrics = metrics.loc[metrics["analysis_set"].eq("primary")].copy()
    axes = [plt.subplot(subgrid[0, index]) for index in range(2)]
    positions = np.arange(len(COHORT_ORDER))[::-1]
    offsets = {"hes1_only": 0.20, "full_panel": 0.0, "remove_hes1": -0.20}
    specifications = (
        ("AUROC", "auroc_ci_low", "auroc_ci_high", "Held-out AUROC"),
        (
            "average_precision",
            "average_precision_ci_low",
            "average_precision_ci_high",
            "Held-out average precision",
        ),
    )
    for axis_index, (metric, low, high, title) in enumerate(specifications):
        ax = axes[axis_index]
        for position, cohort in zip(positions, COHORT_ORDER):
            cohort_rows = metrics.loc[metrics["held_out_cohort"].eq(cohort)]
            for model in MODEL_ORDER:
                row = cohort_rows.loc[cohort_rows["model"].eq(model)].iloc[0]
                value = float(row[metric])
                ax.errorbar(
                    value,
                    position + offsets[model],
                    xerr=[
                        [max(0.0, value - float(row[low]))],
                        [max(0.0, float(row[high]) - value)],
                    ],
                    fmt="o",
                    color=MODEL_COLOR[model],
                    markersize=3.8,
                    capsize=1.8,
                    linewidth=0.85,
                    markeredgecolor="white",
                    markeredgewidth=0.35,
                )
            if metric == "average_precision":
                prevalence = float(cohort_rows.iloc[0]["prevalence_SSc"])
                ax.plot(
                    prevalence,
                    position,
                    marker="|",
                    markersize=8,
                    color="#555555",
                    linestyle="none",
                    markeredgewidth=1.1,
                )
        if metric == "AUROC":
            ax.axvline(0.5, color="#888888", linestyle="--", linewidth=0.7)
        ax.set_xlim(-0.02, 1.03)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks(positions)
        ax.set_title(title, fontsize=8.2)
        ax.set_xlabel("Metric (95% bootstrap CI)", fontsize=6.8)
        ax.tick_params(labelsize=6)
        if axis_index == 0:
            ax.set_yticklabels(COHORT_ORDER, fontsize=6.2)
            fs.panel_label(ax, "d")
        else:
            ax.set_yticklabels([])
        for boundary in np.arange(0.5, len(COHORT_ORDER) - 0.5, 1):
            ax.axhline(boundary, color="#EEEEEE", linewidth=0.5, zorder=0)
        fs.despine(ax)
    model_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=MODEL_COLOR[model],
            linestyle="none",
            label=MODEL_LABEL[model],
        )
        for model in MODEL_ORDER
    ]
    model_handles.append(
        plt.Line2D(
            [0],
            [0],
            marker="|",
            color="#555555",
            linestyle="none",
            markersize=8,
            label="SSc prevalence (AP reference)",
        )
    )
    axes[1].legend(
        handles=model_handles,
        fontsize=5.4,
        frameon=False,
        loc="lower right",
        labelspacing=0.25,
        handletextpad=0.3,
    )


def panel_e(subgrid) -> None:
    calibration = pd.read_csv(SOURCE / "Figure3e_calibration_bins.csv")
    metrics = pd.read_csv(SOURCE / "Figure3d_bulk_ablation_metrics.csv")
    calibration = calibration.loc[calibration["analysis_set"].eq("primary")]
    metrics = metrics.loc[metrics["analysis_set"].eq("primary")]
    axes = [
        plt.subplot(subgrid[row, column])
        for row in range(2)
        for column in range(4)
    ]
    for index, (ax, cohort) in enumerate(zip(axes, COHORT_ORDER)):
        ax.plot([0, 1], [0, 1], color="#999999", linestyle="--", linewidth=0.7)
        for model in ("hes1_only", "full_panel", "remove_hes1"):
            data = calibration.loc[
                calibration["held_out_cohort"].eq(cohort) & calibration["model"].eq(model)
            ].sort_values("bin")
            ax.plot(
                data["mean_prediction"],
                data["observed_fraction"],
                marker="o",
                markersize=2.7,
                linewidth=1.0,
                linestyle="--" if model == "remove_hes1" else "-",
                color=MODEL_COLOR[model],
                label=MODEL_LABEL[model],
            )
        full = metrics.loc[metrics["held_out_cohort"].eq(cohort) & metrics["model"].eq("full_panel")].iloc[0]
        ax.set_title(f"{cohort}\nfull Brier={full.brier_score:.2f}", fontsize=6.6)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.tick_params(labelsize=5.5)
        ax.set_xlabel("Predicted", fontsize=6)
        if index % 4 == 0:
            ax.set_ylabel("Observed", fontsize=6)
        else:
            ax.set_yticklabels([])
        if index == 0:
            fs.panel_label(ax, "e")
        fs.despine(ax)
    axes[0].text(
        2.4,
        1.42,
        "Held-out calibration across primary cohorts",
        transform=axes[0].transAxes,
        ha="center",
        va="bottom",
        fontsize=9.5,
    )
    axes[-1].legend(
        fontsize=5.2,
        frameon=False,
        loc="lower right",
        handlelength=1.4,
        labelspacing=0.2,
    )


def main() -> None:
    required = [
        SOURCE / "Figure3a_donor_activity_results.csv",
        SOURCE / "Figure3b_trajectory_results.csv",
        SOURCE / "Figure3d_bulk_ablation_metrics.csv",
        SOURCE / "Figure3e_calibration_bins.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Run 05_figure3_donor_validation.py first:\n" + "\n".join(missing))

    fig = plt.figure(figsize=(15.8, 13.4))
    grid = GridSpec(
        3,
        6,
        figure=fig,
        height_ratios=[0.82, 1.15, 1.02],
        wspace=1.0,
        hspace=0.72,
    )
    panel_a(fig.add_subplot(grid[0, 0:2]))
    panel_b(fig.add_subplot(grid[0, 2:6]))
    panel_c(fig.add_subplot(grid[1, 0:2]))
    panel_d(grid[1, 2:6].subgridspec(1, 2, wspace=0.22))
    panel_e(grid[2, 0:6].subgridspec(2, 4, wspace=0.28, hspace=0.68))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.055)

    output_pdf = SOURCE / "Figure3_revised.pdf"
    output_png = SOURCE / "Figure3_revised.png"
    fig.savefig(output_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    MANUSCRIPT_FIGURES.mkdir(parents=True, exist_ok=True)
    shutil.copy2(output_pdf, MANUSCRIPT_FIGURES / "Fig3.pdf")
    shutil.copy2(output_png, MANUSCRIPT_FIGURES / "Fig3.png")
    print(output_pdf)


if __name__ == "__main__":
    main()
