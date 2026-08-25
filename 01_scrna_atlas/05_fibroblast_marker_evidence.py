#!/usr/bin/env python
"""Donor-balanced marker evidence for the seven atlas fibroblast states.

The compact figure atlas contains only 129 genes. To evaluate the deposited
subtype labels with their complete literature-derived marker panels, exact cell
identifiers are matched back to the seven locally retained full-gene cohort
objects. Expression and detection are first aggregated within donor and state;
the dot plot then gives each donor equal weight.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "server_archive" / "fig_atlas" / "atlas_fig.h5ad"
DONOR_MAP = ROOT / "04_manuscript" / "revision_20260722" / "donor_metadata" / "outputs" / "sample_to_donor.csv"
ANNOTATION_ROOT = ROOT / "00_pipeline_current" / "local_annotation"
OUT = ROOT / "04_manuscript" / "revision_20260722" / "figure1" / "outputs"

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
    "Myofibroblast": "Myofibroblast",
    "SFRP4_proFib": "SFRP4+ pro-fibrotic",
    "SFRP2_DPP4": "SFRP2+/DPP4+",
    "Adipogenic": "Adipogenic",
    "FMO1_LSP1": "FMO1+/LSP1+",
    "LGR5_Gur": "LGR5+",
    "Inflammatory": "Inflammatory-like*",
}
MARKER_PANELS = {
    "Myofibroblast": ["ACTA2", "TAGLN", "POSTN", "COL11A1", "COMP"],
    "SFRP4_proFib": ["SFRP4", "PRSS23", "TNC"],
    "SFRP2_DPP4": ["SFRP2", "DPP4", "PI16"],
    "Adipogenic": ["APOE", "CFD", "GPX3", "PPARG"],
    "FMO1_LSP1": ["FMO1", "LSP1", "APOD", "CD34"],
    "LGR5_Gur": ["LGR5", "WIF1", "SFRP2"],
    "Inflammatory": ["CCL19", "CXCL12", "CD74", "HLA-DRA"],
}
PANEL_COLORS = {
    "Myofibroblast": "#D55E00",
    "SFRP4_proFib": "#E69F00",
    "SFRP2_DPP4": "#0072B2",
    "Adipogenic": "#CC79A7",
    "FMO1_LSP1": "#009E73",
    "LGR5_Gur": "#56B4E9",
    "Inflammatory": "#777777",
}


def log(message: str) -> None:
    print(f"[figure1-marker] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_exact_labels() -> tuple[pd.DataFrame, dict]:
    donor_map = pd.read_csv(DONOR_MAP, dtype={"sample": str, "donor_id": str})
    donor_map["analysis_keep"] = donor_map["analysis_keep"].astype(str).str.lower().eq("true")
    kept = donor_map[donor_map["analysis_keep"]][
        ["sample", "donor_id", "cohort", "condition"]
    ]
    atlas = ad.read_h5ad(ATLAS, backed="r")
    labels = atlas.obs[["sample", "celltype", "fib_subtype"]].copy().reset_index(names="cell_id")
    atlas.file.close()
    labels["sample"] = labels["sample"].astype(str)
    labels["celltype"] = labels["celltype"].astype(str)
    labels["fib_subtype"] = labels["fib_subtype"].astype(str)
    labels = labels.merge(kept, on="sample", how="inner", validate="many_to_one")
    labels = labels[
        labels["celltype"].eq("Fibroblast") & labels["fib_subtype"].isin(SUBTYPES)
    ].copy()
    labels = labels.set_index("cell_id", verify_integrity=True)
    audit = {
        "exact_label_cells": int(len(labels)),
        "exact_label_donors": int(labels["donor_id"].nunique()),
        "state_cell_counts": labels["fib_subtype"].value_counts().to_dict(),
        "state_donor_counts": labels.groupby("fib_subtype")["donor_id"].nunique().to_dict(),
    }
    return labels, audit


def aggregate_full_gene_expression(labels: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    genes = list(dict.fromkeys(gene for panel in MARKER_PANELS.values() for gene in panel))
    rows: list[dict] = []
    input_files: list[dict] = []
    files = sorted(ANNOTATION_ROOT.glob("*/*_annotated.h5ad"))
    for path in files:
        adata = ad.read_h5ad(path, backed="r")
        common = adata.obs_names.intersection(labels.index)
        if len(common) == 0:
            adata.file.close()
            continue
        if adata.raw is None:
            adata.file.close()
            raise ValueError(f"Full-gene raw matrix is absent: {path}")
        missing = [gene for gene in genes if gene not in adata.raw.var_names]
        if missing:
            adata.file.close()
            raise ValueError(f"Required marker genes missing from {path}: {missing}")
        positions = adata.obs_names.get_indexer(common)
        gene_positions = [adata.raw.var_names.get_loc(gene) for gene in genes]
        matrix = adata.raw.X[positions, :][:, gene_positions].tocsr()
        cell_meta = labels.loc[common, ["donor_id", "cohort", "condition", "fib_subtype"]].reset_index(drop=True)
        grouped = cell_meta.groupby(["donor_id", "cohort", "condition", "fib_subtype"], observed=True).indices
        for (donor_id, cohort, condition, subtype), idx in grouped.items():
            group_matrix = matrix[np.asarray(idx), :]
            sums = np.asarray(group_matrix.sum(axis=0)).ravel()
            positives = np.asarray((group_matrix > 0).sum(axis=0)).ravel()
            n_cells = int(len(idx))
            for gene_i, gene in enumerate(genes):
                rows.append(
                    {
                        "donor_id": donor_id,
                        "cohort": cohort,
                        "condition": condition,
                        "subtype": subtype,
                        "gene": gene,
                        "n_cells": n_cells,
                        "sum_expression": float(sums[gene_i]),
                        "positive_cells": int(positives[gene_i]),
                    }
                )
        input_files.append(
            {
                "path": str(path),
                "size_bytes": int(path.stat().st_size),
                "matched_cells": int(len(common)),
            }
        )
        log(f"{path.parent.name}: matched {len(common):,} cells")
        adata.file.close()
    if not rows:
        raise ValueError("No exact cell identifiers matched the full-gene cohort objects")

    donor = pd.DataFrame(rows)
    donor = (
        donor.groupby(["donor_id", "cohort", "condition", "subtype", "gene"], observed=True)
        .agg(
            n_cells=("n_cells", "sum"),
            sum_expression=("sum_expression", "sum"),
            positive_cells=("positive_cells", "sum"),
        )
        .reset_index()
    )
    donor["mean_expression"] = donor["sum_expression"] / donor["n_cells"]
    donor["fraction_expressing"] = donor["positive_cells"] / donor["n_cells"]
    return donor, input_files


def summarize(donor: pd.DataFrame) -> pd.DataFrame:
    summary = (
        donor.groupby(["subtype", "gene"], observed=True)
        .agg(
            n_donors=("donor_id", "nunique"),
            n_cells=("n_cells", "sum"),
            donor_mean_expression=("mean_expression", "mean"),
            donor_median_expression=("mean_expression", "median"),
            donor_mean_fraction_expressing=("fraction_expressing", "mean"),
            donor_median_fraction_expressing=("fraction_expressing", "median"),
        )
        .reset_index()
    )
    summary["expression_z"] = summary.groupby("gene")["donor_mean_expression"].transform(
        lambda x: (x - x.mean()) / (x.std(ddof=0) if x.std(ddof=0) > 0 else 1.0)
    )
    return summary


def plot_marker_dotplot(summary: pd.DataFrame, donor: pd.DataFrame) -> None:
    columns = [(panel, gene) for panel, genes in MARKER_PANELS.items() for gene in genes]
    state_counts = (
        donor.drop_duplicates(["donor_id", "subtype"])
        .groupby("subtype", observed=True)
        .agg(n_donors=("donor_id", "nunique"), n_cells=("n_cells", "sum"))
    )
    lookup = summary.set_index(["subtype", "gene"])
    fig, ax = plt.subplots(figsize=(12.4, 4.7), constrained_layout=True)
    last_scatter = None
    for y, subtype in enumerate(SUBTYPES[::-1]):
        for x, (_, gene) in enumerate(columns):
            row = lookup.loc[(subtype, gene)]
            fraction = float(row["donor_mean_fraction_expressing"])
            last_scatter = ax.scatter(
                x, y, s=10 + 180 * fraction, c=float(row["expression_z"]),
                cmap="RdBu_r", vmin=-2.2, vmax=2.2, edgecolor="#666666", linewidth=0.35,
            )

    boundaries = np.cumsum([len(genes) for genes in MARKER_PANELS.values()])
    for boundary in boundaries[:-1]:
        ax.axvline(boundary - 0.5, color="#BDBDBD", linewidth=0.7)
    start = 0
    for panel, genes in MARKER_PANELS.items():
        stop = start + len(genes)
        ax.plot(
            [start - 0.42, stop - 0.58], [1.035, 1.035], transform=ax.get_xaxis_transform(),
            color=PANEL_COLORS[panel], linewidth=3.0, clip_on=False,
        )
        ax.text(
            (start + stop - 1) / 2, 1.075, DISPLAY[panel],
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=6.8, color=PANEL_COLORS[panel], fontweight="bold", clip_on=False,
        )
        start = stop

    ylabels = []
    for subtype in SUBTYPES[::-1]:
        counts = state_counts.loc[subtype]
        ylabels.append(f"{DISPLAY[subtype]}  (n={int(counts.n_donors)} donors)")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([gene for _, gene in columns], rotation=55, ha="right", fontsize=7, fontstyle="italic")
    ax.set_yticks(range(len(SUBTYPES)))
    ax.set_yticklabels(ylabels, fontsize=7.4)
    ax.set_xlim(-0.7, len(columns) - 0.3)
    ax.set_ylim(-0.65, len(SUBTYPES) - 0.35)
    ax.grid(color="#E8E8E8", linewidth=0.55, linestyle=":")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(
        "Donor-balanced marker evidence for atlas fibroblast states",
        fontsize=10.5, fontweight="bold", pad=42,
    )
    colorbar = fig.colorbar(last_scatter, ax=ax, fraction=0.018, pad=0.015)
    colorbar.set_label("Mean expression z-score across states", fontsize=7)
    colorbar.ax.tick_params(labelsize=6.5)
    size_handles = [
        ax.scatter([], [], s=10 + 180 * fraction, color="#BDBDBD", edgecolor="#666666", linewidth=0.35)
        for fraction in [0.1, 0.5, 0.9]
    ]
    ax.legend(
        size_handles, ["10%", "50%", "90%"], title="Mean donor-level\ndetection",
        loc="upper left", bbox_to_anchor=(1.015, 0.72), frameon=False,
        fontsize=6.5, title_fontsize=6.8,
    )
    fig.text(
        0.01, 0.005,
        "* Inflammatory-like is a sparse state (54 cells from one HC donor) and is shown descriptively.",
        fontsize=7, ha="left", va="bottom",
    )
    png = OUT / "Figure1_marker_fibroblast_states.png"
    pdf = OUT / "Figure1_marker_fibroblast_states.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    Image.open(png).convert("L").save(OUT / "Figure1_marker_fibroblast_states_grayscale.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels, audit = load_exact_labels()
    donor, input_files = aggregate_full_gene_expression(labels)
    summary = summarize(donor)
    donor.to_csv(OUT / "Figure1_marker_donor_source_data.csv", index=False)
    summary.to_csv(OUT / "Figure1_marker_summary.csv", index=False)
    pd.DataFrame(
        [
            {"subtype": subtype, "display": DISPLAY[subtype], "gene": gene, "panel_order": order}
            for subtype, genes in MARKER_PANELS.items()
            for order, gene in enumerate(genes, start=1)
        ]
    ).to_csv(OUT / "Figure1_marker_panel.csv", index=False)
    support = (
        labels.groupby(["fib_subtype", "cohort"], observed=True)
        .agg(n_cells=("sample", "size"), n_donors=("donor_id", "nunique"))
        .reset_index()
    )
    support.to_csv(OUT / "Figure1_subtype_support_by_cohort.csv", index=False)
    plot_marker_dotplot(summary, donor)
    manifest = {
        **audit,
        "aggregation": "exact cell-id match; donor-state aggregation before equal-donor summary",
        "expression_source": "full-gene raw log-normalized matrices from local annotated cohort H5AD files",
        "marker_panels": MARKER_PANELS,
        "atlas_sha256": sha256(ATLAS),
        "donor_map_sha256": sha256(DONOR_MAP),
        "input_files": input_files,
    }
    (OUT / "Figure1_marker_run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    log(f"completed marker evidence: {OUT}")


if __name__ == "__main__":
    main()
