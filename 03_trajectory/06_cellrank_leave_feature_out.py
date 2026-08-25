#!/usr/bin/env python3
"""Rebuild CellRank graphs under two prespecified HES1 feature ablations.

The graph is independently rebuilt from the full per-cohort expression matrices
for baseline, HES1 removal, and fixed CollecTRI HES1-target removal. HES1
activity, donor identities, cells, root rule, terminal-state rule, and graph
parameters are held fixed. This is a sensitivity analysis, not an independent
validation or a new outcome-optimized trajectory.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
from importlib.metadata import version as package_version
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr, wilcoxon


SEED = 0
N_HVG = 2000
N_PCS = 50
N_NEIGHBORS = 30
N_DCS = 15
HARD_FRAC_TO_KEEP = 0.3
FATE_N_JOBS = 1
N_SCHUR = 20
N_MACROSTATES = 10
MIN_DONOR_CELLS = 30
N_BOOT = 10_000

ROOT_MARKERS = ["SFRP2", "DPP4", "PI16"]
TERM_MARKERS = ["CTHRC1", "SFRP4", "ACTA2", "COMP", "POSTN", "COL11A1"]
MYO_MARKERS = [
    "ACTA2", "TAGLN", "POSTN", "COL1A1", "COL1A2",
    "COMP", "COL11A1", "CTHRC1", "COL3A1", "FN1",
]
ACTIVITY_COLUMNS = ["HES1_regulon", "FOSB_regulon", "SMAD3_regulon", "MEF2C_regulon"]
BASELINE_ARM = "baseline"
REMOVE_HES1_ARM = "remove_hes1"
REMOVE_TARGETS_ARM = "remove_collectri_hes1_targets"
ABLATION_ARMS = [REMOVE_HES1_ARM, REMOVE_TARGETS_ARM]
ARM_LABELS = {
    REMOVE_HES1_ARM: "Remove HES1",
    REMOVE_TARGETS_ARM: "Remove fixed CollecTRI HES1 targets",
}

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCES = [
    PROJECT / "00_pipeline_current/local_annotation/GSE181957/GSE181957_annotated.h5ad",
    PROJECT / "00_pipeline_current/local_annotation/GSE236111/GSE236111_annotated.h5ad",
    PROJECT / "00_pipeline_current/local_annotation/GSE249279/GSE249279_annotated.h5ad",
    PROJECT / "00_pipeline_current/local_annotation/GSE292979/GSE292979_annotated.h5ad",
    PROJECT / "00_pipeline_current/local_annotation/gur/gur_annotated.h5ad",
    PROJECT / "00_pipeline_current/local_annotation/Tabib_GSE138669/Tabib_GSE138669_annotated.h5ad",
]
DEFAULT_ACTIVITY = (
    PROJECT / "03_results_figures/plot_data_local/powered/out_04_trajectory/trajectory_pseudotime.csv"
)
DEFAULT_DONOR_MAP = (
    PROJECT / "04_manuscript/revision_20260722/donor_metadata/outputs/sample_to_donor.csv"
)
DEFAULT_COLLECTRI = PROJECT / "02_analysis_modules/rigor_fixes/collectri_net.csv"
DEFAULT_OUT = PROJECT / "04_manuscript/revision_20260722/figure3/leave_feature_out"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-h5ad", action="append", type=Path, dest="sources")
    parser.add_argument("--activity-csv", type=Path, default=DEFAULT_ACTIVITY)
    parser.add_argument("--donor-map", type=Path, default=DEFAULT_DONOR_MAP)
    parser.add_argument("--collectri-net", type=Path, default=DEFAULT_COLLECTRI)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--max-cells-per-donor",
        type=int,
        default=0,
        help="Outcome-blind deterministic cap; 0 retains every matched cell.",
    )
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, with_hash: bool = False) -> dict[str, object]:
    stat = path.stat()
    record: dict[str, object] = {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if with_hash:
        record["sha256"] = sha256(path)
    return record


def register_anndata_null_reader() -> None:
    """Read AnnData 0.12 ``null`` metadata without altering the H5AD."""
    major_minor = tuple(
        int(value) for value in package_version("anndata").split(".")[:2]
    )
    if major_minor >= (0, 12):
        return
    from h5py import Dataset

    from anndata._io.specs import IOSpec, _REGISTRY

    @_REGISTRY.register_read(Dataset, IOSpec("null", "0.1.0"))
    def read_null(_elem, _reader):
        return None


def p_adjust(values: np.ndarray, method: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    if method == "bh":
        adjusted = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    elif method == "holm":
        adjusted = np.maximum.accumulate(ranked * (m - np.arange(m)))
    else:
        raise ValueError(method)
    restored = np.empty(m)
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return out


def matrix_max(x) -> float:
    if sparse.issparse(x):
        return float(x.max())
    return float(np.nanmax(np.asarray(x)))


def mean_score(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        raise RuntimeError(f"None of the fixed marker genes are present: {genes}")
    values = adata[:, present].X.mean(axis=1)
    return np.asarray(values).ravel().astype(float), present


def load_activity(path: Path) -> pd.DataFrame:
    activity = pd.read_csv(path)
    id_col = "cell_id" if "cell_id" in activity.columns else activity.columns[0]
    missing = [col for col in ACTIVITY_COLUMNS if col not in activity.columns]
    if missing:
        raise RuntimeError(f"Activity file is missing columns: {missing}")
    activity = activity.rename(columns={id_col: "cell_id"}).set_index("cell_id")
    if activity.index.has_duplicates:
        raise RuntimeError("Activity cell identifiers are not unique")
    return activity[ACTIVITY_COLUMNS]


def load_sources(paths: list[Path], activity_ids: set[str]) -> ad.AnnData:
    register_anndata_null_reader()
    parts: list[ad.AnnData] = []
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
        source = ad.read_h5ad(path, backed="r")
        if source.raw is None:
            source.file.close()
            raise RuntimeError(f"{path} has no full-gene .raw matrix")
        if not {"celltype", "cohort", "sample", "fib_subtype"}.issubset(source.obs.columns):
            source.file.close()
            raise RuntimeError(f"{path} is missing required cell metadata")
        keep = source.obs["celltype"].astype(str).eq("Fibroblast")
        keep &= source.obs_names.isin(activity_ids)
        idx = np.flatnonzero(np.asarray(keep))
        if not len(idx):
            source.file.close()
            continue
        x = source.raw.X[idx, :]
        part = ad.AnnData(
            X=x,
            obs=source.obs.iloc[idx].copy(),
            var=source.raw.var.copy(),
        )
        parts.append(part)
        print(f"loaded {path.name}: {part.n_obs:,} matched fibroblasts x {part.n_vars:,} genes")
        source.file.close()
    if not parts:
        raise RuntimeError("No source cells matched the archived trajectory cell identifiers")
    combined = ad.concat(parts, axis=0, join="inner", merge="same", index_unique=None)
    if combined.obs_names.has_duplicates:
        duplicate = combined.obs_names[combined.obs_names.duplicated()].unique()[:10]
        raise RuntimeError(f"Duplicated cell identifiers after concatenation: {list(duplicate)}")
    return combined


def attach_donor_metadata(adata: ad.AnnData, donor_path: Path) -> ad.AnnData:
    donor = pd.read_csv(donor_path)
    columns = ["sample", "cohort", "condition", "donor_id", "analysis_keep", "exclusion_reason"]
    missing = [column for column in columns if column not in donor.columns]
    if missing:
        raise RuntimeError(f"Donor map is missing columns: {missing}")
    donor = donor[columns].drop_duplicates(["sample", "cohort"])
    if donor.duplicated(["sample", "cohort"]).any():
        raise RuntimeError("Donor map is not one-to-one for sample and cohort")
    obs = adata.obs.reset_index(names="cell_id")
    obs["sample"] = obs["sample"].astype(str)
    obs["cohort"] = obs["cohort"].astype(str)
    donor["sample"] = donor["sample"].astype(str)
    donor["cohort"] = donor["cohort"].astype(str)
    obs = obs.merge(donor, on=["sample", "cohort"], how="left", validate="many_to_one")
    if obs["analysis_keep"].isna().any():
        examples = obs.loc[obs["analysis_keep"].isna(), ["sample", "cohort"]].drop_duplicates().head()
        raise RuntimeError(f"Cells lack donor-map rows:\n{examples}")
    keep = obs["analysis_keep"].astype(bool).to_numpy()
    retained = adata[keep].copy()
    retained.obs = obs.loc[keep].set_index("cell_id").loc[retained.obs_names].copy()
    return retained


def deterministic_cap(adata: ad.AnnData, max_cells: int) -> ad.AnnData:
    if max_cells <= 0:
        return adata
    rng = np.random.default_rng(SEED)
    keep: list[str] = []
    for _, frame in adata.obs.groupby("donor_id", observed=True):
        ids = frame.index.to_numpy()
        if len(ids) > max_cells:
            ids = np.sort(rng.choice(ids, size=max_cells, replace=False))
        keep.extend(ids.tolist())
    return adata[keep].copy()


def exclusion_table(net_path: Path, var_names: pd.Index) -> pd.DataFrame:
    net = pd.read_csv(net_path)
    if not {"source", "target"}.issubset(net.columns):
        raise RuntimeError("CollecTRI table must contain source and target columns")
    targets = sorted(set(net.loc[net["source"].astype(str).eq("HES1"), "target"].astype(str)))
    genes = ["HES1"] + [gene for gene in targets if gene != "HES1"]
    return pd.DataFrame(
        {
            "gene": genes,
            "role": ["TF"] + ["fixed_CollecTRI_target"] * (len(genes) - 1),
            "present_in_joint_matrix": [gene in var_names for gene in genes],
        }
    )


def select_root(adata: ad.AnnData) -> str:
    state = adata.obs["fib_subtype"].astype(str)
    preferred = state.str.contains("SFRP2", case=False, na=False) & state.str.contains(
        "DPP4", case=False, na=False
    )
    candidates = np.flatnonzero(preferred.to_numpy())
    if not len(candidates):
        candidates = np.arange(adata.n_obs)
    criterion = adata.obs["root_score"].to_numpy() - adata.obs["term_score"].to_numpy()
    return str(adata.obs_names[candidates[np.argmax(criterion[candidates])]])


def cellrank_fate(adata: ad.AnnData) -> tuple[np.ndarray, dict[str, object]]:
    from cellrank.estimators import GPCCA
    from cellrank.kernels import PseudotimeKernel

    kernel = PseudotimeKernel(adata, time_key="pseudotime")
    kernel.compute_transition_matrix(
        threshold_scheme="hard",
        frac_to_keep=HARD_FRAC_TO_KEEP,
        n_jobs=-1,
    )
    estimator = GPCCA(kernel)
    estimator.compute_schur(n_components=N_SCHUR)
    estimator.compute_macrostates(n_states=N_MACROSTATES)
    macrostates = estimator.macrostates
    labels = np.asarray(macrostates.astype(str))
    valid = labels != "nan"
    if not valid.any():
        raise RuntimeError("GPCCA assigned no macrostate cells")
    by_state = (
        pd.Series(adata.obs["myo"].to_numpy()[valid])
        .groupby(labels[valid], observed=True)
        .mean()
        .sort_values(ascending=False)
    )
    myo_state = str(by_state.index[0])
    estimator.predict_terminal_states(method="top_n", n_states=3)
    predicted = list(map(str, estimator.terminal_states.cat.categories))
    terminal = list(dict.fromkeys([myo_state] + predicted))
    estimator.set_terminal_states(states=terminal)
    terminal = list(map(str, estimator.terminal_states.cat.categories))
    estimator.compute_fate_probabilities(
        solver="gmres",
        use_petsc=True,
        n_jobs=FATE_N_JOBS,
        tol=1e-6,
    )
    probabilities = estimator.fate_probabilities
    probability_matrix = np.asarray(probabilities.X)
    if not np.isfinite(probability_matrix).all():
        raise RuntimeError("CellRank returned non-finite fate probabilities")
    maximum_row_sum_error = float(
        np.max(np.abs(probability_matrix.sum(axis=1) - 1.0))
    )
    if maximum_row_sum_error > 1e-4:
        raise RuntimeError(
            "CellRank fate probabilities do not sum to one; "
            f"maximum error={maximum_row_sum_error:.3g}"
        )
    lineage_names = list(probabilities.names)
    myo_lineage = next((name for name in lineage_names if str(name).startswith(myo_state)), None)
    if myo_lineage is None:
        means = []
        for name in lineage_names:
            values = np.asarray(probabilities[name].X).ravel()
            high = values > np.quantile(values, 0.9)
            means.append(adata.obs.loc[high, "myo"].mean() if high.any() else -np.inf)
        myo_lineage = lineage_names[int(np.argmax(means))]
    fate = np.asarray(probabilities[myo_lineage].X).ravel()
    info = {
        "myo_macrostate": myo_state,
        "myo_lineage": str(myo_lineage),
        "terminal_states": terminal,
        "fate_probability_maximum_row_sum_error": maximum_row_sum_error,
    }
    return fate, info


def run_arm(base: ad.AnnData, features: list[str], arm: str, root_cell: str) -> tuple[pd.DataFrame, dict]:
    try:
        import harmonypy  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "harmonypy is required to rebuild the batch-integrated graph; refusing to use unintegrated PCA"
        ) from exc
    import scanpy.external as sce

    adata = base[:, features].copy()
    sc.pp.pca(adata, n_comps=min(N_PCS, len(features) - 1), random_state=SEED)
    sce.pp.harmony_integrate(
        adata,
        key="cohort",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        random_state=SEED,
        max_iter_harmony=20,
    )
    sc.pp.neighbors(
        adata,
        use_rep="X_pca_harmony",
        n_neighbors=N_NEIGHBORS,
        random_state=SEED,
    )
    n_components, component_labels = connected_components(
        adata.obsp["connectivities"],
        directed=False,
    )
    component_sizes = np.bincount(component_labels)
    if n_components != 1:
        raise RuntimeError(
            f"{arm}: neighbour graph has {n_components} disconnected components; "
            "refusing fate inference until the graph is repaired"
        )
    sc.tl.diffmap(adata, n_comps=N_DCS)
    adata.uns["iroot"] = int(adata.obs_names.get_loc(root_cell))
    sc.tl.dpt(adata, n_dcs=N_DCS)
    pseudotime = adata.obs["dpt_pseudotime"].to_numpy(dtype=float)
    if not np.isfinite(pseudotime).all() or np.nanstd(pseudotime) == 0:
        raise RuntimeError(f"{arm}: non-finite DPT values; graph must be repaired before inference")
    orientation = spearmanr(pseudotime, adata.obs["myo"].to_numpy()).statistic
    if not np.isfinite(orientation):
        raise RuntimeError(f"{arm}: non-finite pseudotime orientation")
    weak_orientation = bool(abs(orientation) < 0.05)
    if weak_orientation:
        print(f"[diagnostic] {arm}: weak pseudotime-marker orientation rho={orientation:.4f}")
    if orientation < 0:
        pseudotime = 1.0 - pseudotime
        orientation = -orientation
    adata.obs["pseudotime"] = pseudotime
    fate, info = cellrank_fate(adata)
    result = adata.obs[
        ["donor_id", "cohort", "condition", "sample", "fib_subtype", *ACTIVITY_COLUMNS]
    ].copy()
    result.insert(0, "cell_id", adata.obs_names)
    result["arm"] = arm
    result["pseudotime"] = pseudotime
    result["myo_fate"] = fate
    info.update(
        {
            "arm": arm,
            "n_cells": adata.n_obs,
            "n_graph_features": adata.n_vars,
            "n_graph_components": int(n_components),
            "largest_graph_component_fraction": float(component_sizes.max() / adata.n_obs),
            "pseudotime_myo_orientation_rho": float(orientation),
            "weak_orientation_diagnostic": weak_orientation,
        }
    )
    del adata
    gc.collect()
    return result, info


def bootstrap_median(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        draws[i] = np.median(rng.choice(values, size=len(values), replace=True))
    return tuple(np.quantile(draws, [0.025, 0.975]))


def donor_statistics(per_cell: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for (arm, donor_id, cohort, condition), group in per_cell.groupby(
        ["arm", "donor_id", "cohort", "condition"], observed=True
    ):
        if len(group) < MIN_DONOR_CELLS or group["myo_fate"].nunique() < 2:
            continue
        for column in ACTIVITY_COLUMNS:
            if group[column].nunique() < 2:
                rho = np.nan
            else:
                rho = spearmanr(group[column], group["myo_fate"]).statistic
            rows.append(
                {
                    "arm": arm,
                    "donor_id": donor_id,
                    "cohort": cohort,
                    "condition": condition,
                    "activity": column,
                    "n_cells": len(group),
                    "rho": rho,
                }
            )
    donor = pd.DataFrame(rows).dropna(subset=["rho"])
    rng = np.random.default_rng(SEED)
    summary_rows = []
    for (arm, activity), group in donor.groupby(["arm", "activity"], observed=True):
        values = group["rho"].to_numpy(float)
        ci_low, ci_high = bootstrap_median(values, rng)
        p_value = float(wilcoxon(values, alternative="two-sided").pvalue)
        summary_rows.append(
            {
                "arm": arm,
                "activity": activity,
                "n_donors": len(values),
                "median_rho": float(np.median(values)),
                "bootstrap_95_ci_low": ci_low,
                "bootstrap_95_ci_high": ci_high,
                "positive_fraction": float(np.mean(values > 0)),
                "wilcoxon_p": p_value,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["bh_q_within_arm"] = np.nan
    summary["holm_p_within_arm"] = np.nan
    for arm, idx in summary.groupby("arm", observed=True).groups.items():
        pvals = summary.loc[idx, "wilcoxon_p"].to_numpy()
        summary.loc[idx, "bh_q_within_arm"] = p_adjust(pvals, "bh")
        summary.loc[idx, "holm_p_within_arm"] = p_adjust(pvals, "holm")

    paired_parts = []
    comparison_rows = []
    paired_wide = donor.pivot_table(
        index=["donor_id", "cohort", "condition", "activity"],
        columns="arm",
        values="rho",
    )
    for ablation_arm in ABLATION_ARMS:
        paired = paired_wide.dropna(subset=[BASELINE_ARM, ablation_arm]).copy()
        if paired.empty:
            raise RuntimeError(f"No paired donor correlations for {ablation_arm}")
        paired["comparison_arm"] = ablation_arm
        paired["delta_ablation_minus_baseline"] = (
            paired[ablation_arm] - paired[BASELINE_ARM]
        )
        paired_parts.append(paired.reset_index())
        for activity, group in paired.groupby("activity", observed=True):
            values = group["delta_ablation_minus_baseline"].to_numpy(float)
            ci_low, ci_high = bootstrap_median(values, rng)
            p_value = 1.0 if np.allclose(values, 0) else float(wilcoxon(values).pvalue)
            comparison_rows.append(
                {
                    "comparison_arm": ablation_arm,
                    "activity": activity,
                    "n_paired_donors": len(values),
                    "median_delta_rho": float(np.median(values)),
                    "bootstrap_95_ci_low": ci_low,
                    "bootstrap_95_ci_high": ci_high,
                    "paired_wilcoxon_p": p_value,
                }
            )
    paired = pd.concat(paired_parts, ignore_index=True)
    comparison = pd.DataFrame(comparison_rows)
    comparison["bh_q_within_comparison"] = np.nan
    comparison["holm_p_within_comparison"] = np.nan
    for arm, idx in comparison.groupby("comparison_arm", observed=True).groups.items():
        pvals = comparison.loc[idx, "paired_wilcoxon_p"].to_numpy()
        comparison.loc[idx, "bh_q_within_comparison"] = p_adjust(pvals, "bh")
        comparison.loc[idx, "holm_p_within_comparison"] = p_adjust(pvals, "holm")

    loo_rows = []
    for arm in donor["arm"].unique():
        arm_data = donor[donor["arm"] == arm]
        for omitted in sorted(arm_data["cohort"].unique()):
            subset = arm_data[arm_data["cohort"] != omitted]
            for activity, group in subset.groupby("activity", observed=True):
                values = group["rho"].to_numpy(float)
                ci_low, ci_high = bootstrap_median(values, rng)
                loo_rows.append(
                    {
                        "arm": arm,
                        "omitted_cohort": omitted,
                        "activity": activity,
                        "n_donors": len(values),
                        "median_rho": float(np.median(values)),
                        "bootstrap_95_ci_low": ci_low,
                        "bootstrap_95_ci_high": ci_high,
                    }
                )
    loo = pd.DataFrame(loo_rows)
    return donor, summary, comparison, loo


def plot_sensitivity(donor: pd.DataFrame, comparison: pd.DataFrame, outdir: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        len(ABLATION_ARMS),
        2,
        figsize=(7.2, 3.0 * len(ABLATION_ARMS)),
        constrained_layout=True,
        squeeze=False,
    )
    hes1_wide = donor[donor["activity"] == "HES1_regulon"].pivot_table(
        index="donor_id", columns="arm", values="rho"
    )
    for row, ablation_arm in enumerate(ABLATION_ARMS):
        hes1 = hes1_wide.dropna(subset=[BASELINE_ARM, ablation_arm])
        if hes1.empty:
            raise RuntimeError(
                f"No paired HES1 donor correlations are available for {ablation_arm}"
            )
        ax = axes[row, 0]
        ax.scatter(
            hes1[BASELINE_ARM],
            hes1[ablation_arm],
            s=14,
            alpha=0.55,
            color="#0072B2",
            edgecolor="white",
            linewidth=0.3,
        )
        limits = [
            float(min(hes1[[BASELINE_ARM, ablation_arm]].min().min(), -1.0)),
            float(max(hes1[[BASELINE_ARM, ablation_arm]].max().max(), 1.0)),
        ]
        ax.plot(limits, limits, color="#666666", lw=0.8, ls="--")
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_xlabel("Baseline donor HES1--fate rho")
        ax.set_ylabel(f"{ARM_LABELS[ablation_arm]} donor rho")
        ax.set_title(f"Paired deposited donors (n={len(hes1)})")

        data = comparison[comparison["comparison_arm"] == ablation_arm].sort_values(
            "median_delta_rho"
        )
        order = list(data["activity"])
        data = data.set_index("activity").loc[order]
        y = np.arange(len(data))
        median = data["median_delta_rho"].to_numpy(float)
        low = data["bootstrap_95_ci_low"].to_numpy(float)
        high = data["bootstrap_95_ci_high"].to_numpy(float)
        ax = axes[row, 1]
        ax.errorbar(
            median,
            y,
            xerr=np.vstack([median - low, high - median]),
            fmt="D",
            color="#D55E00",
            ecolor="#D55E00",
            capsize=2.5,
            ms=4.5,
            lw=1.1,
        )
        ax.axvline(0, color="#666666", lw=0.8, ls="--")
        ax.set_yticks(y, [value.replace("_regulon", "") for value in order])
        ax.set_xlabel("Median donor rho change\n(ablation minus baseline)")
        ax.set_title(f"{ARM_LABELS[ablation_arm]}: bootstrap 95% CI")
    for suffix in ("pdf", "png"):
        fig.savefig(
            outdir / f"Figure3_leave_feature_out_audit.{suffix}",
            dpi=400 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def self_test() -> None:
    p = np.array([0.01, 0.04, 0.03, 0.2])
    assert np.allclose(p_adjust(p, "bh"), [0.04, 0.05333333333333334, 0.05333333333333334, 0.2])
    assert np.allclose(p_adjust(p, "holm"), [0.04, 0.09, 0.09, 0.2])
    print("self-test passed")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    sources = args.sources or DEFAULT_SOURCES
    args.outdir.mkdir(parents=True, exist_ok=True)
    activity = load_activity(args.activity_csv)
    base = load_sources(sources, set(activity.index))
    base = attach_donor_metadata(base, args.donor_map)
    base = deterministic_cap(base, args.max_cells_per_donor)
    activity = activity.loc[base.obs_names]
    for column in ACTIVITY_COLUMNS:
        base.obs[column] = activity[column].to_numpy(float)

    if matrix_max(base.X) > 50:
        sc.pp.normalize_total(base, target_sum=1e4)
        sc.pp.log1p(base)
    base.obs["root_score"], root_genes = mean_score(base, ROOT_MARKERS)
    base.obs["term_score"], term_genes = mean_score(base, TERM_MARKERS)
    base.obs["myo"], myo_genes = mean_score(base, MYO_MARKERS)
    root_cell = select_root(base)

    sc.pp.highly_variable_genes(
        base,
        n_top_genes=N_HVG,
        flavor="seurat",
        batch_key="cohort",
        subset=False,
    )
    baseline_features = list(base.var_names[base.var["highly_variable"].to_numpy(bool)])
    excluded = exclusion_table(args.collectri_net, base.var_names)
    excluded.to_csv(args.outdir / "Figure3_leave_feature_out_genes.csv", index=False)
    present_exclusions = list(
        excluded.loc[excluded["present_in_joint_matrix"], "gene"].astype(str)
    )
    baseline_features = list(dict.fromkeys(baseline_features + present_exclusions))
    hes1_set = {"HES1"} & set(present_exclusions)
    target_set = set(
        excluded.loc[
            excluded["present_in_joint_matrix"]
            & excluded["role"].eq("fixed_CollecTRI_target"),
            "gene",
        ].astype(str)
    )
    if not hes1_set:
        raise RuntimeError("HES1 is absent from the joint fibroblast matrix")
    if not target_set:
        raise RuntimeError("No fixed CollecTRI HES1 targets are present in the joint matrix")
    arm_features = {
        BASELINE_ARM: baseline_features,
        REMOVE_HES1_ARM: [gene for gene in baseline_features if gene not in hes1_set],
        REMOVE_TARGETS_ARM: [gene for gene in baseline_features if gene not in target_set],
    }
    for arm in ABLATION_ARMS:
        if not arm_features[arm] or len(arm_features[arm]) == len(baseline_features):
            raise RuntimeError(f"{arm}: prespecified genes were not removed")

    outputs = []
    arm_info = []
    for arm, features in arm_features.items():
        print(f"running {arm}: {base.n_obs:,} cells x {len(features):,} fixed graph features")
        per_cell, info = run_arm(base, features, arm, root_cell)
        outputs.append(per_cell)
        arm_info.append(info)
        per_cell.to_csv(args.outdir / f"Figure3_{arm}_per_cell.csv.gz", index=False)

    combined = pd.concat(outputs, ignore_index=True)
    donor, summary, comparison, loo = donor_statistics(combined)
    donor.to_csv(args.outdir / "Figure3_leave_feature_out_donor_rho.csv", index=False)
    summary.to_csv(args.outdir / "Figure3_leave_feature_out_summary.csv", index=False)
    comparison.to_csv(args.outdir / "Figure3_leave_feature_out_paired_comparison.csv", index=False)
    loo.to_csv(args.outdir / "Figure3_leave_feature_out_leave_one_cohort_out.csv", index=False)
    plot_sensitivity(donor, comparison, args.outdir)

    manifest = {
        "analysis_unit": "deposited donor",
        "purpose": "CellRank graph leave-feature-out sensitivity",
        "script_sha256": sha256(Path(__file__).resolve()),
        "seed": SEED,
        "parameters": {
            "n_hvg_requested": N_HVG,
            "baseline_feature_rule": (
                "cohort-aware HVGs plus HES1 and all fixed CollecTRI HES1 targets "
                "present in the joint fibroblast matrix"
            ),
            "n_pcs": N_PCS,
            "integration": "Harmony by cohort in each arm",
            "n_neighbors": N_NEIGHBORS,
            "n_diffusion_components": N_DCS,
            "cellrank_kernel": "PseudotimeKernel hard threshold",
            "cellrank_hard_frac_to_keep": HARD_FRAC_TO_KEEP,
            "cellrank_transition_jobs": -1,
            "fate_solver": "PETSc GMRES, tolerance 1e-6",
            "fate_solver_jobs": FATE_N_JOBS,
            "terminal_state_rule": (
                "three GPCCA top-n terminal macrostates with the highest fixed "
                "myofibroblast-marker macrostate forced into the set"
            ),
            "n_schur": N_SCHUR,
            "n_macrostates": N_MACROSTATES,
            "min_cells_per_donor": MIN_DONOR_CELLS,
            "donor_bootstrap_replicates": N_BOOT,
            "max_cells_per_donor": args.max_cells_per_donor,
        },
        "fixed_anchors": {
            "root_cell": root_cell,
            "root_markers_present": root_genes,
            "terminal_markers_present": term_genes,
            "myofibroblast_markers_present": myo_genes,
            "activity_columns_reused_from_archived_pre-ablation_output": ACTIVITY_COLUMNS,
        },
        "exclusion": {
            "definition": "HES1 and fixed CollecTRI HES1 targets are ablated in separate arms",
            "n_network_genes": len(excluded),
            "n_present_in_joint_matrix": int(excluded["present_in_joint_matrix"].sum()),
            "n_removed_remove_hes1": len(baseline_features)
            - len(arm_features[REMOVE_HES1_ARM]),
            "n_removed_remove_collectri_hes1_targets": len(baseline_features)
            - len(arm_features[REMOVE_TARGETS_ARM]),
        },
        "arms": arm_info,
        "software": {
            "python": platform.python_version(),
            "anndata": ad.__version__,
            "cellrank": package_version("cellrank"),
            "harmonypy": package_version("harmonypy"),
            "scanpy": sc.__version__,
            "scvelo": package_version("scvelo"),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": package_version("scipy"),
        },
        "inputs": [file_record(path, with_hash=True) for path in sources]
        + [
            file_record(args.activity_csv, with_hash=True),
            file_record(args.donor_map, with_hash=True),
            file_record(args.collectri_net, with_hash=True),
        ],
        "interpretation_boundary": (
            "These separate rebuilds test graph-feature leakage from HES1 and its fixed "
            "external CollecTRI targets. They remain observational and do not establish "
            "HES1 causality."
        ),
        "outputs": {},
    }
    output_paths = [
        args.outdir / "Figure3_baseline_per_cell.csv.gz",
        args.outdir / "Figure3_remove_hes1_per_cell.csv.gz",
        args.outdir / "Figure3_remove_collectri_hes1_targets_per_cell.csv.gz",
        args.outdir / "Figure3_leave_feature_out_genes.csv",
        args.outdir / "Figure3_leave_feature_out_donor_rho.csv",
        args.outdir / "Figure3_leave_feature_out_summary.csv",
        args.outdir / "Figure3_leave_feature_out_paired_comparison.csv",
        args.outdir / "Figure3_leave_feature_out_leave_one_cohort_out.csv",
        args.outdir / "Figure3_leave_feature_out_audit.pdf",
        args.outdir / "Figure3_leave_feature_out_audit.png",
    ]
    manifest["outputs"] = {
        str(path): sha256(path) for path in output_paths if path.exists()
    }
    (args.outdir / "Figure3_leave_feature_out_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"completed: {args.outdir}")


if __name__ == "__main__":
    main()
