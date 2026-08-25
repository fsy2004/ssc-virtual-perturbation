"""Rebuild Figure 3 statistics at the deposited-donor and held-out-cohort levels.

This script does not recompute the CellRank graph. It corrects the biological
replicate unit in the existing per-cell trajectory output, and rebuilds the bulk
diagnostic analysis with all preprocessing and model fitting confined to training
cohorts. A separate leave-feature-out CellRank rerun is required to test whether
the trajectory graph changes after excluding HES1 and its regulon genes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import expit, logit
from scipy.stats import binomtest, spearmanr, wilcoxon
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


SEED = 20260723
TF_PANEL = ("HES1", "SMAD3", "MEF2C")
REGULON_COLUMNS = {tf: f"{tf}_regulon" for tf in TF_PANEL}
MODEL_ORDER = ("hes1_only", "full_panel", "remove_hes1")
PRIMARY_COHORTS = (
    "GSE130955",
    "GSE58095",
    "GSE249550",
    "GSE181549",
    "GSE9285",
    "GSE32413",
    "GSE76807",
    "GSE231692",
)
SENSITIVITY_COHORTS = ("GSE95065", "GSE125362", "GSE76885")
ALL_BULK_COHORTS = PRIMARY_COHORTS + SENSITIVITY_COHORTS
ANALYSIS_SETS = {
    "primary": PRIMARY_COHORTS,
    "primary_plus_sensitivity": ALL_BULK_COHORTS,
}
C_GRID = np.logspace(-4, 2, 7)

PROJECT = Path(__file__).resolve().parents[2]
TRAJECTORY_CSV = PROJECT / "04_manuscript/plot_data_local/powered/out_04_trajectory/trajectory_pseudotime.csv"
DONOR_MAP_CSV = PROJECT / "04_manuscript/revision_20260722/donor_metadata/outputs/sample_to_donor.csv"
FEATURE_DIR = PROJECT / "04_manuscript/plot_data_local/out_05_bulk/features"
META_DIR = PROJECT / "server_archive/tier1/out_05_bulk/features"
DATA_AUDIT_OUT = PROJECT / "04_manuscript/revision_20260722/data_audit/outputs"
GSE58095_SOFT = PROJECT / "02_analysis_modules/ml_validation/geo/GSE58095_family.soft.gz"
EXTERNAL_DIR = PROJECT / "04_manuscript/revision_20260722/figure3/external_validation/outputs"
OUT = PROJECT / "04_manuscript/revision_20260722/figure3/outputs"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def adjust_pvalues(frame: pd.DataFrame, p_col: str = "p_value") -> pd.DataFrame:
    out = frame.copy()
    valid = out[p_col].notna()
    out["q_bh"] = np.nan
    out["p_holm"] = np.nan
    if valid.any():
        out.loc[valid, "q_bh"] = multipletests(out.loc[valid, p_col], method="fdr_bh")[1]
        out.loc[valid, "p_holm"] = multipletests(out.loc[valid, p_col], method="holm")[1]
    return out


def bootstrap_interval(values: np.ndarray, statistic, n_boot: int, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    draws = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        draws[i] = statistic(sample)
    return tuple(np.quantile(draws[np.isfinite(draws)], [0.025, 0.975]))


def load_true_donor_cells() -> pd.DataFrame:
    trajectory = pd.read_csv(TRAJECTORY_CSV)
    donor_map = pd.read_csv(DONOR_MAP_CSV)
    keep_cols = ["sample", "cohort", "condition", "donor_id", "analysis_keep", "exclusion_reason"]
    cells = trajectory.merge(
        donor_map[keep_cols], on=["sample", "cohort"], how="left", validate="many_to_one", indicator=True
    )
    if not cells["_merge"].eq("both").all():
        missing = cells.loc[~cells["_merge"].eq("both"), ["sample", "cohort"]].drop_duplicates()
        raise RuntimeError(f"Trajectory samples missing from donor map:\n{missing.to_string(index=False)}")
    cells = cells.loc[cells["analysis_keep"].eq(True)].copy()
    if cells["condition"].isna().any() or cells["donor_id"].isna().any():
        raise RuntimeError("Retained cells contain missing condition or donor_id")
    return cells.drop(columns="_merge")


def donor_activity(cells: pd.DataFrame, n_boot: int, rng: np.random.Generator) -> None:
    myo = cells.loc[cells["fib_subtype"].eq("Myofibroblast")].copy()
    rows = []
    for tf, column in REGULON_COLUMNS.items():
        agg = (
            myo.groupby(["donor_id", "cohort", "condition"], observed=True)
            .agg(activity=(column, "mean"), n_myofibroblasts=(column, "size"))
            .reset_index()
        )
        agg["TF"] = tf
        rows.append(agg)
    source = pd.concat(rows, ignore_index=True)
    source.to_csv(OUT / "Figure3a_donor_activity_source_data.csv", index=False)

    results = []
    loco_rows = []
    for tf in TF_PANEL:
        d = source.loc[source["TF"].eq(tf)].copy()
        d["is_ssc"] = d["condition"].eq("SSc").astype(int)
        fit = smf.ols("activity ~ is_ssc + C(cohort)", data=d).fit(cov_type="HC3")
        beta = float(fit.params["is_ssc"])
        ci_lo, ci_hi = map(float, fit.conf_int().loc["is_ssc"])
        resid_sd = float(np.sqrt(np.sum(fit.resid**2) / fit.df_resid))
        ssc = d.loc[d["is_ssc"].eq(1), "activity"].to_numpy()
        hc = d.loc[d["is_ssc"].eq(0), "activity"].to_numpy()
        pooled_sd = np.sqrt(((len(ssc) - 1) * np.var(ssc, ddof=1) + (len(hc) - 1) * np.var(hc, ddof=1)) / (len(d) - 2))
        correction = 1 - 3 / (4 * (len(d) - 2) - 1)
        hedges_g = correction * (np.mean(ssc) - np.mean(hc)) / pooled_sd
        median_diff = float(np.median(ssc) - np.median(hc))
        boot_diff = np.empty(n_boot)
        for i in range(n_boot):
            boot_diff[i] = np.median(rng.choice(ssc, len(ssc), replace=True)) - np.median(
                rng.choice(hc, len(hc), replace=True)
            )
        results.append(
            {
                "TF": tf,
                "n_donors": len(d),
                "n_SSc": len(ssc),
                "n_HC": len(hc),
                "median_SSc": np.median(ssc),
                "median_HC": np.median(hc),
                "median_difference": median_diff,
                "median_difference_ci_low": np.quantile(boot_diff, 0.025),
                "median_difference_ci_high": np.quantile(boot_diff, 0.975),
                "adjusted_mean_difference": beta,
                "adjusted_ci_low": ci_lo,
                "adjusted_ci_high": ci_hi,
                "standardized_adjusted_difference": beta / resid_sd,
                "hedges_g_unadjusted": hedges_g,
                "p_value": float(fit.pvalues["is_ssc"]),
                "r_squared": float(fit.rsquared),
                "model": "OLS activity ~ condition + cohort; HC3 robust SE",
            }
        )
        for omitted in sorted(d["cohort"].unique()):
            sub = d.loc[~d["cohort"].eq(omitted)].copy()
            sub["is_ssc"] = sub["condition"].eq("SSc").astype(int)
            if sub["is_ssc"].nunique() < 2:
                continue
            loco = smf.ols("activity ~ is_ssc + C(cohort)", data=sub).fit(cov_type="HC3")
            lo, hi = map(float, loco.conf_int().loc["is_ssc"])
            loco_rows.append(
                {
                    "TF": tf,
                    "omitted_cohort": omitted,
                    "n_donors": len(sub),
                    "adjusted_mean_difference": float(loco.params["is_ssc"]),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_value": float(loco.pvalues["is_ssc"]),
                }
            )
    adjust_pvalues(pd.DataFrame(results)).to_csv(OUT / "Figure3a_donor_activity_results.csv", index=False)
    pd.DataFrame(loco_rows).to_csv(OUT / "Figure3a_leave_one_cohort_out.csv", index=False)

    cohort_summary = (
        source.groupby(["TF", "cohort", "condition"], observed=True)
        .agg(n_donors=("donor_id", "nunique"), median_activity=("activity", "median"), mean_activity=("activity", "mean"))
        .reset_index()
    )
    cohort_summary.to_csv(OUT / "Figure3a_cohort_summary.csv", index=False)


def donor_trajectory(cells: pd.DataFrame, n_boot: int, rng: np.random.Generator) -> None:
    cell_counts = cells.groupby(["donor_id", "cohort", "condition"], observed=True).size().rename("n_cells").reset_index()
    quartiles = cell_counts["n_cells"].quantile([0.25, 0.50, 0.75]).round().astype(int).tolist()
    thresholds = sorted(set([20, 30] + quartiles))
    audit = pd.DataFrame(
        {
            "threshold": thresholds,
            "source": [
                "inherited minimum" if x == 30 else ("descriptive lower sensitivity" if x == 20 else "outcome-blind donor cell-count quartile")
                for x in thresholds
            ],
            "n_donors": [int((cell_counts["n_cells"] >= x).sum()) for x in thresholds],
        }
    )
    audit.to_csv(OUT / "Figure3b_cell_count_threshold_audit.csv", index=False)

    rho_rows = []
    for donor_id, d in cells.groupby("donor_id", observed=True):
        n_cells = len(d)
        meta = d.iloc[0]
        for tf, column in REGULON_COLUMNS.items():
            ok = np.isfinite(d[column]) & np.isfinite(d["myo_fate"])
            rho = spearmanr(d.loc[ok, column], d.loc[ok, "myo_fate"]).statistic if ok.sum() >= 3 else np.nan
            rho_rows.append(
                {
                    "donor_id": donor_id,
                    "cohort": meta["cohort"],
                    "condition": meta["condition"],
                    "TF": tf,
                    "rho": rho,
                    "n_cells": int(ok.sum()),
                }
            )
    source = pd.DataFrame(rho_rows)
    source.to_csv(OUT / "Figure3b_donor_trajectory_source_data.csv", index=False)

    primary = source.loc[source["n_cells"].ge(30) & source["rho"].notna()].copy()
    results = []
    for tf in TF_PANEL:
        d = primary.loc[primary["TF"].eq(tf)]
        vals = d["rho"].to_numpy()
        lo, hi = bootstrap_interval(vals, np.median, n_boot, rng)
        w = wilcoxon(vals, alternative="two-sided", zero_method="wilcox")
        positive = int((vals > 0).sum())
        sign_ci = binomtest(positive, len(vals), 0.5).proportion_ci(confidence_level=0.95, method="wilson")
        results.append(
            {
                "TF": tf,
                "n_donors": len(vals),
                "n_SSc": int(d["condition"].eq("SSc").sum()),
                "n_HC": int(d["condition"].eq("HC").sum()),
                "minimum_cells": 30,
                "median_rho": np.median(vals),
                "median_rho_ci_low": lo,
                "median_rho_ci_high": hi,
                "positive_fraction": positive / len(vals),
                "positive_fraction_ci_low": sign_ci.low,
                "positive_fraction_ci_high": sign_ci.high,
                "p_value": float(w.pvalue),
                "sign_test_p": float(binomtest(positive, len(vals), 0.5).pvalue),
            }
        )
    adjust_pvalues(pd.DataFrame(results)).to_csv(OUT / "Figure3b_trajectory_results.csv", index=False)

    sensitivity = []
    for threshold in thresholds:
        for tf in TF_PANEL:
            d = source.loc[source["TF"].eq(tf) & source["n_cells"].ge(threshold) & source["rho"].notna()]
            vals = d["rho"].to_numpy()
            if len(vals) == 0:
                continue
            sensitivity.append(
                {
                    "minimum_cells": threshold,
                    "TF": tf,
                    "n_donors": len(vals),
                    "median_rho": np.median(vals),
                    "positive_fraction": np.mean(vals > 0),
                    "wilcoxon_p": wilcoxon(vals, alternative="two-sided", zero_method="wilcox").pvalue,
                }
            )
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity["q_bh_within_threshold"] = sensitivity.groupby("minimum_cells")["wilcoxon_p"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1]
    )
    sensitivity.to_csv(OUT / "Figure3b_cell_count_sensitivity.csv", index=False)

    cohort_summary = (
        primary.groupby(["TF", "cohort", "condition"], observed=True)
        .agg(n_donors=("donor_id", "nunique"), median_rho=("rho", "median"), positive_fraction=("rho", lambda x: np.mean(x > 0)))
        .reset_index()
    )
    cohort_summary.to_csv(OUT / "Figure3b_cohort_summary.csv", index=False)


def labels_from_meta(meta: pd.DataFrame) -> pd.Series:
    label_col = next(c for c in ("grp", "SSc", "status", "group", "condition", "disease") if c in meta.columns)

    def parse(value) -> float:
        text = str(value).strip().upper()
        if "SSC" in text or text in {"1", "TRUE", "CASE"}:
            return 1.0
        if any(token in text for token in ("HC", "CTRL", "CONTROL", "NORMAL", "HEALTHY")) or text in {"0", "FALSE"}:
            return 0.0
        return np.nan

    return meta[label_col].map(parse)


def parse_soft_characteristics(values: list[str]) -> dict[str, str]:
    parsed = {}
    for value in values:
        if ":" not in value:
            continue
        key, item = value.split(":", 1)
        parsed[key.strip().lower()] = item.strip()
    return parsed


def gse58095_subject_map() -> pd.DataFrame:
    gse = GEOparse.get_GEO(filepath=str(GSE58095_SOFT), silent=True)
    rows = []
    for geo_accession, sample in gse.gsms.items():
        fields = parse_soft_characteristics(
            sample.metadata.get("characteristics_ch1", [])
        )
        raw_subject = fields.get("individual") or fields.get("patient") or ""
        subject = re.sub(
            r"^PATIENT0*",
            "PATIENT",
            raw_subject.strip().upper(),
        )
        timepoint = fields.get("time point", "").strip().lower()
        rows.append(
            {
                "cohort": "GSE58095",
                "sample": geo_accession,
                "raw_subject_id": raw_subject,
                "subject_id": f"GSE58095:{subject}",
                "timepoint": timepoint or "deposited_initial",
                "time_rank": 1 if timepoint == "late" else 0,
            }
        )
    mapping = pd.DataFrame(rows)
    if mapping["subject_id"].str.endswith(":").any():
        raise RuntimeError("GSE58095 SOFT contains an unresolved participant ID")
    minimum_rank = mapping.groupby("subject_id")["time_rank"].transform("min")
    mapping["analysis_keep"] = mapping["time_rank"].eq(minimum_rank)
    mapping["exclusion_reason"] = np.where(
        mapping["analysis_keep"],
        "",
        "later longitudinal biopsy; earliest deposited time retained",
    )
    mapping["aggregation_rule"] = (
        "mean activity across retained same-time biopsies within participant"
    )
    mapping.to_csv(DATA_AUDIT_OUT / "GSE58095_sample_to_subject.csv", index=False)
    return mapping


def legacy_subject_map(
    cohort: str,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    if cohort == "GSE130955":
        mapping = pd.read_csv(DATA_AUDIT_OUT / "GSE130955_sample_to_subject.csv")
        mapping = mapping.rename(
            columns={"sample": "sample", "subject_id": "raw_subject_id"}
        )
        mapping["cohort"] = cohort
        mapping["subject_id"] = cohort + ":" + mapping["raw_subject_id"].astype(str)
        mapping["analysis_keep"] = True
        mapping["exclusion_reason"] = ""
        mapping["aggregation_rule"] = (
            "mean activity across all deposited biopsies within participant"
        )
        return mapping[
            [
                "cohort",
                "sample",
                "raw_subject_id",
                "subject_id",
                "analysis_keep",
                "exclusion_reason",
                "aggregation_rule",
            ]
        ]
    if cohort == "GSE58095":
        return gse58095_subject_map()
    if cohort == "GSE249550":
        if "subject" not in meta or meta["subject"].isna().any():
            raise RuntimeError("GSE249550 lacks complete subject IDs")
        raw_subject = meta["subject"].astype(str)
        if raw_subject.duplicated().any():
            raise RuntimeError("GSE249550 baseline participant IDs are not unique")
        return pd.DataFrame(
            {
                "cohort": cohort,
                "sample": meta.index.astype(str),
                "raw_subject_id": raw_subject.to_numpy(),
                "subject_id": (cohort + ":" + raw_subject).to_numpy(),
                "analysis_keep": True,
                "exclusion_reason": "",
                "aggregation_rule": "one baseline sample per participant",
            }
        )
    raise KeyError(cohort)


def load_legacy_cohort(
    cohort: str,
) -> tuple[pd.DataFrame, pd.Series, dict[str, object], pd.DataFrame]:
    activity = pd.read_csv(
        FEATURE_DIR / f"{cohort}_regulon_activity.csv", index_col=0
    )
    activity.index = activity.index.astype(str)
    meta = pd.read_csv(META_DIR / f"{cohort}_meta.csv", index_col=0)
    meta.index = meta.index.astype(str)
    raw_labels = labels_from_meta(meta).dropna().astype(int)
    mapping = legacy_subject_map(cohort, meta)
    mapping["sample"] = mapping["sample"].astype(str)
    available_samples = activity.index.intersection(raw_labels.index)
    mapping = mapping.loc[mapping["sample"].isin(available_samples)]
    if mapping["sample"].duplicated().any():
        raise RuntimeError(f"{cohort}: sample-to-subject mapping is not unique")
    retained = mapping.loc[mapping["analysis_keep"]].set_index("sample")
    sample_ids = activity.index.intersection(retained.index).intersection(raw_labels.index)
    if len(sample_ids) != len(retained):
        missing = sorted(set(retained.index) - set(sample_ids))
        raise RuntimeError(f"{cohort}: mapped retained samples missing from features: {missing}")
    sample_activity = activity.loc[sample_ids].copy()
    sample_activity["subject_id"] = retained.loc[sample_ids, "subject_id"]
    subject_activity = sample_activity.groupby("subject_id", sort=True).mean()
    sample_labels = raw_labels.loc[sample_ids].rename("label").to_frame()
    sample_labels["subject_id"] = retained.loc[sample_ids, "subject_id"]
    conflicts = sample_labels.groupby("subject_id")["label"].nunique()
    if conflicts.gt(1).any():
        raise RuntimeError(f"{cohort}: a participant has conflicting disease labels")
    subject_labels = sample_labels.groupby("subject_id")["label"].first()
    source_counts = sample_labels.groupby("subject_id").size()
    if not subject_activity.index.equals(subject_labels.index):
        raise RuntimeError(f"{cohort}: participant activity/label mismatch")
    prevalence = {
        "cohort": cohort,
        "analysis_role": "primary",
        "n_source_samples": len(sample_ids),
        "n_subjects": len(subject_labels),
        "n_SSc": int(subject_labels.sum()),
        "n_HC": int((1 - subject_labels).sum()),
        "prevalence_SSc": float(subject_labels.mean()),
        "participants_with_multiple_retained_samples": int(source_counts.gt(1).sum()),
        "subject_ids_available": True,
        "mapping_source": (
            "NCBI GEO family SOFT"
            if cohort in {"GSE130955", "GSE58095"}
            else "deposited baseline metadata"
        ),
    }
    return subject_activity, subject_labels, prevalence, mapping


def load_external_cohort(
    cohort: str,
) -> tuple[pd.DataFrame, pd.Series, dict[str, object], pd.DataFrame]:
    activity = pd.read_csv(
        EXTERNAL_DIR / f"{cohort}_subject_regulon_activity.csv", index_col=0
    )
    meta = pd.read_csv(
        EXTERNAL_DIR / f"{cohort}_subject_metadata.csv", index_col=0
    )
    common = activity.index.intersection(meta.index)
    if len(common) != len(meta):
        raise RuntimeError(f"{cohort}: external participant activity/metadata mismatch")
    activity = activity.loc[common].copy()
    labels = meta.loc[common, "condition"].map({"HC": 0, "SSc": 1})
    if labels.isna().any():
        raise RuntimeError(f"{cohort}: unexpected external condition")
    prefixed = pd.Index([f"{cohort}:{value}" for value in common], name="subject_id")
    activity.index = prefixed
    labels.index = prefixed
    mapping = pd.read_csv(EXTERNAL_DIR / f"{cohort}_sample_to_subject.csv")
    mapping.insert(0, "cohort", cohort)
    mapping["raw_subject_id"] = mapping["subject_id"]
    mapping["subject_id"] = cohort + ":" + mapping["subject_id"].astype(str)
    mapping = mapping.rename(columns={"geo_accession": "sample"})
    mapping["aggregation_rule"] = (
        "mean normalized expression across retained same-time samples before ULM"
    )
    retained_counts = meta.loc[common, "n_source_samples"]
    prevalence = {
        "cohort": cohort,
        "analysis_role": (
            "primary" if cohort in PRIMARY_COHORTS else "sensitivity"
        ),
        "n_source_samples": int(retained_counts.sum()),
        "n_subjects": len(labels),
        "n_SSc": int(labels.sum()),
        "n_HC": int((1 - labels).sum()),
        "prevalence_SSc": float(labels.mean()),
        "participants_with_multiple_retained_samples": int(
            retained_counts.gt(1).sum()
        ),
        "subject_ids_available": True,
        "mapping_source": "audited NCBI GEO family SOFT",
    }
    return activity, labels.astype(int), prevalence, mapping


def load_bulk() -> tuple[
    dict[str, pd.DataFrame],
    dict[str, pd.Series],
    pd.DataFrame,
    pd.DataFrame,
]:
    activities: dict[str, pd.DataFrame] = {}
    labels: dict[str, pd.Series] = {}
    prevalence_rows = []
    mapping_rows = []
    for cohort in ALL_BULK_COHORTS:
        if cohort in {"GSE130955", "GSE58095", "GSE249550"}:
            activity, y, prevalence, mapping = load_legacy_cohort(cohort)
        else:
            activity, y, prevalence, mapping = load_external_cohort(cohort)
        if "HES1" not in activity:
            raise RuntimeError(f"{cohort}: HES1 absent from regulon features")
        activities[cohort] = activity
        labels[cohort] = y
        prevalence_rows.append(prevalence)
        mapping_rows.append(mapping)
    return (
        activities,
        labels,
        pd.DataFrame(prevalence_rows),
        pd.concat(mapping_rows, ignore_index=True, sort=False),
    )


def train_predict(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    c_value: float,
) -> tuple[np.ndarray, dict[str, object]]:
    imputer = SimpleImputer(strategy="median").fit(train_x)
    train_i = imputer.transform(train_x)
    test_i = imputer.transform(test_x)
    scaler = StandardScaler().fit(train_i)
    model = LogisticRegression(
        C=float(c_value), solver="liblinear", max_iter=5000, random_state=SEED
    ).fit(scaler.transform(train_i), train_y)
    prediction = model.predict_proba(scaler.transform(test_i))[:, 1]
    return prediction, {"imputer": imputer, "scaler": scaler, "model": model}


def select_c(
    train_cohorts: list[str],
    features: list[str],
    activities: dict[str, pd.DataFrame],
    labels: dict[str, pd.Series],
) -> tuple[float, list[dict[str, object]]]:
    rows = []
    for c_value in C_GRID:
        scores = []
        for inner_held in train_cohorts:
            inner_train = [c for c in train_cohorts if c != inner_held]
            x_train = np.vstack([activities[c][features].to_numpy() for c in inner_train])
            y_train = np.concatenate([labels[c].to_numpy() for c in inner_train])
            pred, _ = train_predict(x_train, y_train, activities[inner_held][features].to_numpy(), c_value)
            scores.append(roc_auc_score(labels[inner_held], pred))
        rows.append({"C": float(c_value), "mean_inner_cohort_AUROC": float(np.mean(scores)), "inner_scores": scores})
    best = sorted(rows, key=lambda row: (-row["mean_inner_cohort_AUROC"], row["C"]))[0]
    return float(best["C"]), rows


def stratified_bootstrap_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    return np.concatenate([rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)])


def calibration_statistics(y: np.ndarray, prediction: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(prediction, 1e-6, 1 - 1e-6)
    design = sm.add_constant(logit(clipped))
    try:
        fit = sm.GLM(y, design, family=sm.families.Binomial()).fit()
        return float(fit.params[0]), float(fit.params[1])
    except Exception:
        return np.nan, np.nan


def bootstrap_prediction_metrics(
    y: np.ndarray, prediction: np.ndarray, n_boot: int, rng: np.random.Generator
) -> dict[str, float]:
    values = np.empty((n_boot, 3), dtype=float)
    for i in range(n_boot):
        idx = stratified_bootstrap_indices(y, rng)
        values[i, 0] = roc_auc_score(y[idx], prediction[idx])
        values[i, 1] = average_precision_score(y[idx], prediction[idx])
        values[i, 2] = brier_score_loss(y[idx], prediction[idx])
    quantiles = np.quantile(values, [0.025, 0.975], axis=0)
    return {
        "auroc_ci_low": quantiles[0, 0],
        "auroc_ci_high": quantiles[1, 0],
        "average_precision_ci_low": quantiles[0, 1],
        "average_precision_ci_high": quantiles[1, 1],
        "brier_ci_low": quantiles[0, 2],
        "brier_ci_high": quantiles[1, 2],
    }


def paired_ablation(
    predictions: pd.DataFrame,
    n_boot: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    contrasts = (("full_panel", "hes1_only"), ("full_panel", "remove_hes1"))
    for analysis_set, cohort_order in ANALYSIS_SETS.items():
        for cohort in cohort_order:
            cohort_predictions = predictions.loc[
                predictions["analysis_set"].eq(analysis_set)
                & predictions["held_out_cohort"].eq(cohort)
            ]
            wide = cohort_predictions.pivot(
                index="sample", columns="model", values="prediction"
            )
            y = (
                cohort_predictions.drop_duplicates("sample")
                .set_index("sample")
                .loc[wide.index, "y"]
                .to_numpy()
            )
            for model_a, model_b in contrasts:
                pa = wide[model_a].to_numpy()
                pb = wide[model_b].to_numpy()
                for metric_name, metric in (
                    ("AUROC", roc_auc_score),
                    ("average_precision", average_precision_score),
                ):
                    observed = float(metric(y, pa) - metric(y, pb))
                    draws = np.empty(n_boot)
                    for i in range(n_boot):
                        idx = stratified_bootstrap_indices(y, rng)
                        draws[i] = metric(y[idx], pa[idx]) - metric(y[idx], pb[idx])
                    p_value = (
                        2 * min(np.sum(draws <= 0), np.sum(draws >= 0)) + 1
                    ) / (n_boot + 1)
                    rows.append(
                        {
                            "analysis_set": analysis_set,
                            "held_out_cohort": cohort,
                            "metric": metric_name,
                            "model_a": model_a,
                            "model_b": model_b,
                            "difference_a_minus_b": observed,
                            "ci_low": np.quantile(draws, 0.025),
                            "ci_high": np.quantile(draws, 0.975),
                            "p_value": min(float(p_value), 1.0),
                        }
                    )
    return adjust_pvalues(pd.DataFrame(rows))


def bulk_ablation(n_boot: int, rng: np.random.Generator) -> None:
    activities, labels, prevalence, subject_mapping = load_bulk()
    prevalence.to_csv(OUT / "Figure3de_bulk_cohort_prevalence.csv", index=False)
    subject_mapping.to_csv(
        OUT / "Figure3de_bulk_sample_to_subject.csv", index=False
    )
    prediction_rows = []
    metric_rows = []
    tuning_rows = []
    calibration_rows = []
    feature_set_rows = []

    for analysis_set, cohort_order in ANALYSIS_SETS.items():
        shared = sorted(
            set.intersection(*(set(activities[c].columns) for c in cohort_order))
        )
        if "HES1" not in shared:
            raise RuntimeError(
                f"{analysis_set}: HES1 absent from shared regulon feature panel"
            )
        feature_sets = {
            "hes1_only": ["HES1"],
            "full_panel": shared,
            "remove_hes1": [feature for feature in shared if feature != "HES1"],
        }
        for model_name, features in feature_sets.items():
            feature_set_rows.append(
                {
                    "analysis_set": analysis_set,
                    "model": model_name,
                    "n_features": len(features),
                    "features": ";".join(features),
                }
            )

        for held_out in cohort_order:
            train_cohorts = [cohort for cohort in cohort_order if cohort != held_out]
            for model_name in MODEL_ORDER:
                features = feature_sets[model_name]
                selected_c, tuning = select_c(
                    train_cohorts, features, activities, labels
                )
                for row in tuning:
                    tuning_rows.append(
                        {
                            "analysis_set": analysis_set,
                            "held_out_cohort": held_out,
                            "model": model_name,
                            "n_features": len(features),
                            "C": row["C"],
                            "mean_inner_cohort_AUROC": row[
                                "mean_inner_cohort_AUROC"
                            ],
                            "inner_AUROCs": ";".join(
                                f"{score:.8g}" for score in row["inner_scores"]
                            ),
                        }
                    )
                x_train = np.vstack(
                    [activities[c][features].to_numpy() for c in train_cohorts]
                )
                y_train = np.concatenate(
                    [labels[c].to_numpy() for c in train_cohorts]
                )
                x_test = activities[held_out][features].to_numpy()
                y_test = labels[held_out].to_numpy()
                prediction, _ = train_predict(
                    x_train, y_train, x_test, selected_c
                )
                boot = bootstrap_prediction_metrics(
                    y_test, prediction, n_boot, rng
                )
                calibration_intercept, calibration_slope = (
                    calibration_statistics(y_test, prediction)
                )
                metric_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "held_out_cohort": held_out,
                        "model": model_name,
                        "n_features": len(features),
                        "selected_C": selected_c,
                        "n_test": len(y_test),
                        "n_SSc": int(y_test.sum()),
                        "n_HC": int((1 - y_test).sum()),
                        "prevalence_SSc": float(y_test.mean()),
                        "AUROC": roc_auc_score(y_test, prediction),
                        "average_precision": average_precision_score(
                            y_test, prediction
                        ),
                        "brier_score": brier_score_loss(y_test, prediction),
                        "calibration_intercept": calibration_intercept,
                        "calibration_slope": calibration_slope,
                        **boot,
                    }
                )
                sample_names = activities[held_out].index.to_numpy()
                for sample, y_value, probability in zip(
                    sample_names, y_test, prediction
                ):
                    prediction_rows.append(
                        {
                            "sample": sample,
                            "analysis_set": analysis_set,
                            "held_out_cohort": held_out,
                            "model": model_name,
                            "y": int(y_value),
                            "prediction": float(probability),
                            "selected_C": selected_c,
                            "n_features": len(features),
                        }
                    )
                ranked = pd.Series(prediction).rank(method="first")
                bins = pd.qcut(
                    ranked, q=min(5, len(prediction)), labels=False
                )
                cal = pd.DataFrame(
                    {"y": y_test, "prediction": prediction, "bin": bins}
                )
                cal = cal.groupby("bin", observed=True).agg(
                    n=("y", "size"),
                    mean_prediction=("prediction", "mean"),
                    observed_fraction=("y", "mean"),
                )
                for bin_id, row in cal.iterrows():
                    calibration_rows.append(
                        {
                            "analysis_set": analysis_set,
                            "held_out_cohort": held_out,
                            "model": model_name,
                            "bin": int(bin_id),
                            "n": int(row["n"]),
                            "mean_prediction": float(row["mean_prediction"]),
                            "observed_fraction": float(
                                row["observed_fraction"]
                            ),
                        }
                    )

    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(OUT / "Figure3de_bulk_predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(OUT / "Figure3d_bulk_ablation_metrics.csv", index=False)
    pd.DataFrame(tuning_rows).to_csv(OUT / "Figure3d_nested_tuning_audit.csv", index=False)
    pd.DataFrame(calibration_rows).to_csv(OUT / "Figure3e_calibration_bins.csv", index=False)
    paired_ablation(predictions, n_boot, rng).to_csv(
        OUT / "Figure3d_paired_ablation_tests.csv", index=False
    )
    pd.DataFrame(feature_set_rows).to_csv(
        OUT / "Figure3d_feature_sets.csv", index=False
    )


def write_manifest(n_boot: int) -> None:
    input_paths = [TRAJECTORY_CSV, DONOR_MAP_CSV]
    for cohort in ("GSE130955", "GSE58095", "GSE249550"):
        input_paths.extend(
            [FEATURE_DIR / f"{cohort}_regulon_activity.csv", META_DIR / f"{cohort}_meta.csv"]
        )
    input_paths.extend(
        [
            DATA_AUDIT_OUT / "GSE130955_sample_to_subject.csv",
            DATA_AUDIT_OUT / "GSE58095_sample_to_subject.csv",
            GSE58095_SOFT,
            EXTERNAL_DIR / "external_feature_manifest.json",
        ]
    )
    for cohort in ALL_BULK_COHORTS:
        if cohort in {"GSE130955", "GSE58095", "GSE249550"}:
            continue
        input_paths.extend(
            [
                EXTERNAL_DIR / f"{cohort}_subject_regulon_activity.csv",
                EXTERNAL_DIR / f"{cohort}_subject_metadata.csv",
                EXTERNAL_DIR / f"{cohort}_processing_qc.json",
            ]
        )
    manifest = {
        "seed": SEED,
        "bootstrap_replicates": n_boot,
        "trajectory_primary_minimum_cells": 30,
        "bulk_outer_validation": "leave one complete cohort out",
        "bulk_inner_tuning": (
            "leave one complete training cohort out; cohort AUROCs averaged equally"
        ),
        "bulk_analysis_sets": {
            name: list(cohorts) for name, cohorts in ANALYSIS_SETS.items()
        },
        "bulk_independent_unit": "participant",
        "preprocessing": "median imputation and standardization fitted only in each training split",
        "models": {
            "hes1_only": "L2 logistic regression using HES1 activity only",
            "full_panel": "L2 logistic regression using all shared regulon activities",
            "remove_hes1": "same full panel after removing HES1 activity",
        },
        "C_grid": C_GRID.tolist(),
        "input_sha256": {str(path.relative_to(PROJECT)): sha256(path) for path in input_paths},
        "known_boundary": (
            "The deposited local trajectory output contains fate probabilities but not the expression matrix or graph. "
            "Excluding HES1 and its regulon genes from trajectory construction therefore requires a separate CellRank rerun."
        ),
    }
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=10_000)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cells = load_true_donor_cells()
    donor_activity(cells, args.bootstrap, rng)
    donor_trajectory(cells, args.bootstrap, rng)
    bulk_ablation(args.bootstrap, rng)
    write_manifest(args.bootstrap)
    print(f"Figure 3 validation outputs written to {OUT}")


if __name__ == "__main__":
    main()
