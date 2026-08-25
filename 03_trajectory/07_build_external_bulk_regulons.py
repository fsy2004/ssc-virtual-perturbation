#!/usr/bin/env python
"""Build donor-level CollecTRI/ULM features for external SSc skin cohorts.

This extends the existing bulk workflow with outcome-blind cohort screening and
one row per participant. GEOparse performs deposited matrix/platform parsing;
decoupler's OmniPath CollecTRI resource and ULM implementation provide the TF
activity model. No custom rank score is used.

Adapted from the repository's ``04_bulk_regulon_ml/05_bulk_regulon_severity.py``
and reusable modules 008/009 (GEO annotation) and 076 (decoupler activity).
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
import re
from pathlib import Path

import anndata as ad
import decoupler as dc
import GEOparse
import mygene
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
AUDIT = PROJECT / "04_manuscript/revision_20260722/data_audit"
SOFT = AUDIT / "geo_soft_cache"
SAMPLE_AUDIT = AUDIT / "outputs/candidate_sample_metadata.csv"
OUT = PROJECT / "04_manuscript/revision_20260722/figure3/external_validation/outputs"
GSE231692_COUNTS = PROJECT / ".codex_tmp/external_bulk/GSE231692_count.tsv.gz"

COHORTS = {
    "GSE181549": {
        "platform": "GPL13497",
        "map_column": "GENE_SYMBOL",
        "role": "primary_large_adult_whole_skin",
    },
    "GSE9285": {
        "platform": "GPL5981",
        "map_column": "GENE",
        "role": "primary_adult_whole_skin",
    },
    "GSE32413": {
        "platform": "GPL4133",
        "map_column": "GENE_SYMBOL",
        "role": "primary_adult_whole_skin_baseline",
    },
    "GSE76807": {
        "platform": "GPL6480",
        "map_column": "GENE_SYMBOL",
        "role": "primary_adult_limited_SSc_whole_skin",
    },
    "GSE76885": {
        "platform": "GPL6480",
        "map_column": "GENE_SYMBOL",
        "role": "sensitivity_adult_longitudinal_baseline",
    },
    "GSE95065": {
        "platform": "GPL23080",
        "map_column": "ENTREZ_GENE_ID",
        "role": "sensitivity_sparse_publication_and_batch_metadata",
    },
    "GSE125362": {
        "platform": "GPL6480",
        "map_column": "GENE_SYMBOL",
        "role": "sensitivity_small_paired_platform_cohort",
    },
    "GSE231692": {
        "platform": "RNA-seq raw counts",
        "map_column": "Ensembl_via_org.Hs.eg.db",
        "role": "primary_adult_whole_skin_baseline",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_symbol(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "null", "---"}:
        return None
    text = re.split(r"\s*(?:///|//|;|\|)\s*", text, maxsplit=1)[0].strip()
    return text or None


def entrez_map(values: pd.Series) -> pd.Series:
    cache = OUT / "GPL23080_entrez_to_symbol.csv"
    if cache.exists():
        table = pd.read_csv(cache, dtype=str)
    else:
        ids = (
            pd.to_numeric(values, errors="coerce")
            .dropna()
            .astype("int64")
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        query = mygene.MyGeneInfo().querymany(
            ids,
            scopes="entrezgene",
            fields="symbol",
            species="human",
            as_dataframe=False,
            verbose=False,
        )
        rows = [
            {"entrez_id": str(item["query"]), "symbol": item.get("symbol", "")}
            for item in query
            if item.get("symbol") and not item.get("notfound", False)
        ]
        table = pd.DataFrame(rows).drop_duplicates("entrez_id")
        table.to_csv(cache, index=False)
    mapping = dict(zip(table["entrez_id"], table["symbol"]))
    numeric = pd.to_numeric(values, errors="coerce")
    return numeric.map(lambda value: mapping.get(str(int(value))) if pd.notna(value) else None)


def coalesce(frame: pd.DataFrame, pattern: str) -> pd.Series:
    columns = [column for column in frame.columns if re.search(pattern, column, re.I)]
    if not columns:
        return pd.Series(index=frame.index, dtype=object)
    return frame[columns].bfill(axis=1).iloc[:, 0]


def audited_metadata(accession: str) -> pd.DataFrame:
    rows = pd.read_csv(SAMPLE_AUDIT, dtype=str).fillna("")
    rows = rows.loc[rows["accession"].eq(accession)].copy()
    return rows.set_index("geo_accession", drop=False)


def subject_map(accession: str) -> pd.DataFrame:
    rows = audited_metadata(accession)
    out = pd.DataFrame(index=rows.index)
    out["geo_accession"] = rows["geo_accession"]
    out["title"] = rows["title"]
    out["source_name"] = rows["source_name"]
    out["condition"] = ""
    out["subject_id"] = ""
    out["time_rank"] = 0
    out["analysis_keep"] = True
    out["exclusion_reason"] = ""

    if accession == "GSE181549":
        is_ssc = rows["characteristics"].str.contains(
            "condition: systemic sclerosis patient", case=False, regex=False
        )
        is_hc = rows["characteristics"].str.contains(
            "condition: Control", case=False, regex=False
        )
        out.loc[is_ssc, "condition"] = "SSc"
        out.loc[is_hc, "condition"] = "HC"
        out.loc[is_ssc, "subject_id"] = rows.loc[is_ssc, "subject_id"]
        out.loc[is_hc, "subject_id"] = "HC_" + rows.loc[is_hc, "title"]
        sample_number = pd.to_numeric(
            rows["characteristics"].str.extract(
                r"sample number within subject:\s*([^|]+)", expand=False
            ),
            errors="coerce",
        )
        out.loc[is_ssc & sample_number.ne(1), "analysis_keep"] = False
        out.loc[is_ssc & sample_number.ne(1), "exclusion_reason"] = (
            "follow-up biopsy; first deposited biopsy is primary"
        )
    elif accession == "GSE9285":
        title = rows["title"]
        is_ssc = title.str.match(r"^(?:dSSc|lSSc)\d+", case=False)
        is_hc = title.str.match(r"^Nor\d+", case=False)
        excluded_disease = rows["characteristics"].str.contains(
            "morphea|eosinophilic fasciitis", case=False, regex=True
        )
        out.loc[is_ssc & ~excluded_disease, "condition"] = "SSc"
        out.loc[is_hc, "condition"] = "HC"
        out["subject_id"] = title.str.extract(
            r"^((?:dSSc|lSSc|Nor)\d+)", flags=re.I, expand=False
        ).fillna("")
        invalid = (~is_ssc & ~is_hc) | excluded_disease
        out.loc[invalid, "analysis_keep"] = False
        out.loc[invalid, "exclusion_reason"] = "morphea, eosinophilic fasciitis, or non-SSc comparator"
    elif accession == "GSE32413":
        source = rows["source_name"]
        is_hc = source.str.contains("Normal Non-SSc", case=False, regex=False)
        is_ssc = source.str.contains("Rituximab Patient", case=False, regex=False)
        out.loc[is_hc, "condition"] = "HC"
        out.loc[is_ssc, "condition"] = "SSc"
        normal = source.str.extract(r"(NormalRIT_\d+)", flags=re.I, expand=False)
        dssc = source.str.extract(r"(dSSc_\d+)", flags=re.I, expand=False)
        rit = source.str.extract(r"\b(RIT\d+)_", flags=re.I, expand=False)
        out["subject_id"] = normal.fillna(dssc).fillna(rit).fillna("")
        title = rows["title"]
        exact_base = title.str.contains(r"(?<!PRE)BASE", case=False, regex=True)
        untimed = ~title.str.contains(r"(?:PREBASE|\d+\s*MOS)", case=False, regex=True)
        prebase = title.str.contains("PREBASE", case=False, regex=False)
        out["time_rank"] = np.select(
            [is_hc | exact_base | untimed, prebase],
            [0, 1],
            default=2,
        )
        invalid = ~(is_hc | is_ssc) | out["subject_id"].eq("")
        out.loc[invalid, "analysis_keep"] = False
        out.loc[invalid, "exclusion_reason"] = "unresolved disease or participant identity"
        for subject_id, group in out.loc[out["analysis_keep"]].groupby("subject_id"):
            minimum = group["time_rank"].min()
            later = group.index[group["time_rank"].gt(minimum)]
            out.loc[later, "analysis_keep"] = False
            out.loc[later, "exclusion_reason"] = "post-baseline biopsy"
    elif accession == "GSE76885":
        out["condition"] = rows["condition_screen"]
        out["subject_id"] = rows["subject_id"]
        is_hc = out["condition"].eq("HC")
        is_baseline_ssc = out["condition"].eq("SSc") & rows["title"].str.contains(
            r"(?:^|_)base(?:_|$)", case=False, regex=True
        )
        out["analysis_keep"] = is_hc | is_baseline_ssc
        out.loc[
            out["condition"].eq("localized_scleroderma"), "exclusion_reason"
        ] = "morphea/localized scleroderma"
        out.loc[
            ~out["analysis_keep"] & out["exclusion_reason"].eq(""),
            "exclusion_reason",
        ] = "post-baseline SSc biopsy"
    elif accession == "GSE231692":
        out["condition"] = rows["condition_screen"]
        out["subject_id"] = rows["subject_id"]
        is_hc = out["condition"].eq("HC")
        is_baseline_ssc = out["condition"].eq("SSc") & rows["timepoint_screen"].eq(
            "baseline"
        )
        out["analysis_keep"] = is_hc | is_baseline_ssc
        out.loc[~out["analysis_keep"], "exclusion_reason"] = (
            "post-baseline SSc biopsy"
        )
    elif accession in {"GSE76807", "GSE95065", "GSE125362"}:
        out["condition"] = rows["condition_screen"]
        out["subject_id"] = rows["title"]
        invalid = ~out["condition"].isin(["SSc", "HC"])
        out.loc[invalid, "analysis_keep"] = False
        out.loc[invalid, "exclusion_reason"] = "not adult SSc or healthy control"
    else:
        raise KeyError(accession)

    unresolved = out["condition"].eq("") | out["subject_id"].eq("")
    out.loc[unresolved, "analysis_keep"] = False
    out.loc[unresolved & out["exclusion_reason"].eq(""), "exclusion_reason"] = (
        "unresolved condition or participant"
    )
    kept = out.loc[out["analysis_keep"]]
    conflict = kept.groupby("subject_id")["condition"].nunique()
    if conflict.gt(1).any():
        raise RuntimeError(f"{accession}: subject has conflicting conditions")
    return out.reset_index(drop=True)


def load_gene_expression(accession: str) -> tuple[pd.DataFrame, dict[str, int | float | str]]:
    config = COHORTS[accession]
    if accession == "GSE231692":
        path = OUT / "GSE231692_edgeR_logCPM.csv"
        gene_qc_path = OUT / "GSE231692_edgeR_gene_qc.csv"
        if not path.exists() or not gene_qc_path.exists():
            raise FileNotFoundError(
                "Prepare GSE231692 with 07_prepare_gse231692_edgeR.R before "
                "running CollecTRI/ULM"
            )
        sample_by_gene = pd.read_csv(path, index_col=0)
        matrix = sample_by_gene.T
        gene_qc = pd.read_csv(gene_qc_path).set_index("metric")["value"]
        return matrix, {
            "accession": accession,
            "platform": config["platform"],
            "source_samples": 76,
            "source_probes": int(gene_qc["source_ensembl_rows"]),
            "mapped_genes": int(gene_qc["symbols_after_filterByExpr"]),
            "missing_values_filled": 0,
            "input_min": float(np.nanmin(matrix.to_numpy())),
            "input_max": float(np.nanmax(matrix.to_numpy())),
            "log2_x_plus_1_applied": False,
            "normalization": (
                "edgeR 4.4.2 filterByExpr(group=condition), "
                "TMM, logCPM(prior.count=2)"
            ),
            "ensembl_mapping": (
                "AnnotationDbi 1.68.0 and org.Hs.eg.db 3.20.0; "
                "counts summed by HGNC symbol before normalization"
            ),
            "count_matrix_sha256": sha256(GSE231692_COUNTS),
        }
    path = SOFT / f"{accession}_family.soft.gz"
    gse = GEOparse.get_GEO(filepath=str(path), silent=True)
    if list(gse.gpls) != [config["platform"]]:
        raise RuntimeError(
            f"{accession}: expected {config['platform']}, found {list(gse.gpls)}"
        )
    matrix = gse.pivot_samples("VALUE").apply(pd.to_numeric, errors="coerce")
    source_probes = int(matrix.shape[0])
    platform_table = gse.gpls[config["platform"]].table.copy()
    id_column = "ID" if "ID" in platform_table else platform_table.columns[0]
    annotation = platform_table.set_index(id_column)[config["map_column"]]
    annotation = annotation.reindex(matrix.index)
    if config["map_column"] == "ENTREZ_GENE_ID":
        symbols = entrez_map(annotation)
    else:
        symbols = annotation.map(clean_symbol)
    matrix.index = symbols
    matrix = matrix.loc[matrix.index.notna()]
    matrix = matrix.replace([np.inf, -np.inf], np.nan)
    matrix = matrix.groupby(level=0).mean()
    missing_before_fill = int(matrix.isna().sum().sum())
    matrix = matrix.T.fillna(matrix.mean(axis=1)).T.fillna(0.0)
    before_log_max = float(np.nanmax(matrix.to_numpy()))
    log_transform = bool(float(np.nanmin(matrix.to_numpy())) >= 0 and before_log_max > 100)
    if log_transform:
        matrix = np.log2(matrix.clip(lower=0) + 1)
    qc = {
        "accession": accession,
        "platform": config["platform"],
        "source_samples": len(gse.gsms),
        "source_probes": source_probes,
        "mapped_genes": int(matrix.shape[0]),
        "missing_values_filled": missing_before_fill,
        "input_min": float(np.nanmin(matrix.to_numpy())),
        "input_max": float(np.nanmax(matrix.to_numpy())),
        "log2_x_plus_1_applied": log_transform,
    }
    return matrix, qc


def collectri_network() -> pd.DataFrame:
    path = OUT / "collectri_human_decoupler_2_1_6.csv"
    if path.exists():
        return pd.read_csv(path)
    network = dc.op.collectri(organism="human")
    network.to_csv(path, index=False)
    return network


def run_cohort(
    accession: str,
    network: pd.DataFrame,
    force: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    activity_path = OUT / f"{accession}_subject_regulon_activity.csv"
    metadata_path = OUT / f"{accession}_subject_metadata.csv"
    qc_path = OUT / f"{accession}_processing_qc.json"
    if activity_path.exists() and metadata_path.exists() and qc_path.exists() and not force:
        return (
            pd.read_csv(activity_path, index_col=0),
            pd.read_csv(metadata_path, index_col=0),
            json.loads(qc_path.read_text(encoding="utf-8")),
        )

    mapping = subject_map(accession)
    mapping.to_csv(OUT / f"{accession}_sample_to_subject.csv", index=False)
    expression, qc = load_gene_expression(accession)
    kept = mapping.loc[mapping["analysis_keep"]].copy()
    common = kept["geo_accession"].loc[
        kept["geo_accession"].isin(expression.columns)
    ]
    if len(common) != len(kept):
        missing = sorted(set(kept["geo_accession"]) - set(expression.columns))
        raise RuntimeError(f"{accession}: kept samples absent from matrix: {missing}")
    kept = kept.set_index("geo_accession").loc[common]
    sample_by_gene = expression[common].T
    sample_by_gene["subject_id"] = kept.loc[sample_by_gene.index, "subject_id"]
    subject_by_gene = sample_by_gene.groupby("subject_id", sort=True).mean()
    subject_meta = (
        kept.groupby("subject_id")
        .agg(
            condition=("condition", "first"),
            n_source_samples=("condition", "size"),
        )
        .loc[subject_by_gene.index]
    )
    subject_meta["cohort"] = accession
    subject_meta["role"] = COHORTS[accession]["role"]

    usable_network = network.loc[
        network["target"].isin(subject_by_gene.columns),
        ["source", "target", "weight"],
    ].dropna()
    data = ad.AnnData(subject_by_gene.to_numpy(dtype=np.float32))
    data.obs_names = subject_by_gene.index.astype(str)
    data.var_names = subject_by_gene.columns.astype(str)
    dc.mt.ulm(data=data, net=usable_network, tmin=5, verbose=False)
    score_keys = [
        key
        for key in data.obsm
        if "ulm" in key.lower() and "pval" not in key.lower() and "padj" not in key.lower()
    ]
    if len(score_keys) != 1:
        raise RuntimeError(f"{accession}: ambiguous ULM score keys: {score_keys}")
    activity = pd.DataFrame(data.obsm[score_keys[0]], index=data.obs_names)
    if "HES1" not in activity.columns:
        raise RuntimeError(f"{accession}: HES1 absent after CollecTRI/ULM")
    activity.to_csv(activity_path)
    subject_meta.to_csv(metadata_path)

    qc.update(
        {
            "role": COHORTS[accession]["role"],
            "retained_subjects": int(len(subject_meta)),
            "retained_SSc": int(subject_meta["condition"].eq("SSc").sum()),
            "retained_HC": int(subject_meta["condition"].eq("HC").sum()),
            "retained_source_samples": int(subject_meta["n_source_samples"].sum()),
            "collectri_edges_after_gene_intersection": int(len(usable_network)),
            "collectri_sources_after_tmin": int(activity.shape[1]),
            "ulm_tmin": 5,
            "ulm_score_key": score_keys[0],
            "hes1_present": True,
            "source_sha256": (
                sha256(GSE231692_COUNTS)
                if accession == "GSE231692"
                else sha256(SOFT / f"{accession}_family.soft.gz")
            ),
        }
    )
    qc_path.write_text(json.dumps(qc, indent=2), encoding="utf-8")
    return activity, subject_meta, qc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cohorts",
        nargs="+",
        choices=sorted(COHORTS),
        default=list(COHORTS),
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    network = collectri_network()
    activities: list[pd.DataFrame] = []
    metadata: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    for accession in args.cohorts:
        print(f"[external-bulk] {accession}", flush=True)
        activity, meta, qc = run_cohort(accession, network, args.force)
        activity = activity.copy()
        activity.insert(0, "subject_id", activity.index)
        activity.insert(1, "cohort", accession)
        activities.append(activity.reset_index(drop=True))
        metadata.append(meta.reset_index(names="subject_id"))
        qc_rows.append(qc)
        gc.collect()

    pd.concat(activities, ignore_index=True).to_csv(
        OUT / "external_subject_regulon_activity.csv", index=False
    )
    pd.concat(metadata, ignore_index=True).to_csv(
        OUT / "external_subject_metadata.csv", index=False
    )
    pd.DataFrame(qc_rows).to_csv(OUT / "external_processing_qc.csv", index=False)
    manifest = {
        "python": platform.python_version(),
        "GEOparse": getattr(GEOparse, "__version__", "unknown"),
        "decoupler": dc.__version__,
        "anndata": ad.__version__,
        "mygene": getattr(mygene, "__version__", "unknown"),
        "cohorts": args.cohorts,
        "independent_unit": "participant",
        "within_participant_aggregation": (
            "mean normalized expression across technical replicates or biopsy sites "
            "within the selected baseline/initial time point"
        ),
        "network": "decoupler.op.collectri(organism='human')",
        "activity_method": "decoupler.mt.ulm with tmin=5",
        "outcome_blinding": (
            "cohort, tissue, condition, time point and participant eligibility fixed "
            "before HES1 or model outcomes"
        ),
    }
    (OUT / "external_feature_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"Wrote external donor-level regulons to {OUT}", flush=True)


if __name__ == "__main__":
    main()
