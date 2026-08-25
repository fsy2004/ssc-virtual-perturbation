#!/usr/bin/env python3
"""Link archived SCENIC+ eRegulon regions to condition-aware donor DAR.

The archived eRegulons and the new DAR share the GSE312129 multiome source.
This script therefore quantifies condition-aware accessibility support; it
does not label the result as independent validation.

For every TF and edge scope, the primary competition statistic is the median
SSc-vs-HC logFC across tested eRegulon regions. Its empirical one-sided P value
is calculated against random region sets matched on logCPM and region length.
Matching bins and resampling counts are fixed before TF labels are inspected.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import re
import zlib
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260726
PRIORITY_TFS = ("HES1", "RBPJ", "SMAD3")
REGION_PATTERN = re.compile(
    r"^(?P<chrom>chr[^:\s_-]+)[:_](?P<start>\d+)[-_](?P<end>\d+)$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-meta", type=Path, required=True)
    parser.add_argument("--extended-meta", type=Path, required=True)
    parser.add_argument("--dar", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--scopes",
        nargs="+",
        choices=("direct", "extended"),
        default=("direct", "extended"),
    )
    parser.add_argument("--n-permutations", type=int, default=9_999)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--min-regions", type=int, default=10)
    parser.add_argument("--match-bins", type=int, default=10)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_region(value: object) -> tuple[str, int]:
    text = str(value).strip()
    match = REGION_PATTERN.match(text)
    if match is None:
        fields = re.split(r"\s+", text)
        if len(fields) >= 3 and fields[0].startswith("chr"):
            chrom, start_text, end_text = fields[:3]
            start, end = int(start_text), int(end_text)
        else:
            raise ValueError(f"Cannot parse genomic region: {text!r}")
    else:
        chrom = match.group("chrom")
        start, end = int(match.group("start")), int(match.group("end"))
    if start < 0 or end <= start:
        raise ValueError(f"Invalid genomic interval: {text!r}")
    return f"{chrom}:{start}-{end}", end - start


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    ordered_index = valid.index.to_numpy()[order]
    result.loc[ordered_index] = adjusted
    return result


def quantile_bins(values: pd.Series, n_bins: int) -> pd.Series:
    unique = int(values.nunique(dropna=True))
    bins = min(n_bins, unique)
    if bins < 2:
        return pd.Series(0, index=values.index, dtype=int)
    ranked = values.rank(method="first")
    return pd.qcut(ranked, q=bins, labels=False, duplicates="drop").astype(int)


def load_dar(path: Path, match_bins: int) -> pd.DataFrame:
    frame = pd.read_csv(require_file(path))
    required = {"region", "logFC", "logCPM", "FDR"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"DAR table lacks columns: {sorted(missing)}")
    parsed = frame["region"].map(canonical_region)
    frame["region"] = parsed.map(lambda value: value[0])
    frame["region_length"] = parsed.map(lambda value: value[1])
    if frame["region"].duplicated().any():
        raise ValueError("DAR table contains duplicate canonical regions")
    numeric = ["logFC", "logCPM", "FDR", "region_length"]
    frame[numeric] = frame[numeric].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame[["logFC", "logCPM"]].to_numpy()).all():
        raise ValueError("DAR logFC/logCPM contains non-finite values")
    if not frame["FDR"].between(0, 1).all():
        raise ValueError("DAR FDR values fall outside [0,1]")
    frame["abundance_bin"] = quantile_bins(frame["logCPM"], match_bins)
    frame["length_bin"] = quantile_bins(frame["region_length"], match_bins)
    frame["match_bin"] = (
        frame["abundance_bin"].astype(str)
        + "_"
        + frame["length_bin"].astype(str)
    )
    return frame


def load_scope(path: Path, scope: str) -> pd.DataFrame:
    frame = pd.read_csv(require_file(path))
    required = {"TF", "Region"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {sorted(missing)}")
    parsed = frame["Region"].map(canonical_region)
    frame["region"] = parsed.map(lambda value: value[0])
    frame["scope"] = scope
    keep = ["scope", "TF", "region"]
    if "eRegulon_name" in frame.columns:
        keep.append("eRegulon_name")
    return frame[keep].drop_duplicates()


def matched_null(
    linked: pd.DataFrame,
    background: pd.DataFrame,
    n_permutations: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    strata_counts = linked["match_bin"].value_counts().to_dict()
    pools = {
        stratum: background.loc[
            background["match_bin"].eq(stratum), "logFC"
        ].to_numpy(dtype=float)
        for stratum in strata_counts
    }
    null = np.empty(n_permutations, dtype=float)
    chunk_size = 500
    for start in range(0, n_permutations, chunk_size):
        stop = min(start + chunk_size, n_permutations)
        n_draws = stop - start
        draws = [
            rng.choice(
                pools[stratum],
                size=(n_draws, count),
                replace=True,
            )
            for stratum, count in strata_counts.items()
        ]
        null[start:stop] = np.median(np.concatenate(draws, axis=1), axis=1)
    return null


def summarize_tf(
    task: tuple[str, str, pd.DataFrame, pd.DataFrame, int, int, int]
) -> dict[str, object]:
    scope, tf, linked, background, n_permutations, min_regions, seed = task
    observed = float(linked["logFC"].median())
    result: dict[str, object] = {
        "scope": scope,
        "TF": tf,
        "n_eRegulon_regions": int(linked["n_eRegulon_regions"].iloc[0]),
        "n_DAR_tested_regions": int(len(linked)),
        "tested_region_fraction": float(
            len(linked) / linked["n_eRegulon_regions"].iloc[0]
        ),
        "SSc_positive_region_fraction": float((linked["logFC"] > 0).mean()),
        "median_logFC": observed,
        "mean_logFC": float(linked["logFC"].mean()),
        "n_positive_logFC": int((linked["logFC"] > 0).sum()),
        "n_FDR_lt_0_05": int((linked["FDR"] < 0.05).sum()),
        "empirical_competition_p": np.nan,
        "matched_null_median": np.nan,
        "matched_null_q025": np.nan,
        "matched_null_q975": np.nan,
        "median_logFC_minus_null": np.nan,
        "competition_status": "insufficient_tested_regions",
    }
    if len(linked) < min_regions:
        return result
    null = matched_null(linked, background, n_permutations, seed)
    result.update(
        {
            "empirical_competition_p": float(
                (1 + np.count_nonzero(null >= observed))
                / (n_permutations + 1)
            ),
            "matched_null_median": float(np.median(null)),
            "matched_null_q025": float(np.quantile(null, 0.025)),
            "matched_null_q975": float(np.quantile(null, 0.975)),
            "median_logFC_minus_null": float(observed - np.median(null)),
            "competition_status": "tested",
        }
    )
    return result


def main() -> None:
    args = parse_args()
    if args.n_permutations < 999:
        raise ValueError("At least 999 matched permutations are required")
    if args.n_jobs < 1 or args.min_regions < 2 or args.match_bins < 2:
        raise ValueError("n-jobs, min-regions and match-bins are invalid")

    args.outdir.mkdir(parents=True, exist_ok=True)
    dar = load_dar(args.dar, args.match_bins)
    scope_paths = {
        "direct": args.direct_meta,
        "extended": args.extended_meta,
    }
    egrn = pd.concat(
        [load_scope(scope_paths[scope], scope) for scope in args.scopes],
        ignore_index=True,
    )
    region_counts = (
        egrn.groupby(["scope", "TF"], observed=True)["region"]
        .nunique()
        .rename("n_eRegulon_regions")
        .reset_index()
    )
    egrn_unique = egrn[["scope", "TF", "region"]].drop_duplicates()
    linked = egrn_unique.merge(
        dar,
        on="region",
        how="inner",
        validate="many_to_one",
    ).merge(
        region_counts,
        on=["scope", "TF"],
        how="left",
        validate="many_to_one",
    )
    if linked.empty:
        raise RuntimeError("No SCENIC+ regions overlap the donor-DAR regions")

    tasks = []
    for (scope, tf), group in linked.groupby(["scope", "TF"], observed=True):
        tf_seed = (
            args.seed
            + zlib.crc32(f"{scope}:{tf}".encode("utf-8"))
        ) % (2**32 - 1)
        tasks.append(
            (
                str(scope),
                str(tf),
                group.copy(),
                dar,
                args.n_permutations,
                args.min_regions,
                tf_seed,
            )
        )

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.n_jobs
    ) as executor:
        rows = list(executor.map(summarize_tf, tasks))
    summary = pd.DataFrame(rows)
    summary["empirical_competition_bh_q"] = bh_adjust(
        summary["empirical_competition_p"]
    )
    summary["priority_TF"] = summary["TF"].isin(PRIORITY_TFS)
    summary = summary.sort_values(
        [
            "scope",
            "empirical_competition_bh_q",
            "median_logFC_minus_null",
            "TF",
        ],
        ascending=[True, True, False, True],
        na_position="last",
    )
    summary.to_csv(args.outdir / "all_TF_condition_accessibility_support.csv", index=False)
    summary.loc[summary["priority_TF"]].to_csv(
        args.outdir / "HES1_RBPJ_SMAD3_condition_accessibility_support.csv",
        index=False,
    )

    linked_detail = linked.loc[
        linked["TF"].isin(PRIORITY_TFS)
    ].merge(
        egrn,
        on=["scope", "TF", "region"],
        how="left",
    )
    linked_detail.to_csv(
        args.outdir / "HES1_RBPJ_SMAD3_linked_DAR_regions.csv.gz",
        index=False,
        compression="gzip",
    )
    manifest = {
        "analysis": "condition_aware_accessibility_support_for_archived_eRegulons",
        "shared_source_boundary": (
            "Archived eRegulons and donor DAR use GSE312129; this is "
            "condition-aware shared-source support, not independent validation."
        ),
        "competition_statistic": "median_SSc_vs_HC_logFC",
        "competition_null": (
            "bootstrap region draws with replacement, matched on logCPM and "
            "region-length quantile bins"
        ),
        "competition_alternative": "TF_region_median_logFC_greater_than_matched_null",
        "n_permutations": args.n_permutations,
        "match_bins": args.match_bins,
        "min_tested_regions": args.min_regions,
        "multiple_testing": "BH_across_all_TF_scope_tests",
        "priority_TFs_prespecified": list(PRIORITY_TFS),
        "seed": args.seed,
        "n_DAR_regions": int(len(dar)),
        "n_TF_scope_tests": int(len(summary)),
        "inputs": {
            str(path.resolve()): sha256(path)
            for path in (args.direct_meta, args.extended_meta, args.dar)
        },
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
