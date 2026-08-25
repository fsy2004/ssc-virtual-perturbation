#!/usr/bin/env python3
"""Audit the prespecified one-factor CellOracle sensitivity matrix.

The primary configuration remains the inferential reference. This script
summarises all TFs under every frozen sensitivity setting and never selects a
setting from the HES1 result.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


VARIANTS = {
    "propagation2": "n_propagation=2",
    "propagation4": "n_propagation=4",
    "neighbours100": "vf_n_neighbors=100",
    "neighbours300": "vf_n_neighbors=300",
    "hvg3000": "n_top_genes=3000",
    "links1000": "link_topn=1000",
    "links3000": "link_topn=3000",
}
PRESPECIFIED_GENES = ("HES1", "SMAD3", "FOSL2", "RUNX1")


def read_ranking(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    ranking = pd.read_csv(path)
    required = {"gene", "rank", "n_ranked", "rank_fraction", "ps_score", "ps_qbh"}
    missing = required.difference(ranking.columns)
    if missing:
        raise RuntimeError(f"{path} lacks columns: {sorted(missing)}")
    if ranking["gene"].duplicated().any():
        raise RuntimeError(f"{path} contains duplicate genes")
    if not (ranking["n_ranked"] == len(ranking)).all():
        raise RuntimeError(f"{path} has inconsistent n_ranked values")
    if len(ranking) <= 43:
        raise RuntimeError(f"{path} did not retain an all-active denominator")
    absent = {"HES1", "SMAD3"}.difference(ranking["gene"])
    if absent:
        raise RuntimeError(f"{path} lacks required genes: {sorted(absent)}")
    return ranking.sort_values("rank").reset_index(drop=True)


def gene_values(ranking: pd.DataFrame, gene: str) -> dict[str, object]:
    match = ranking.loc[ranking["gene"] == gene]
    if match.empty:
        return {
            f"{gene}_rank": np.nan,
            f"{gene}_rank_fraction": np.nan,
            f"{gene}_ps_score": np.nan,
            f"{gene}_ps_qbh": np.nan,
            f"{gene}_negative_fdr": False,
        }
    row = match.iloc[0]
    return {
        f"{gene}_rank": int(row["rank"]),
        f"{gene}_rank_fraction": float(row["rank_fraction"]),
        f"{gene}_ps_score": float(row["ps_score"]),
        f"{gene}_ps_qbh": float(row["ps_qbh"]),
        f"{gene}_negative_fdr": bool(
            float(row["ps_score"]) < 0 and float(row["ps_qbh"]) < 0.05
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--sensitivity-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    primary = read_ranking(args.primary)
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_long = [primary.assign(variant="primary", parameter="frozen_primary")]
    summary_rows: list[dict[str, object]] = []

    primary_top_n = max(1, math.ceil(0.1 * len(primary)))
    primary_top = set(primary.nsmallest(primary_top_n, "rank")["gene"])

    for variant, parameter in VARIANTS.items():
        variant_root = args.sensitivity_root / variant
        ranking_path = variant_root / "KO_ranking_all_active_skinatac.csv"
        ranking = read_ranking(ranking_path)
        all_long.append(ranking.assign(variant=variant, parameter=parameter))

        common = primary[["gene", "ps_score"]].merge(
            ranking[["gene", "ps_score"]],
            on="gene",
            suffixes=("_primary", "_variant"),
            validate="one_to_one",
        )
        rho = float(
            spearmanr(
                common["ps_score_primary"],
                common["ps_score_variant"],
            ).statistic
        )
        sign_agreement = float(
            np.mean(
                np.sign(common["ps_score_primary"])
                == np.sign(common["ps_score_variant"])
            )
        )
        variant_top_n = max(1, math.ceil(0.1 * len(ranking)))
        variant_top = set(ranking.nsmallest(variant_top_n, "rank")["gene"])
        union = primary_top | variant_top
        overlap = primary_top & variant_top

        failure_path = variant_root / "KO_failures_all_active_skinatac.csv"
        failures = pd.read_csv(failure_path) if failure_path.is_file() else pd.DataFrame()
        row: dict[str, object] = {
            "variant": variant,
            "parameter": parameter,
            "n_ranked": int(len(ranking)),
            "n_common_with_primary": int(len(common)),
            "n_failed": int(len(failures)),
            "all_TF_score_spearman_vs_primary": rho,
            "score_sign_agreement_vs_primary": sign_agreement,
            "primary_top_decile_n": primary_top_n,
            "variant_top_decile_n": variant_top_n,
            "top_decile_overlap_n": int(len(overlap)),
            "top_decile_jaccard": float(len(overlap) / len(union)),
            "top_decile_shared_genes": ";".join(sorted(overlap)),
        }
        for gene in PRESPECIFIED_GENES:
            row.update(gene_values(ranking, gene))

        recovery_path = (
            variant_root / "baseline_recovery_summary_all_active_skinatac.json"
        )
        if recovery_path.is_file():
            recovery = json.loads(recovery_path.read_text(encoding="utf-8"))
            for key in (
                "ps_auroc_recover_drivers",
                "de_baseline_auroc",
                "ps_precision_at_10",
                "de_precision_at_10",
                "perturbation_beats_DE_baseline",
            ):
                row[key] = recovery.get(key)
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    if not summary["SMAD3_negative_fdr"].all():
        failed = summary.loc[~summary["SMAD3_negative_fdr"], "variant"].tolist()
        raise RuntimeError(f"SMAD3 positive-control gate failed for: {failed}")

    long_frame = pd.concat(all_long, ignore_index=True)
    long_frame.to_csv(
        args.outdir / "CellOracle_sensitivity_all_TF_long.csv",
        index=False,
    )
    summary.to_csv(
        args.outdir / "CellOracle_sensitivity_summary.csv",
        index=False,
    )
    matrix = long_frame.pivot(
        index="gene",
        columns="variant",
        values=["rank", "rank_fraction", "ps_score", "ps_qbh"],
    )
    matrix.columns = [f"{metric}__{variant}" for metric, variant in matrix.columns]
    matrix.reset_index().to_csv(
        args.outdir / "CellOracle_sensitivity_all_TF_matrix.csv",
        index=False,
    )
    manifest = {
        "primary": str(args.primary.resolve()),
        "sensitivity_root": str(args.sensitivity_root.resolve()),
        "variants": VARIANTS,
        "prespecified_genes": list(PRESPECIFIED_GENES),
        "selection_rule": "primary configuration retained; no HES1-guided selection",
        "all_variants_passed_smad3_control": True,
    }
    (args.outdir / "CellOracle_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
