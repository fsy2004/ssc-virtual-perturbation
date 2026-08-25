#!/usr/bin/env python3
"""Aggregate prespecified repeated CellOracle nulls and apply BH correction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    summaries = sorted(
        args.root.glob("*/*/repeated_randomized_grn_summary.json")
    )
    if not summaries:
        raise FileNotFoundError(f"No repeated-null summaries under {args.root}")
    rows = []
    for path in summaries:
        row = json.loads(path.read_text(encoding="utf-8"))
        row["summary_path"] = str(path.resolve())
        rows.append(row)
    result = pd.DataFrame(rows)
    if result.duplicated(["arm", "gene"]).any():
        raise RuntimeError("Duplicate repeated-null summary for the same arm/gene")
    _, q_values, _, _ = multipletests(
        result["empirical_plus_one_p"].to_numpy(dtype=float),
        method="fdr_bh",
    )
    result["empirical_bh_q"] = q_values
    result = result.sort_values(
        ["empirical_bh_q", "empirical_plus_one_p", "arm", "gene"]
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
