#!/usr/bin/env python3
"""Condition-aware SSc-skin extension for the completed SCENIC+ analysis.

This script creates a new analysis root. It never writes into the completed
``/data/ssc/scenicplus`` archive.

The biological adaptation is fixed before results are inspected:

1. Peaks are called for cell-type-by-condition pseudobulks and then unioned.
2. A region is retained when it is accessible in at least two donors within
   either SSc or HC. Requiring support in both conditions would remove genuine
   disease-restricted accessibility.
3. Efficiency subsampling is stratified by condition, cell type and donor.
4. Fibroblast differential accessibility is exported from every fibroblast
   barcode in every donor. The 10,000-cell efficiency subset is never used
   for donor differential-accessibility counts.

MACS2 and pycisTopic numerical settings remain the SCENIC+ method defaults.
They are technical settings, not tuned SSc parameters.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SEED = 20260726
SEP = "___"
FROZEN = {
    "design_version": "ssc_skin_condition_aware_20260726_v2",
    "peak_grouping": ["cell_type", "condition"],
    "donor_support_rule": "accessible_in_at_least_2_donors_within_either_condition",
    "min_donors_within_condition": 2,
    "subsample_strata": ["condition", "cell_type", "sample_id"],
    "target_cells": 10_000,
    "fibroblast_label": "Fibroblast",
    "dar_cell_rule": "all_fibroblast_barcodes_per_donor",
    "dar_region_rule": (
        "condition_union_peak_accessible_in_at_least_2_fibroblast_donors_"
        "within_either_condition"
    ),
    "macs2": {
        "input_format": "BEDPE",
        "shift": 73,
        "ext_size": 146,
        "keep_dup": "all",
        "q_value": 0.05,
        "genome_size": "hs",
    },
    "consensus_peak_half_width": 250,
    "cistopic_partition": 20,
    "seed": SEED,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("manifest", "peaks", "cistopic", "dar-input", "all"),
        default="manifest",
    )
    parser.add_argument(
        "--adata",
        type=Path,
        default=Path("/data/ssc/scenicplus/rna/adata_annot.h5ad"),
    )
    parser.add_argument(
        "--fragments-dir",
        type=Path,
        default=Path("/data/ssc/basegrn/frags"),
    )
    parser.add_argument(
        "--chromsizes",
        type=Path,
        default=Path("/data/ssc/scenicplus/ref/hg38.chrom.sizes"),
    )
    parser.add_argument(
        "--blacklist",
        type=Path,
        default=Path("/data/ssc/scenicplus/ref/hg38-blacklist.v2.bed"),
    )
    parser.add_argument("--macs2", default="macs2")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("/data/ssc/scenicplus_ssc_skin_adapted_20260726_v2"),
    )
    parser.add_argument("--n-cpu", type=int, default=6)
    parser.add_argument(
        "--sample-workers",
        type=int,
        default=1,
        help="Concurrent donor-level cisTopic builders; use a memory-safe value.",
    )
    parser.add_argument(
        "--sample-cpu",
        type=int,
        default=1,
        help="CPUs passed to each donor-level cisTopic builder.",
    )
    parser.add_argument("--target-cells", type=int, default=FROZEN["target_cells"])
    parser.add_argument(
        "--min-donors-within-condition",
        type=int,
        default=FROZEN["min_donors_within_condition"],
    )
    parser.add_argument("--fibroblast-label", default=FROZEN["fibroblast_label"])
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--hash-inputs", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, hash_input: bool) -> dict[str, object]:
    require_file(path)
    stat = path.stat()
    record: dict[str, object] = {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if hash_input:
        record["sha256"] = sha256(path)
    return record


def condition_from_sample(sample_id: str) -> str:
    value = str(sample_id).strip().upper()
    if value.startswith("SSC"):
        return "SSc"
    if value.startswith("HC"):
        return "HC"
    raise ValueError(f"Cannot derive SSc/HC condition from sample_id={sample_id!r}")


def safe_group(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def load_metadata(adata_path: Path) -> pd.DataFrame:
    require_file(adata_path)
    import scanpy as sc

    adata = sc.read_h5ad(adata_path, backed="r")
    required = {"sample_id", "cell_type"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise KeyError(f"{adata_path} lacks required obs columns: {sorted(missing)}")
    meta = adata.obs[["sample_id", "cell_type"]].copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta["cell_type"] = meta["cell_type"].astype(str)
    meta["condition"] = meta["sample_id"].map(condition_from_sample)
    meta["peak_group"] = [
        safe_group(f"{cell_type}__{condition}")
        for cell_type, condition in zip(meta["cell_type"], meta["condition"])
    ]
    if getattr(adata, "file", None) is not None:
        adata.file.close()
    if set(meta["condition"]) != {"HC", "SSc"}:
        raise ValueError(
            f"Both HC and SSc are required; observed {sorted(meta['condition'].unique())}"
        )
    if not meta.index.is_unique:
        raise ValueError("AnnData cell identifiers are not unique")
    return meta


def discover_fragments(fragments_dir: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    pattern = re.compile(r"GSM\d+_([A-Za-z]+\d+)_atac_fragments\.tsv\.gz$")
    for path in sorted(fragments_dir.glob("GSM*_atac_fragments.tsv.gz")):
        match = pattern.search(path.name)
        if match:
            paths[match.group(1)] = str(path.resolve())
    if not paths:
        raise FileNotFoundError(f"No multiome fragment files found under {fragments_dir}")
    return paths


def largest_remainder_alloc(counts: pd.Series, target: int) -> pd.Series:
    if target <= 0:
        raise ValueError("target-cells must be positive")
    counts = counts.astype(int)
    target = min(target, int(counts.sum()))
    raw = counts / counts.sum() * target
    allocation = np.floor(raw).astype(int)
    allocation[(counts > 0) & (allocation == 0)] = 1
    allocation = np.minimum(allocation, counts)

    while int(allocation.sum()) > target:
        candidates = allocation[allocation > 1]
        if candidates.empty:
            raise RuntimeError("Cannot reduce stratified allocation to requested target")
        idx = (allocation[candidates.index] - raw[candidates.index]).idxmax()
        allocation.loc[idx] -= 1

    remainder = raw - allocation
    while int(allocation.sum()) < target:
        candidates = counts.index[allocation < counts]
        if len(candidates) == 0:
            break
        idx = remainder[candidates].idxmax()
        allocation.loc[idx] += 1
        remainder.loc[idx] -= 1
    return allocation.astype(int)


def write_manifest(args: argparse.Namespace, meta: pd.DataFrame) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_dir = args.outdir / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    frozen_path = manifest_dir / "frozen_design.json"
    cells_path = manifest_dir / f"stratified_cells_{args.target_cells}.tsv"
    counts_path = manifest_dir / f"stratified_cell_counts_{args.target_cells}.tsv"
    peak_counts_path = manifest_dir / "peak_group_counts.tsv"
    if frozen_path.exists():
        existing = json.loads(frozen_path.read_text(encoding="utf-8"))
        expected = {
            "design_version": FROZEN["design_version"],
            "dar_cell_rule": FROZEN["dar_cell_rule"],
            "dar_region_rule": FROZEN["dar_region_rule"],
            "target_cells": args.target_cells,
            "min_donors_within_condition": args.min_donors_within_condition,
            "fibroblast_label": args.fibroblast_label,
            "seed": args.seed,
            "output_root": str(args.outdir.resolve()),
        }
        conflict = {
            key: (existing.get(key), value)
            for key, value in expected.items()
            if existing.get(key) != value
        }
        existing_input = existing.get("input", {}).get("path")
        expected_input = str(args.adata.resolve())
        if existing_input != expected_input:
            conflict["input.path"] = (existing_input, expected_input)
        if conflict:
            raise RuntimeError(
                "Existing frozen manifest conflicts with requested settings: "
                f"{conflict}. Use a new output root; do not overwrite the frozen run."
            )
        if (
            not args.force
            and cells_path.exists()
            and counts_path.exists()
            and peak_counts_path.exists()
        ):
            print(f"[manifest] exists, keeping frozen design: {frozen_path}")
            return
        print("[manifest] frozen settings match but manifest tables are incomplete; "
              "regenerating the missing tables")

    group_cols = ["condition", "cell_type", "sample_id"]
    counts = meta.groupby(group_cols, observed=True).size().rename("n_cells")
    allocation = largest_remainder_alloc(counts, args.target_cells)
    rng = np.random.default_rng(args.seed)
    selected: list[str] = []
    for stratum, n_select in allocation.items():
        mask = np.ones(len(meta), dtype=bool)
        for column, value in zip(group_cols, stratum):
            mask &= meta[column].to_numpy() == value
        cells = meta.index.to_numpy()[mask]
        selected.extend(rng.choice(cells, size=int(n_select), replace=False).tolist())

    selected_meta = meta.loc[selected].copy()
    selected_meta.to_csv(cells_path, sep="\t")
    (
        selected_meta.groupby(group_cols, observed=True)
        .size()
        .rename("n_selected")
        .reset_index()
        .merge(counts.reset_index(), on=group_cols, how="left")
        .to_csv(counts_path, sep="\t", index=False)
    )
    (
        meta.groupby(["peak_group", "cell_type", "condition"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(peak_counts_path, sep="\t", index=False)
    )

    frozen = dict(FROZEN)
    frozen.update(
        {
            "target_cells": args.target_cells,
            "min_donors_within_condition": args.min_donors_within_condition,
            "fibroblast_label": args.fibroblast_label,
            "seed": args.seed,
            "input": file_record(args.adata, args.hash_inputs),
            "output_root": str(args.outdir.resolve()),
        }
    )
    with frozen_path.open("w", encoding="utf-8") as handle:
        json.dump(frozen, handle, indent=2, sort_keys=True)
    print(f"[manifest] wrote frozen design and {len(selected_meta)} selected cells")


def load_stratified_meta(
    args: argparse.Namespace, full_meta: pd.DataFrame
) -> pd.DataFrame:
    cells_path = (
        args.outdir / "manifest" / f"stratified_cells_{args.target_cells}.tsv"
    )
    require_file(cells_path)
    selected = pd.read_table(cells_path, index_col=0)
    if not selected.index.is_unique:
        raise ValueError(f"Duplicate cell identifiers in {cells_path}")
    missing = selected.index.difference(full_meta.index)
    if len(missing):
        raise KeyError(
            f"{len(missing)} frozen cells are absent from the current AnnData; "
            f"examples: {missing[:5].tolist()}"
        )
    selected_meta = full_meta.loc[selected.index].copy()
    if len(selected_meta) != min(args.target_cells, len(full_meta)):
        raise RuntimeError(
            f"Frozen cell manifest has {len(selected_meta)} cells, expected "
            f"{min(args.target_cells, len(full_meta))}"
        )
    return selected_meta


def chromsizes_as_pyranges(path: Path):
    import pyranges as pr

    frame = pd.read_table(path, header=None, names=["Chromosome", "End"])
    frame["Start"] = 0
    return pr.PyRanges(frame[["Chromosome", "Start", "End"]])


def stage_peaks(args: argparse.Namespace, meta: pd.DataFrame) -> Path:
    from pycisTopic.iterative_peak_calling import get_consensus_peaks
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling

    out = args.outdir / "condition_peak_calling"
    consensus_path = out / "consensus_celltype_by_condition.bed"
    if consensus_path.exists() and not args.force:
        print(f"[peaks] exists, keeping: {consensus_path}")
        return consensus_path

    fragments = discover_fragments(args.fragments_dir)
    missing = sorted(set(meta["sample_id"]) - set(fragments))
    if missing:
        raise FileNotFoundError(f"Missing fragment files for samples: {missing}")

    out.mkdir(parents=True, exist_ok=True)
    bed_dir = out / "pseudobulk_bed"
    bw_dir = out / "pseudobulk_bw"
    macs_dir = out / "MACS"
    for path in (bed_dir, bw_dir, macs_dir):
        path.mkdir(parents=True, exist_ok=True)

    chromsizes = chromsizes_as_pyranges(args.chromsizes)
    _, bed_paths = export_pseudobulk(
        input_data=meta,
        variable="peak_group",
        sample_id_col="sample_id",
        chromsizes=chromsizes,
        bed_path=str(bed_dir),
        bigwig_path=str(bw_dir),
        path_to_fragments=fragments,
        n_cpu=args.n_cpu,
        normalize_bigwig=True,
        split_pattern=SEP,
        temp_dir=str(out / "tmp"),
    )
    narrow = peak_calling(
        macs_path=args.macs2,
        bed_paths=bed_paths,
        outdir=str(macs_dir),
        genome_size=FROZEN["macs2"]["genome_size"],
        n_cpu=args.n_cpu,
        input_format=FROZEN["macs2"]["input_format"],
        shift=FROZEN["macs2"]["shift"],
        ext_size=FROZEN["macs2"]["ext_size"],
        keep_dup=FROZEN["macs2"]["keep_dup"],
        q_value=FROZEN["macs2"]["q_value"],
    )
    consensus = get_consensus_peaks(
        narrow_peaks_dict=narrow,
        peak_half_width=FROZEN["consensus_peak_half_width"],
        chromsizes=chromsizes,
        path_to_blacklist=str(args.blacklist),
    )
    consensus.to_bed(str(consensus_path))
    print(f"[peaks] {len(consensus)} union peaks -> {consensus_path}")
    return consensus_path


def _build_sample_cistopic(task: dict[str, object]) -> tuple[str, int]:
    """Build one donor object in an isolated process and checkpoint it to disk."""
    from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

    output_path = Path(str(task["output_path"]))
    if output_path.exists() and not bool(task["force"]):
        with output_path.open("rb") as handle:
            obj = pickle.load(handle)
        return str(output_path), len(obj.cell_names)

    obj = create_cistopic_object_from_fragments(
        path_to_fragments=str(task["fragment_path"]),
        path_to_regions=str(task["consensus_path"]),
        path_to_blacklist=str(task["blacklist"]),
        valid_bc=list(task["raw_barcodes"]),
        n_cpu=int(task["n_cpu"]),
        project=str(task["sample_id"]),
        split_pattern=SEP,
        partition=FROZEN["cistopic_partition"],
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".part")
    with temporary.open("wb") as handle:
        pickle.dump(obj, handle)
    os.replace(temporary, output_path)
    return str(output_path), len(obj.cell_names)


def stage_cistopic(
    args: argparse.Namespace, meta: pd.DataFrame, consensus_path: Path
) -> Path:
    from pycisTopic.cistopic_class import merge

    out = args.outdir / "cistopic"
    output_path = out / "cistopic_condition_union_donor_supported.pkl"
    if output_path.exists() and not args.force:
        with output_path.open("rb") as handle:
            existing = pickle.load(handle)
        existing_cells = getattr(existing, "cell_names", None)
        existing_regions = getattr(existing, "region_names", None)
        if (
            existing_cells is None
            or existing_regions is None
            or len(existing_cells) == 0
            or len(existing_regions) == 0
        ):
            raise RuntimeError(f"Existing cisTopic object is incomplete: {output_path}")
        print(f"[cistopic] exists, keeping: {output_path}")
        return output_path

    out.mkdir(parents=True, exist_ok=True)
    fragments = discover_fragments(args.fragments_dir)
    sample_dir = out / "per_sample"
    tasks = []
    for sample_id in sorted(meta["sample_id"].unique()):
        sample_meta = meta[meta["sample_id"] == sample_id]
        raw_barcodes = [name.split(SEP)[0] for name in sample_meta.index]
        tasks.append(
            {
                "sample_id": sample_id,
                "fragment_path": fragments[sample_id],
                "consensus_path": str(consensus_path),
                "blacklist": str(args.blacklist),
                "raw_barcodes": raw_barcodes,
                "n_cpu": args.sample_cpu,
                "output_path": str(sample_dir / f"{safe_group(sample_id)}.pkl"),
                "force": args.force,
            }
        )

    completed: dict[str, tuple[str, int]] = {}
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.sample_workers
    ) as executor:
        future_to_sample = {
            executor.submit(_build_sample_cistopic, task): str(task["sample_id"])
            for task in tasks
        }
        for future in concurrent.futures.as_completed(future_to_sample):
            sample_id = future_to_sample[future]
            path, n_cells = future.result()
            completed[sample_id] = (path, n_cells)
            print(f"[cistopic] {sample_id}: {n_cells} cells -> {path}", flush=True)

    objects = []
    for sample_id in sorted(completed):
        with Path(completed[sample_id][0]).open("rb") as handle:
            objects.append(pickle.load(handle))

    cistopic = merge(objects, split_pattern=SEP) if len(objects) > 1 else objects[0]
    object_cells = pd.Index(cistopic.cell_names)
    absent_from_meta = object_cells.difference(meta.index)
    absent_from_object = meta.index.difference(object_cells)
    if len(absent_from_meta) or len(absent_from_object):
        raise RuntimeError(
            "cisTopic/manifest cell identifiers do not match exactly: "
            f"object_only={len(absent_from_meta)}, manifest_only={len(absent_from_object)}"
        )
    cistopic.add_cell_data(
        meta.loc[
            cistopic.cell_names,
            ["sample_id", "condition", "cell_type", "peak_group"],
        ],
        split_pattern=SEP,
    )

    matrix = cistopic.fragment_matrix
    cell_data = cistopic.cell_data.loc[cistopic.cell_names]
    support: dict[str, np.ndarray] = {}
    for condition in ("HC", "SSc"):
        donor_counts = np.zeros(matrix.shape[0], dtype=np.int16)
        donors = sorted(
            cell_data.loc[cell_data["condition"].astype(str) == condition, "sample_id"]
            .astype(str)
            .unique()
        )
        for donor in donors:
            cols = np.flatnonzero(cell_data["sample_id"].astype(str).to_numpy() == donor)
            donor_counts += (
                np.asarray((matrix[:, cols] > 0).sum(axis=1)).ravel() > 0
            ).astype(np.int16)
        support[condition] = donor_counts

    keep = np.maximum(support["HC"], support["SSc"]) >= args.min_donors_within_condition
    support_frame = pd.DataFrame(
        {
            "region": cistopic.region_names,
            "n_HC_donors_accessible": support["HC"],
            "n_SSc_donors_accessible": support["SSc"],
            "retained": keep,
        }
    )
    support_frame.to_csv(out / "region_donor_support.tsv.gz", sep="\t", index=False)
    retained = [region for region, flag in zip(cistopic.region_names, keep) if flag]
    if not retained:
        raise RuntimeError("Donor-support rule removed every condition-union region")
    cistopic = cistopic.subset(cells=None, regions=retained, copy=True)
    temporary = output_path.with_suffix(output_path.suffix + ".part")
    with temporary.open("wb") as handle:
        pickle.dump(cistopic, handle)
    os.replace(temporary, output_path)
    print(
        f"[cistopic] donor support retained {len(retained)}/{len(keep)} regions "
        f"-> {output_path}"
    )
    return output_path


def stage_dar_input(
    args: argparse.Namespace,
    full_meta: pd.DataFrame,
    cistopic_path: Path,
    consensus_path: Path,
) -> None:
    out = args.outdir / "fibroblast_donor_pseudobulk"
    matrix_path = out / "fibroblast_peak_by_donor.mtx.gz"
    regions_path = out / "regions.tsv.gz"
    samples_path = out / "samples.tsv"
    cells_path = out / "full_fibroblast_cell_manifest.tsv.gz"
    full_support_path = out / "full_fibroblast_region_donor_support.tsv.gz"
    manifest_path = out / "dar_input_manifest.json"
    complete_outputs = all(
        path.is_file() and path.stat().st_size > 0
        for path in (
            matrix_path,
            regions_path,
            samples_path,
            cells_path,
            full_support_path,
            manifest_path,
        )
    )
    if complete_outputs and not args.force:
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            existing.get("dar_cell_rule") == FROZEN["dar_cell_rule"]
            and existing.get("dar_region_rule") == FROZEN["dar_region_rule"]
        ):
            print(f"[dar-input] full-fibroblast outputs exist, keeping: {matrix_path}")
            return
        raise RuntimeError(
            "Existing DAR input was not built from every fibroblast barcode. "
            "Use the v2 output root or --force; do not reuse the sampled matrix."
        )
    if matrix_path.exists() and not complete_outputs and not args.force:
        print("[dar-input] incomplete checkpoint detected; regenerating all outputs")

    out.mkdir(parents=True, exist_ok=True)
    require_file(cistopic_path)

    fibroblast = full_meta["cell_type"].astype(str) == args.fibroblast_label
    if not fibroblast.any():
        raise ValueError(
            f"No cells match fibroblast label {args.fibroblast_label!r}"
        )
    fibro_meta = full_meta.loc[fibroblast].copy()
    donors = sorted(fibro_meta["sample_id"].astype(str).unique())
    donor_conditions = (
        fibro_meta[["sample_id", "condition"]]
        .drop_duplicates()
        .set_index("sample_id")["condition"]
    )
    if donor_conditions.index.duplicated().any():
        raise RuntimeError("A donor maps to more than one condition")
    condition_counts = donor_conditions.value_counts()
    if any(condition_counts.get(condition, 0) < 2 for condition in ("HC", "SSc")):
        raise RuntimeError("At least two fibroblast donors per condition are required")

    fragments = discover_fragments(args.fragments_dir)
    missing_fragments = sorted(set(donors) - set(fragments))
    if missing_fragments:
        raise FileNotFoundError(
            f"Missing fragment files for fibroblast donors: {missing_fragments}"
        )

    sample_dir = out / "per_sample_full_fibroblast_cistopic"
    tasks: list[dict[str, object]] = []
    expected_raw_barcodes: dict[str, set[str]] = {}
    for donor in donors:
        donor_cells = fibro_meta.index[
            fibro_meta["sample_id"].astype(str).to_numpy() == donor
        ]
        raw_barcodes = [str(name).split(SEP)[0] for name in donor_cells]
        if len(raw_barcodes) != len(set(raw_barcodes)):
            raise RuntimeError(
                f"Fibroblast raw barcodes are duplicated within donor {donor}"
            )
        expected_raw_barcodes[donor] = set(raw_barcodes)
        tasks.append(
            {
                "sample_id": donor,
                "fragment_path": fragments[donor],
                "consensus_path": str(consensus_path),
                "blacklist": str(args.blacklist),
                "raw_barcodes": raw_barcodes,
                "n_cpu": args.sample_cpu,
                "output_path": str(sample_dir / f"{safe_group(donor)}.pkl"),
                "force": args.force,
            }
        )

    completed: dict[str, tuple[str, int]] = {}
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.sample_workers
    ) as executor:
        future_to_sample = {
            executor.submit(_build_sample_cistopic, task): str(task["sample_id"])
            for task in tasks
        }
        for future in concurrent.futures.as_completed(future_to_sample):
            donor = future_to_sample[future]
            path, n_cells = future.result()
            completed[donor] = (path, n_cells)
            print(
                f"[dar-input] full fibroblast {donor}: {n_cells} cells -> {path}",
                flush=True,
            )

    consensus_bed = pd.read_table(
        consensus_path,
        header=None,
        usecols=[0, 1, 2],
        names=["chromosome", "start", "end"],
    )
    if (
        (consensus_bed["start"] < 0).any()
        or (consensus_bed["end"] <= consensus_bed["start"]).any()
    ):
        raise RuntimeError(f"Invalid interval in consensus BED: {consensus_path}")
    common_regions = pd.Index(
        consensus_bed.apply(
            lambda row: (
                f"{row['chromosome']}:{int(row['start'])}-{int(row['end'])}"
            ),
            axis=1,
        ),
        dtype=str,
    )
    if common_regions.empty or not common_regions.is_unique:
        raise RuntimeError("Condition-union consensus regions are empty or duplicated")

    donor_region_counts: dict[str, np.ndarray] = {}
    sample_rows = []
    support = {
        "HC": np.zeros(len(common_regions), dtype=np.int16),
        "SSc": np.zeros(len(common_regions), dtype=np.int16),
    }
    for donor in donors:
        with Path(completed[donor][0]).open("rb") as handle:
            donor_object = pickle.load(handle)
        observed_raw_barcodes = {
            str(name).split(SEP)[0] for name in donor_object.cell_names
        }
        missing_cells = expected_raw_barcodes[donor] - observed_raw_barcodes
        extra_cells = observed_raw_barcodes - expected_raw_barcodes[donor]
        if missing_cells or extra_cells:
            raise RuntimeError(
                f"Full-fibroblast barcode mismatch for {donor}: "
                f"missing={len(missing_cells)}, extra={len(extra_cells)}"
            )
        donor_regions = pd.Index(donor_object.region_names, dtype=str)
        if donor_regions.empty or not donor_regions.is_unique:
            raise RuntimeError(
                f"Donor regions are empty or duplicated for donor {donor}"
            )
        extra_regions = donor_regions.difference(common_regions)
        if len(extra_regions):
            raise RuntimeError(
                f"{len(extra_regions)} donor regions are absent from the "
                f"condition-union BED for {donor}; "
                f"examples: {extra_regions[:5].tolist()}"
            )
        positions = donor_regions.get_indexer(common_regions)
        donor_counts = np.zeros(len(common_regions), dtype=np.int64)
        present = np.flatnonzero(positions >= 0)
        donor_counts[present] = np.asarray(
            donor_object.fragment_matrix[positions[present], :].sum(axis=1)
        ).ravel()
        donor_region_counts[donor] = donor_counts
        condition = str(donor_conditions.loc[donor])
        support[condition] += (donor_counts > 0).astype(np.int16)
        sample_rows.append(
            {
                "sample_id": donor,
                "condition": condition,
                "n_fibroblast_cells": int(len(observed_raw_barcodes)),
            }
        )
    keep = (
        np.maximum(support["HC"], support["SSc"])
        >= args.min_donors_within_condition
    )
    if not keep.any():
        raise RuntimeError(
            "Full-fibroblast donor-support rule removed every condition-union peak"
        )
    retained_regions = common_regions[keep]
    pd.DataFrame(
        {
            "region": common_regions,
            "n_HC_fibroblast_donors_accessible": support["HC"],
            "n_SSc_fibroblast_donors_accessible": support["SSc"],
            "retained": keep,
        }
    ).to_csv(
        full_support_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    columns = [
        sparse.csc_matrix(donor_region_counts[donor][keep, None])
        for donor in donors
    ]
    pseudobulk = sparse.hstack(columns, format="csc")
    if np.any(pseudobulk.data < 0):
        raise ValueError("Negative fragment counts in donor pseudobulk matrix")
    matrix_tmp = matrix_path.with_suffix(matrix_path.suffix + ".part")
    regions_tmp = regions_path.with_suffix(regions_path.suffix + ".part")
    samples_tmp = samples_path.with_suffix(samples_path.suffix + ".part")
    with gzip.open(matrix_tmp, "wb") as handle:
        mmwrite(handle, pseudobulk, symmetry="general")
    pd.Series(retained_regions).to_csv(
        regions_tmp,
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
    pd.DataFrame(sample_rows).to_csv(
        samples_tmp, sep="\t", index=False
    )
    os.replace(matrix_tmp, matrix_path)
    os.replace(regions_tmp, regions_path)
    os.replace(samples_tmp, samples_path)
    fibro_meta[
        ["sample_id", "condition", "cell_type", "peak_group"]
    ].to_csv(cells_path, sep="\t", compression="gzip")
    manifest = {
        "design_version": FROZEN["design_version"],
        "dar_cell_rule": FROZEN["dar_cell_rule"],
        "dar_region_rule": FROZEN["dar_region_rule"],
        "n_full_fibroblast_cells": int(len(fibro_meta)),
        "n_donors": int(len(donors)),
        "n_HC_donors": int(condition_counts.get("HC", 0)),
        "n_SSc_donors": int(condition_counts.get("SSc", 0)),
        "n_condition_union_regions": int(len(common_regions)),
        "n_full_fibroblast_donor_supported_regions": int(len(retained_regions)),
        "per_donor_cell_counts": {
            row["sample_id"]: row["n_fibroblast_cells"] for row in sample_rows
        },
        "consensus_regions": str(consensus_path.resolve()),
        "efficiency_cistopic_object": str(cistopic_path.resolve()),
        "parameter_selection": "outcome_blind",
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(
        f"[dar-input] all {len(fibro_meta)} fibroblasts: "
        f"{pseudobulk.shape[0]} regions x {pseudobulk.shape[1]} donors -> "
        f"{matrix_path}"
    )


def main() -> None:
    args = parse_args()
    if args.target_cells < 100:
        raise ValueError("target-cells must be at least 100")
    if args.min_donors_within_condition < 2:
        raise ValueError("At least two donors within a condition are required")
    if args.n_cpu < 1 or args.sample_workers < 1 or args.sample_cpu < 1:
        raise ValueError("CPU and sample-worker settings must be positive")
    if args.sample_workers * args.sample_cpu > args.n_cpu:
        raise ValueError(
            "sample-workers * sample-cpu exceeds the declared CPU budget"
        )
    meta = load_metadata(args.adata)
    write_manifest(args, meta)
    resource_path = args.outdir / "manifest" / "execution_resources.json"
    resource_path.write_text(
        json.dumps(
            {
                "n_cpu": args.n_cpu,
                "sample_workers": args.sample_workers,
                "sample_cpu": args.sample_cpu,
                "total_sample_cpu_request": args.sample_workers * args.sample_cpu,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    if args.stage == "manifest":
        return

    consensus = args.outdir / "condition_peak_calling" / (
        "consensus_celltype_by_condition.bed"
    )
    if args.stage in {"peaks", "all"}:
        consensus = stage_peaks(args, meta)
        if args.stage == "peaks":
            return
    else:
        require_file(consensus)

    cistopic_path = args.outdir / "cistopic" / (
        "cistopic_condition_union_donor_supported.pkl"
    )
    if args.stage in {"cistopic", "all"}:
        selected_meta = load_stratified_meta(args, meta)
        cistopic_path = stage_cistopic(args, selected_meta, consensus)
        if args.stage == "cistopic":
            return
    else:
        require_file(cistopic_path)

    if args.stage in {"dar-input", "all"}:
        stage_dar_input(args, meta, cistopic_path, consensus)


if __name__ == "__main__":
    main()
