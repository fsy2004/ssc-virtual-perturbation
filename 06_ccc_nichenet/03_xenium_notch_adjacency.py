#!/usr/bin/env python3
"""Cell-boundary adjacency test for perivascular Notch signalling in Xenium.

Primary analysis
----------------
* Ten GSE312932 SSc Xenium sections.
* Cell types from RCTD, maximum weight >= 0.50.
* JAG1 senders: Endothelial or Pericyte cells with >=1 JAG1 transcript.
* DLL4 senders: Endothelial cells with >=1 DLL4 transcript.
* Receivers: Fib_proFibrotic cells with >=1 NOTCH2 or NOTCH3 transcript.
* Two cells are adjacent in the primary graph when their deposited Xenium
  cell-segmentation masks share an edge. The actual masks, not the
  visualization-only exported polygon approximations, define the main graph.
* The spatial graph is fixed. Expression-positive labels are permuted within
  cell-type x 500-micrometre spatial blocks, preserving the observed number of
  positive cells in every stratum and therefore the coarse block-level
  prevalence, not fine-scale expression autocorrelation.
* Empirical one-sided P values compare observed positive-positive edges with
  the permutation distribution. Four ligand-receptor pairs are corrected by
  Benjamini-Hochberg, Holm and a single-step studentized maxT procedure. Exact
  binomial Monte Carlo intervals quantify permutation uncertainty. Donor-level
  bootstrap intervals (one section per SSc donor) describe uncertainty in the
  pooled log2(observed/expected) enrichment; leave-one-donor-out estimates are
  also exported.

Sensitivity analyses vary the mask-gap tolerance, corner-touch inclusion,
annotation confidence, spatial-block size and graph definition. Approximate
exported polygons and a 20-micrometre centroid Delaunay graph are retained only
as graph-definition sensitivities. Grid-origin shifts test the sensitivity of
the study-specific spatial-block restriction.

Methodological precedents for fixed-graph label permutation include histoCAT
(doi:10.1038/nmeth.4391), Giotto (doi:10.1186/s13059-021-02286-2) and Squidpy
(doi:10.1038/s41592-021-01358-2). Context-restricted permutation is supported
by the comparative audit in doi:10.1038/s41467-026-71699-z. The 500-micrometre
restriction remains a study-specific design choice rather than a named or
universal null model; its outcome-independent parameter audit is implemented
in ``03a_xenium_parameter_audit.py``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import sys
from importlib.metadata import version as package_version
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import shapely
import zarr
from PIL import Image
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import beta as beta_distribution
from shapely.geometry import Polygon
from shapely.strtree import STRtree


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = (
    ROOT / ".codex_tmp" / "xenium_raw" / "raw" / "spatial" / "GSE312932"
    / "xenium_all.h5ad"
)
DEFAULT_BOUNDARIES = ROOT / ".codex_tmp" / "xenium_boundaries"
DEFAULT_MASKS = ROOT / ".codex_tmp" / "xenium_masks"
DEFAULT_RCTD = ROOT / "server_archive" / "tier1" / "powered" / "spatial" / "rctd"
DEFAULT_SPOTLIGHT = (
    ROOT / "server_archive" / "tier1" / "powered" / "spatial" / "spotlight"
)
DEFAULT_OUT = (
    ROOT / "04_manuscript" / "revision_20260722" / "figure6_spatial" / "outputs"
)
DEFAULT_PLOT_DATA = (
    ROOT / "04_manuscript" / "plot_data_local" / "ccc_nichenet"
    / "xenium_adjacency"
)

GENES = ("JAG1", "DLL4", "NOTCH2", "NOTCH3")
PAIR_SPECS = (
    ("JAG1", "NOTCH2"),
    ("JAG1", "NOTCH3"),
    ("DLL4", "NOTCH2"),
    ("DLL4", "NOTCH3"),
)
PAIR_LABELS = {
    ("JAG1", "NOTCH2"): "JAG1–NOTCH2",
    ("JAG1", "NOTCH3"): "JAG1–NOTCH3",
    ("DLL4", "NOTCH2"): "DLL4–NOTCH2",
    ("DLL4", "NOTCH3"): "DLL4–NOTCH3",
}
SENDER_TYPES = {
    "JAG1": ("Endothelial", "Pericyte"),
    "DLL4": ("Endothelial",),
}
RECEIVER_TYPE = "Fib_proFibrotic"


@dataclass(frozen=True)
class AnalysisConfig:
    name: str
    graph: str
    annotation: str = "rctd"
    rctd_threshold: float = 0.50
    spotlight_threshold: float = 0.50
    block_size_um: float | None = 500.0
    block_offset_xy_um: tuple[float, float] = (0.0, 0.0)
    primary: bool = False


@dataclass
class SectionData:
    sample: str
    cell_ids: np.ndarray
    xy: np.ndarray
    expression: dict[str, np.ndarray]
    rctd_label: np.ndarray
    rctd_max: np.ndarray
    spotlight_label: np.ndarray
    spotlight_max: np.ndarray
    graphs: dict[str, np.ndarray]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--boundaries", type=Path, default=DEFAULT_BOUNDARIES)
    parser.add_argument("--masks", type=Path, default=DEFAULT_MASKS)
    parser.add_argument("--rctd", type=Path, default=DEFAULT_RCTD)
    parser.add_argument("--spotlight", type=Path, default=DEFAULT_SPOTLIGHT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--plot-data-dir", type=Path, default=DEFAULT_PLOT_DATA)
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--sensitivity-permutations", type=int, default=5000)
    parser.add_argument("--whole-section-permutations", type=int, default=50000)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260722)
    parser.add_argument(
        "--max-sections",
        type=int,
        default=0,
        help="Limit sections for smoke tests only; 0 runs all sections.",
    )
    parser.add_argument(
        "--refresh-primary-statistics",
        action="store_true",
        help=(
            "Recompute primary multiplicity and Monte Carlo statistics from "
            "the saved joint null without rebuilding mask graphs."
        ),
    )
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return out


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = np.maximum.accumulate(ranked * (len(ranked) - np.arange(len(ranked))))
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return out


def monte_carlo_interval(
    exceedances: int,
    n_permutations: int,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Clopper-Pearson interval for the permutation-tail probability."""
    alpha = 1.0 - confidence
    lower = (
        0.0
        if exceedances == 0
        else float(
            beta_distribution.ppf(
                alpha / 2.0,
                exceedances,
                n_permutations - exceedances + 1,
            )
        )
    )
    upper = (
        1.0
        if exceedances == n_permutations
        else float(
            beta_distribution.ppf(
                1.0 - alpha / 2.0,
                exceedances + 1,
                n_permutations - exceedances,
            )
        )
    )
    return lower, upper


def read_boundary_polygons(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a gzipped Parquet boundary file in its deposited vertex order."""
    with gzip.open(path, "rb") as handle:
        frame = pd.read_parquet(io.BytesIO(handle.read()))
    required = {"cell_id", "vertex_x", "vertex_y"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    cell_id = frame["cell_id"].astype(str).to_numpy()
    coords = frame[["vertex_x", "vertex_y"]].to_numpy(dtype=float)
    starts = np.r_[0, np.flatnonzero(cell_id[1:] != cell_id[:-1]) + 1]
    ends = np.r_[starts[1:], len(frame)]
    ids = cell_id[starts]
    polygons: list[object] = []
    for start, end in zip(starts, ends):
        geom = Polygon(coords[start:end])
        if not geom.is_valid:
            geom = geom.buffer(0)
        if geom.is_empty:
            raise ValueError(f"Empty polygon in {path}: {cell_id[start]}")
        if geom.geom_type == "MultiPolygon":
            geom = max(geom.geoms, key=lambda item: item.area)
        polygons.append(geom)
    return ids, np.asarray(polygons, dtype=object)


def build_boundary_graphs(
    sample: str,
    local_cell_ids: np.ndarray,
    boundary_path: Path,
) -> dict[str, np.ndarray]:
    raw_ids, polygons = read_boundary_polygons(boundary_path)
    full_ids = np.asarray([f"{sample}_{item}" for item in raw_ids], dtype=object)
    local_lookup = {cell_id: pos for pos, cell_id in enumerate(local_cell_ids)}
    try:
        local_positions = np.asarray([local_lookup[item] for item in full_ids], dtype=int)
    except KeyError as exc:
        raise ValueError(
            f"Boundary cell {exc.args[0]} is absent from the H5AD sample {sample}"
        ) from exc
    if len(np.unique(local_positions)) != len(local_cell_ids):
        raise ValueError(
            f"Boundary/H5AD mismatch for {sample}: {len(local_positions)} polygons, "
            f"{len(local_cell_ids)} H5AD cells"
        )

    tree = STRtree(polygons)
    pairs = tree.query(polygons, predicate="dwithin", distance=2.0)
    left, right = pairs
    keep = left < right
    left = left[keep]
    right = right[keep]
    distances = shapely.distance(polygons[left], polygons[right])
    mapped = np.column_stack((local_positions[left], local_positions[right]))
    graphs = {}
    for tolerance in (0.5, 1.0, 2.0):
        key = f"boundary_{tolerance:g}um"
        graphs[key] = mapped[distances <= tolerance + 1e-9]
    return graphs


def build_mask_graphs(
    sample: str,
    xy: np.ndarray,
    rctd_label: np.ndarray,
    rctd_max: np.ndarray,
    mask_path: Path,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Build sender--receiver graphs from the deposited pixel-level cell mask.

    Xenium ``masks/1`` stores cell labels as one-based row indices.  We map
    those rows to H5AD cells by a one-to-one nearest match of the deposited
    cell centroids, then calculate the minimum Euclidean distance between the
    axis-aligned pixel squares of every perivascular-sender/profibrotic-receiver
    mask pair.  Only boundary pixels of the prespecified cell types enter the
    distance calculation.
    """
    store = zarr.storage.ZipStore(mask_path, mode="r")
    try:
        root = zarr.open_group(store=store, mode="r")
        attrs = dict(root.attrs)
        if attrs.get("name") != "CellSegmentationDataset":
            raise ValueError(f"Unexpected Zarr dataset in {mask_path}: {attrs}")
        centroids = np.asarray(root["cell_summary"][:, :2], dtype=float)
        if len(centroids) != len(xy):
            raise ValueError(
                f"Mask/H5AD cell-count mismatch for {sample}: "
                f"{len(centroids)} versus {len(xy)}"
            )
        centroid_distance, local_position = cKDTree(xy).query(centroids, k=1)
        if len(np.unique(local_position)) != len(xy) or centroid_distance.max() > 1e-3:
            raise ValueError(
                f"Mask/H5AD centroid mapping failed for {sample}: "
                f"max distance={centroid_distance.max():.6g} um; "
                f"unique matches={len(np.unique(local_position))}/{len(xy)}"
            )

        transform = np.asarray(root["masks/homogeneous_transform"])
        pixels_per_um_x = float(transform[0, 0])
        pixels_per_um_y = float(transform[1, 1])
        if not np.isclose(pixels_per_um_x, pixels_per_um_y, rtol=1e-6):
            raise ValueError(f"Anisotropic Xenium mask pixels in {mask_path}")
        pixel_size_um = 1.0 / pixels_per_um_x
        mask = root["masks/1"]
        n_cells = len(xy)
        mask_to_local = np.full(n_cells + 1, -1, dtype=np.int32)
        mask_to_local[1:] = local_position.astype(np.int32)

        labels_in_mask_order = rctd_label[local_position]
        confidence_in_mask_order = rctd_max[local_position]
        category = np.zeros(n_cells + 1, dtype=np.uint8)
        sender = (
            (confidence_in_mask_order >= 0.50)
            & np.isin(labels_in_mask_order, ("Endothelial", "Pericyte"))
        )
        receiver = (
            (confidence_in_mask_order >= 0.50)
            & (labels_in_mask_order == RECEIVER_TYPE)
        )
        category[1:][sender] = 1
        category[1:][receiver] = 2

        sender_points: list[np.ndarray] = []
        sender_mask_labels: list[np.ndarray] = []
        receiver_points: list[np.ndarray] = []
        receiver_mask_labels: list[np.ndarray] = []
        chunk_y, chunk_x = mask.chunks
        height, width = mask.shape
        directions = (
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1), (0, 1),
            (1, -1), (1, 0), (1, 1),
        )
        for y0 in range(0, height, chunk_y):
            y1 = min(y0 + chunk_y, height)
            ylo = max(0, y0 - 1)
            yhi = min(height, y1 + 1)
            for x0 in range(0, width, chunk_x):
                x1 = min(x0 + chunk_x, width)
                xlo = max(0, x0 - 1)
                xhi = min(width, x1 + 1)
                tile = np.asarray(mask[ylo:yhi, xlo:xhi])
                cy0, cy1 = y0 - ylo, y1 - ylo
                cx0, cx1 = x0 - xlo, x1 - xlo
                padded = np.pad(tile, pad_width=1, mode="constant", constant_values=0)
                cy0p, cy1p = cy0 + 1, cy1 + 1
                cx0p, cx1p = cx0 + 1, cx1 + 1
                central = padded[cy0p:cy1p, cx0p:cx1p]
                relevant = category[central] > 0
                if not relevant.any():
                    continue
                boundary = np.zeros_like(relevant)
                for dy, dx in directions:
                    neighbor = padded[
                        cy0p + dy:cy1p + dy,
                        cx0p + dx:cx1p + dx,
                    ]
                    boundary |= relevant & (neighbor != central)
                yy, xx = np.nonzero(boundary)
                values = central[yy, xx]
                classes = category[values]
                global_points = np.column_stack((yy + y0, xx + x0)).astype(
                    np.float32
                )
                if (classes == 1).any():
                    keep = classes == 1
                    sender_points.append(global_points[keep])
                    sender_mask_labels.append(values[keep])
                if (classes == 2).any():
                    keep = classes == 2
                    receiver_points.append(global_points[keep])
                    receiver_mask_labels.append(values[keep])

        if not sender_points or not receiver_points:
            empty = np.empty((0, 2), dtype=int)
            return (
                {
                    key: empty
                    for key in (
                        "mask_shared_edge",
                        "mask_corner_touch",
                        "mask_0.5um",
                        "mask_1um",
                        "mask_2um",
                    )
                },
                {
                    "mask_pixel_size_um": pixel_size_um,
                    "mask_max_centroid_mapping_error_um": float(centroid_distance.max()),
                    "n_mask_sender_boundary_pixels": 0,
                    "n_mask_receiver_boundary_pixels": 0,
                },
            )

        sender_xy = np.concatenate(sender_points)
        receiver_xy = np.concatenate(receiver_points)
        sender_values = np.concatenate(sender_mask_labels)
        receiver_values = np.concatenate(receiver_mask_labels)
        max_center_distance_pixels = 2.0 / pixel_size_um + np.sqrt(2.0)
        nearby = cKDTree(sender_xy).sparse_distance_matrix(
            cKDTree(receiver_xy),
            max_distance=max_center_distance_pixels,
            output_type="coo_matrix",
        )
        sxy = sender_xy[nearby.row]
        rxy = receiver_xy[nearby.col]
        delta = np.abs(sxy - rxy)
        shared_edge_touch = (
            ((delta[:, 0] == 1.0) & (delta[:, 1] == 0.0))
            | ((delta[:, 0] == 0.0) & (delta[:, 1] == 1.0))
        )
        pixel_square_gap = np.hypot(
            np.maximum(delta[:, 0] - 1.0, 0.0),
            np.maximum(delta[:, 1] - 1.0, 0.0),
        ) * pixel_size_um
        sender_local = mask_to_local[sender_values[nearby.row]]
        receiver_local = mask_to_local[receiver_values[nearby.col]]
        pair_code = sender_local.astype(np.int64) * n_cells + receiver_local
        order = np.argsort(pair_code, kind="stable")
        sorted_code = pair_code[order]
        starts = np.r_[0, np.flatnonzero(sorted_code[1:] != sorted_code[:-1]) + 1]
        unique_code = sorted_code[starts]
        minimum_gap = np.minimum.reduceat(pixel_square_gap[order], starts)
        has_shared_edge = np.maximum.reduceat(
            shared_edge_touch[order].astype(np.uint8), starts
        ).astype(bool)
        cell_pairs = np.column_stack(
            (unique_code // n_cells, unique_code % n_cells)
        ).astype(int)

        graphs = {
            "mask_shared_edge": cell_pairs[has_shared_edge],
            "mask_corner_touch": cell_pairs[minimum_gap <= 1e-9],
            "mask_0.5um": cell_pairs[minimum_gap <= 0.5 + 1e-9],
            "mask_1um": cell_pairs[minimum_gap <= 1.0 + 1e-9],
            "mask_2um": cell_pairs[minimum_gap <= 2.0 + 1e-9],
        }
        audit = {
            "mask_pixel_size_um": pixel_size_um,
            "mask_max_centroid_mapping_error_um": float(centroid_distance.max()),
            "n_mask_sender_boundary_pixels": len(sender_xy),
            "n_mask_receiver_boundary_pixels": len(receiver_xy),
            "n_mask_sender_receiver_pairs_within_2um": len(cell_pairs),
            "n_mask_shared_edge_pairs": len(graphs["mask_shared_edge"]),
            "n_mask_corner_touch_pairs": len(graphs["mask_corner_touch"]),
            "n_mask_pairs_within_0_5um": len(graphs["mask_0.5um"]),
            "n_mask_pairs_within_1um": len(graphs["mask_1um"]),
            "n_mask_pairs_within_2um": len(graphs["mask_2um"]),
        }
        return graphs, audit
    finally:
        store.close()


def build_delaunay_graph(xy: np.ndarray, max_distance: float = 20.0) -> np.ndarray:
    triangles = Delaunay(xy).simplices
    edges = np.vstack(
        (triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [0, 2]])
    )
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    distance = np.linalg.norm(xy[edges[:, 0]] - xy[edges[:, 1]], axis=1)
    return edges[distance <= max_distance]


def load_weight_labels(path: Path, cell_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    weights = pd.read_csv(path, index_col=0)
    maxima = weights.max(axis=1)
    labels = weights.idxmax(axis=1)
    return (
        labels.reindex(cell_ids).fillna("").to_numpy(dtype=object),
        maxima.reindex(cell_ids).fillna(0.0).to_numpy(dtype=float),
    )


def load_sections(args: argparse.Namespace) -> tuple[list[SectionData], pd.DataFrame]:
    adata = ad.read_h5ad(args.h5ad, backed="r")
    missing_genes = sorted(set(GENES).difference(adata.var_names))
    if missing_genes:
        raise ValueError(f"Genes absent from Xenium panel: {missing_genes}")

    obs_sample = adata.obs["sample"].astype(str).to_numpy()
    all_xy = np.asarray(adata.obsm["spatial"])
    sample_order = list(dict.fromkeys(obs_sample))
    if args.max_sections:
        sample_order = sample_order[: args.max_sections]
    sections: list[SectionData] = []
    concordance_rows: list[dict[str, object]] = []

    for sample in sample_order:
        sample_index = np.flatnonzero(obs_sample == sample)
        cell_ids = adata.obs_names[sample_index].astype(str).to_numpy()
        xy = all_xy[sample_index]
        expression = {}
        for gene in GENES:
            matrix = adata[sample_index, adata.var_names.get_loc(gene)].X
            values = matrix.toarray().ravel() if hasattr(matrix, "toarray") else np.asarray(matrix).ravel()
            expression[gene] = values > 0

        rctd_path = args.rctd / f"GSE312932_Xenium__{sample}_weights.csv"
        spotlight_path = args.spotlight / f"GSE312932_Xenium__{sample}_spotlight.csv"
        boundary_path = args.boundaries / f"{sample}_cell_boundaries.parquet.gz"
        gsm = sample.split("_", 1)[0]
        mask_path = args.masks / f"{gsm}_cells.zarr.zip"
        for required in (rctd_path, spotlight_path, boundary_path, mask_path):
            if not required.exists():
                raise FileNotFoundError(required)

        rctd_label, rctd_max = load_weight_labels(rctd_path, cell_ids)
        spotlight_label, spotlight_max = load_weight_labels(spotlight_path, cell_ids)
        both = (rctd_max >= 0.50) & (spotlight_max > 0)
        both_confident = both & (spotlight_max >= 0.50)
        concordance_rows.append(
            {
                "sample": sample,
                "n_cells": len(cell_ids),
                "n_rctd_any": int((rctd_max > 0).sum()),
                "n_rctd_ge_0_5": int((rctd_max >= 0.50).sum()),
                "n_spotlight_any": int((spotlight_max > 0).sum()),
                "n_both_rctd_ge_0_5": int(both.sum()),
                "agreement_rctd_ge_0_5": float(
                    (rctd_label[both] == spotlight_label[both]).mean()
                ) if both.any() else np.nan,
                "n_both_ge_0_5": int(both_confident.sum()),
                "agreement_both_ge_0_5": float(
                    (rctd_label[both_confident] == spotlight_label[both_confident]).mean()
                ) if both_confident.any() else np.nan,
            }
        )

        graphs, mask_audit = build_mask_graphs(
            sample, xy, rctd_label, rctd_max, mask_path
        )
        graphs.update(build_boundary_graphs(sample, cell_ids, boundary_path))
        graphs["delaunay_20um"] = build_delaunay_graph(xy, 20.0)
        concordance_rows[-1].update(mask_audit)
        sections.append(
            SectionData(
                sample=sample,
                cell_ids=cell_ids,
                xy=xy,
                expression=expression,
                rctd_label=rctd_label,
                rctd_max=rctd_max,
                spotlight_label=spotlight_label,
                spotlight_max=spotlight_max,
                graphs=graphs,
            )
        )

    adata.file.close()
    return sections, pd.DataFrame(concordance_rows)


def labels_for_config(section: SectionData, config: AnalysisConfig) -> np.ndarray:
    labels = np.full(len(section.cell_ids), "Unassigned", dtype=object)
    if config.annotation == "rctd":
        keep = section.rctd_max >= config.rctd_threshold
        labels[keep] = section.rctd_label[keep]
    elif config.annotation == "consensus":
        keep = (
            (section.rctd_max >= config.rctd_threshold)
            & (section.spotlight_max >= config.spotlight_threshold)
            & (section.rctd_label == section.spotlight_label)
        )
        labels[keep] = section.rctd_label[keep]
    else:
        raise ValueError(f"Unknown annotation mode: {config.annotation}")
    return labels


def spatial_blocks(
    xy: np.ndarray,
    block_size_um: float | None,
    offset_xy_um: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    if block_size_um is None:
        return np.zeros(len(xy), dtype=np.int64)
    shifted = (
        xy
        - xy.min(axis=0, keepdims=True)
        + np.asarray(offset_xy_um, dtype=float)[None, :]
    )
    bins = np.floor(shifted / block_size_um).astype(np.int64)
    return (bins[:, 0] << 32) + bins[:, 1]


def make_strata(
    labels: np.ndarray,
    blocks: np.ndarray,
    eligible_types: tuple[str, ...],
    positive: np.ndarray,
) -> list[tuple[np.ndarray, int]]:
    groups: list[tuple[np.ndarray, int]] = []
    for cell_type in eligible_types:
        type_index = np.flatnonzero(labels == cell_type)
        if not len(type_index):
            continue
        for block in np.unique(blocks[type_index]):
            group = type_index[blocks[type_index] == block]
            groups.append((group, int(positive[group].sum())))
    return groups


def permutation_matrix(
    groups: list[tuple[np.ndarray, int]],
    needed_cells: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
    chunk_size: int = 250,
) -> tuple[np.ndarray, np.ndarray]:
    """Return permuted positivity only for cells used by compatible graph edges."""
    needed = np.asarray(np.unique(needed_cells), dtype=int)
    output = np.zeros((n_permutations, len(needed)), dtype=bool)
    needed_column = {cell: col for col, cell in enumerate(needed)}
    for group, n_positive in groups:
        group_needed = np.asarray(
            [cell for cell in group if cell in needed_column], dtype=int
        )
        if not len(group_needed) or n_positive == 0:
            continue
        columns = np.asarray([needed_column[cell] for cell in group_needed], dtype=int)
        if n_positive == len(group):
            output[:, columns] = True
            continue
        local_position = {cell: pos for pos, cell in enumerate(group)}
        selected_positions = np.asarray(
            [local_position[cell] for cell in group_needed], dtype=int
        )
        for start in range(0, n_permutations, chunk_size):
            end = min(start + chunk_size, n_permutations)
            random_scores = rng.random((end - start, len(group)), dtype=np.float32)
            threshold = np.partition(
                random_scores, n_positive - 1, axis=1
            )[:, n_positive - 1]
            output[start:end, columns] = (
                random_scores[:, selected_positions] <= threshold[:, None]
            )
    return output, needed


def orient_compatible_edges(
    edges: np.ndarray,
    sender_mask: np.ndarray,
    receiver_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    forward = sender_mask[edges[:, 0]] & receiver_mask[edges[:, 1]]
    reverse = receiver_mask[edges[:, 0]] & sender_mask[edges[:, 1]]
    sender = np.r_[edges[forward, 0], edges[reverse, 1]]
    receiver = np.r_[edges[forward, 1], edges[reverse, 0]]
    return sender, receiver


def analyse_section(
    section: SectionData,
    config: AnalysisConfig,
    n_permutations: int,
    seed: int,
) -> tuple[list[dict[str, object]], np.ndarray]:
    labels = labels_for_config(section, config)
    blocks = spatial_blocks(
        section.xy, config.block_size_um, config.block_offset_xy_um
    )
    edges = section.graphs[config.graph]
    sender_masks = {
        gene: np.isin(labels, SENDER_TYPES[gene]) for gene in ("JAG1", "DLL4")
    }
    receiver_mask = labels == RECEIVER_TYPE

    oriented: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    needed: dict[str, list[np.ndarray]] = {gene: [] for gene in GENES}
    for ligand, receptor in PAIR_SPECS:
        sender, receiver = orient_compatible_edges(
            edges, sender_masks[ligand], receiver_mask
        )
        oriented[(ligand, receptor)] = (sender, receiver)
        needed[ligand].append(sender)
        needed[receptor].append(receiver)

    strata = {
        "JAG1": make_strata(
            labels, blocks, SENDER_TYPES["JAG1"], section.expression["JAG1"]
        ),
        "DLL4": make_strata(
            labels, blocks, SENDER_TYPES["DLL4"], section.expression["DLL4"]
        ),
        "NOTCH2": make_strata(
            labels, blocks, (RECEIVER_TYPE,), section.expression["NOTCH2"]
        ),
        "NOTCH3": make_strata(
            labels, blocks, (RECEIVER_TYPE,), section.expression["NOTCH3"]
        ),
    }

    rng = np.random.default_rng(seed)
    permuted: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for gene in GENES:
        cells = np.unique(np.concatenate(needed[gene])) if needed[gene] else np.array([], int)
        permuted[gene] = permutation_matrix(
            strata[gene], cells, n_permutations, rng
        )

    null_counts = np.zeros((n_permutations, len(PAIR_SPECS)), dtype=np.int32)
    rows: list[dict[str, object]] = []
    for pair_index, (ligand, receptor) in enumerate(PAIR_SPECS):
        sender, receiver = oriented[(ligand, receptor)]
        observed = int(
            (
                section.expression[ligand][sender]
                & section.expression[receptor][receiver]
            ).sum()
        )
        ligand_perm, ligand_cells = permuted[ligand]
        receptor_perm, receptor_cells = permuted[receptor]
        ligand_col = np.full(len(section.cell_ids), -1, dtype=int)
        receptor_col = np.full(len(section.cell_ids), -1, dtype=int)
        ligand_col[ligand_cells] = np.arange(len(ligand_cells))
        receptor_col[receptor_cells] = np.arange(len(receptor_cells))
        if len(sender):
            if (ligand_col[sender] < 0).any() or (receptor_col[receiver] < 0).any():
                raise RuntimeError("Permutation cell-index map is incomplete")
            null = (
                ligand_perm[:, ligand_col[sender]]
                & receptor_perm[:, receptor_col[receiver]]
            ).sum(axis=1)
        else:
            null = np.zeros(n_permutations, dtype=int)
        null_counts[:, pair_index] = null
        expected = float(null.mean())
        exceedances = int((null >= observed).sum())
        p_empirical = float((1 + exceedances) / (n_permutations + 1))
        mc_low, mc_high = monte_carlo_interval(exceedances, n_permutations)
        estimable = bool(
            len(sender)
            and (sender_masks[ligand] & section.expression[ligand]).any()
            and (receiver_mask & section.expression[receptor]).any()
        )
        rows.append(
            {
                "config": config.name,
                "sample": section.sample,
                "ligand": ligand,
                "receptor": receptor,
                "pair": PAIR_LABELS[(ligand, receptor)],
                "n_cells": len(section.cell_ids),
                "n_graph_edges": len(edges),
                "n_compatible_edges": len(sender),
                "n_sender_cells": int(sender_masks[ligand].sum()),
                "n_ligand_positive_senders": int(
                    (sender_masks[ligand] & section.expression[ligand]).sum()
                ),
                "n_receiver_cells": int(receiver_mask.sum()),
                "n_receptor_positive_receivers": int(
                    (receiver_mask & section.expression[receptor]).sum()
                ),
                "observed_positive_edges": observed,
                "expected_positive_edges": expected,
                "null_sd": float(null.std(ddof=1)),
                "log2_enrichment": float(np.log2((observed + 0.5) / (expected + 0.5)))
                if estimable else np.nan,
                "p_empirical_one_sided": p_empirical if estimable else np.nan,
                "p_mc_95ci_low": mc_low if estimable else np.nan,
                "p_mc_95ci_high": mc_high if estimable else np.nan,
                "estimable": estimable,
                "n_permutations": n_permutations,
            }
        )
    return rows, null_counts


def bootstrap_interval(
    observed: np.ndarray,
    expected: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.empty(n_bootstrap, dtype=float)
    n = len(observed)
    for start in range(0, n_bootstrap, 1000):
        end = min(start + 1000, n_bootstrap)
        index = rng.integers(0, n, size=(end - start, n))
        obs_sum = observed[index].sum(axis=1)
        exp_sum = expected[index].sum(axis=1)
        values[start:end] = np.log2((obs_sum + 0.5) / (exp_sum + 0.5))
    return tuple(np.quantile(values, [0.025, 0.975]))


def summarise_config(
    config: AnalysisConfig,
    per_sample: pd.DataFrame,
    null_by_sample: list[np.ndarray],
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    rows = []
    pooled_null = np.stack(null_by_sample).sum(axis=0)
    null_mean = pooled_null.mean(axis=0)
    null_sd = pooled_null.std(axis=0, ddof=1)
    if np.any(null_sd <= 0):
        failed = [PAIR_LABELS[PAIR_SPECS[i]] for i in np.flatnonzero(null_sd <= 0)]
        raise RuntimeError(f"maxT cannot studentize constant null counts: {failed}")
    standardized_null = (pooled_null - null_mean) / null_sd
    max_standardized_null = standardized_null.max(axis=1)
    for pair_index, (ligand, receptor) in enumerate(PAIR_SPECS):
        subset = per_sample[
            (per_sample["ligand"] == ligand)
            & (per_sample["receptor"] == receptor)
        ].copy()
        observed = subset["observed_positive_edges"].to_numpy(dtype=float)
        expected = subset["expected_positive_edges"].to_numpy(dtype=float)
        observed_total = int(observed.sum())
        expected_total = float(expected.sum())
        effect = float(np.log2((observed_total + 0.5) / (expected_total + 0.5)))
        lower, upper = bootstrap_interval(
            observed, expected, n_bootstrap, seed + pair_index
        )
        null = pooled_null[:, pair_index]
        exceedances = int((null >= observed_total).sum())
        p_empirical = float((1 + exceedances) / (len(null) + 1))
        mc_low, mc_high = monte_carlo_interval(exceedances, len(null))
        observed_studentized = float(
            (observed_total - null_mean[pair_index]) / null_sd[pair_index]
        )
        max_t_exceedances = int(
            (max_standardized_null >= observed_studentized).sum()
        )
        p_max_t = float((1 + max_t_exceedances) / (len(null) + 1))
        max_t_mc_low, max_t_mc_high = monte_carlo_interval(
            max_t_exceedances, len(null)
        )
        rows.append(
            {
                "config": config.name,
                "ligand": ligand,
                "receptor": receptor,
                "pair": PAIR_LABELS[(ligand, receptor)],
                "n_sections": len(subset),
                "n_estimable_sections": int(subset["estimable"].sum()),
                "observed_positive_edges": observed_total,
                "expected_positive_edges": expected_total,
                "log2_enrichment": effect,
                "bootstrap_95ci_low": lower,
                "bootstrap_95ci_high": upper,
                "p_empirical_one_sided": p_empirical,
                "p_mc_95ci_low": mc_low,
                "p_mc_95ci_high": mc_high,
                "studentized_observed_excess": observed_studentized,
                "p_maxT_single_step": p_max_t,
                "p_maxT_mc_95ci_low": max_t_mc_low,
                "p_maxT_mc_95ci_high": max_t_mc_high,
                "n_permutations": len(null),
                "n_bootstrap": n_bootstrap,
            }
        )
    result = pd.DataFrame(rows)
    result["q_bh_four_pairs"] = bh_adjust(result["p_empirical_one_sided"])
    result["p_holm_four_pairs"] = holm_adjust(result["p_empirical_one_sided"])
    return result


def leave_one_section_out_summary(
    config: AnalysisConfig,
    per_sample: pd.DataFrame,
) -> pd.DataFrame:
    """Export pooled effects after omitting each independent donor/section."""
    rows: list[dict[str, object]] = []
    for ligand, receptor in PAIR_SPECS:
        subset = per_sample[
            (per_sample["ligand"] == ligand)
            & (per_sample["receptor"] == receptor)
        ].copy()
        for omitted_sample in subset["sample"]:
            kept = subset[subset["sample"] != omitted_sample]
            observed = float(kept["observed_positive_edges"].sum())
            expected = float(kept["expected_positive_edges"].sum())
            rows.append(
                {
                    "config": config.name,
                    "ligand": ligand,
                    "receptor": receptor,
                    "pair": PAIR_LABELS[(ligand, receptor)],
                    "omitted_sample": omitted_sample,
                    "n_sections_retained": len(kept),
                    "observed_positive_edges": observed,
                    "expected_positive_edges": expected,
                    "log2_enrichment": float(
                        np.log2((observed + 0.5) / (expected + 0.5))
                    ),
                }
            )
    return pd.DataFrame(rows)


def run_config(
    sections: list[SectionData],
    config: AnalysisConfig,
    n_permutations: int,
    n_bootstrap: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    rows: list[dict[str, object]] = []
    null_by_sample: list[np.ndarray] = []
    for section_index, section in enumerate(sections):
        section_rows, null = analyse_section(
            section,
            config,
            n_permutations,
            seed + 10000 * section_index,
        )
        rows.extend(section_rows)
        null_by_sample.append(null)
    per_sample = pd.DataFrame(rows)
    summary = summarise_config(
        config,
        per_sample,
        null_by_sample,
        n_bootstrap,
        seed + 900000,
    )
    return per_sample, summary, np.stack(null_by_sample)


def set_plot_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    matplotlib.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def q_label(value: float) -> str:
    if value < 0.001:
        return "q<0.001"
    return f"q={value:.3f}"


def plot_primary(
    per_sample: pd.DataFrame,
    summary: pd.DataFrame,
    outdir: Path,
) -> None:
    set_plot_style()
    order = [PAIR_LABELS[pair] for pair in PAIR_SPECS]
    colors = {
        "JAG1–NOTCH2": "#0072B2",
        "JAG1–NOTCH3": "#56B4E9",
        "DLL4–NOTCH2": "#D55E00",
        "DLL4–NOTCH3": "#E69F00",
    }
    fig, ax = plt.subplots(figsize=(4.8, 3.0), constrained_layout=True)
    rng = np.random.default_rng(1402)
    for y, pair in enumerate(order):
        section = per_sample[(per_sample["pair"] == pair) & per_sample["estimable"]]
        jitter = rng.normal(0, 0.055, len(section))
        ax.scatter(
            section["log2_enrichment"],
            y + jitter,
            s=20,
            facecolor="white",
            edgecolor=colors[pair],
            linewidth=0.8,
            alpha=0.9,
            zorder=2,
        )
        pooled = summary[summary["pair"] == pair].iloc[0]
        x = float(pooled["log2_enrichment"])
        lo = float(pooled["bootstrap_95ci_low"])
        hi = float(pooled["bootstrap_95ci_high"])
        ax.errorbar(
            x,
            y,
            xerr=np.asarray([[x - lo], [hi - x]]),
            fmt="D",
            ms=5.5,
            color=colors[pair],
            ecolor=colors[pair],
            elinewidth=1.4,
            capsize=2.5,
            markeredgecolor="white",
            markeredgewidth=0.5,
            zorder=4,
        )
        ax.text(
            1.01,
            y,
            q_label(float(pooled["q_bh_four_pairs"])),
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=6.8,
            color="#333333",
            clip_on=False,
        )
    ax.axvline(0, color="#666666", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlim(-2.5, 2.5)
    ax.set_xlabel("Mask shared-edge adjacency enrichment, log$_2$(observed/expected)")
    ax.set_title("Xenium mask shared-edge adjacency (10 SSc sections)", pad=15)
    ax.text(
        0.5,
        1.015,
        "open circles: sections; diamonds/lines: pooled effect and section-bootstrap 95% CI",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.3,
        color="#444444",
    )
    ax.grid(axis="y", visible=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for extension in ("pdf", "png"):
        fig.savefig(
            outdir / f"Figure6e_Xenium_mask_shared_edge_adjacency.{extension}",
            dpi=600 if extension == "png" else None,
            bbox_inches="tight",
        )
    Image.open(outdir / "Figure6e_Xenium_mask_shared_edge_adjacency.png").convert("L").save(
        outdir / "Figure6e_Xenium_mask_shared_edge_adjacency_grayscale.png"
    )
    plt.close(fig)


def plot_sensitivity(summary: pd.DataFrame, outdir: Path) -> None:
    set_plot_style()
    config_order = list(dict.fromkeys(summary["config"]))
    pair_order = [PAIR_LABELS[pair] for pair in PAIR_SPECS]
    matrix = summary.pivot(index="config", columns="pair", values="log2_enrichment")
    matrix = matrix.reindex(index=config_order, columns=pair_order)
    q_matrix = summary.pivot(index="config", columns="pair", values="q_bh_four_pairs")
    q_matrix = q_matrix.reindex(index=config_order, columns=pair_order)
    annotations = matrix.copy().astype(object)
    for row in matrix.index:
        for col in matrix.columns:
            effect = float(matrix.loc[row, col])
            q_value = float(q_matrix.loc[row, col])
            annotations.loc[row, col] = f"{effect:.2f}" + ("*" if q_value < 0.05 else "")
    vmax = max(1.0, float(np.nanmax(np.abs(matrix.to_numpy(dtype=float)))))
    fig, ax = plt.subplots(figsize=(6.8, 4.1), constrained_layout=True)
    sns.heatmap(
        matrix,
        annot=annotations,
        fmt="",
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "log$_2$(observed/expected)"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Xenium adjacency sensitivity analyses (* BH q<0.05)")
    ax.tick_params(axis="x", rotation=0)
    ax.tick_params(axis="y", rotation=0)
    for extension in ("pdf", "png"):
        fig.savefig(
            outdir / f"Figure6e_Xenium_adjacency_sensitivity.{extension}",
            dpi=600 if extension == "png" else None,
            bbox_inches="tight",
        )
    Image.open(outdir / "Figure6e_Xenium_adjacency_sensitivity.png").convert("L").save(
        outdir / "Figure6e_Xenium_adjacency_sensitivity_grayscale.png"
    )
    plt.close(fig)


def primary_source_data(
    per_sample: pd.DataFrame,
    summary: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "pair",
        "log2_enrichment",
        "bootstrap_95ci_low",
        "bootstrap_95ci_high",
        "p_empirical_one_sided",
        "p_mc_95ci_low",
        "p_mc_95ci_high",
        "q_bh_four_pairs",
        "p_holm_four_pairs",
        "p_maxT_single_step",
        "p_maxT_mc_95ci_low",
        "p_maxT_mc_95ci_high",
    ]
    return per_sample.merge(
        summary[columns].rename(
            columns={
                "log2_enrichment": "pooled_log2_enrichment",
                "p_empirical_one_sided": "pooled_p_empirical_one_sided",
            }
        ),
        on="pair",
        how="left",
    )


def refresh_primary_statistics(args: argparse.Namespace) -> None:
    """Refresh primary inference from the retained joint permutation draws."""
    per_sample_path = args.outdir / "primary_per_section.csv"
    null_path = args.outdir / "primary_null_counts.npz"
    if not per_sample_path.exists() or not null_path.exists():
        raise FileNotFoundError(
            "Primary cache is incomplete; expected primary_per_section.csv and "
            "primary_null_counts.npz"
        )
    per_sample = pd.read_csv(per_sample_path)
    with np.load(null_path) as archive:
        null_counts = np.asarray(archive["null_counts"])
        samples = np.asarray(archive["samples"]).astype(str)
        pairs = np.asarray(archive["pairs"]).astype(str)
    expected_pairs = np.asarray([PAIR_LABELS[pair] for pair in PAIR_SPECS])
    if null_counts.ndim != 3 or null_counts.shape[2] != len(PAIR_SPECS):
        raise RuntimeError(f"Unexpected cached null shape: {null_counts.shape}")
    if not np.array_equal(pairs, expected_pairs):
        raise RuntimeError("Cached pair order does not match the prespecified family")
    observed_samples = per_sample["sample"].drop_duplicates().astype(str).to_numpy()
    if not np.array_equal(samples, observed_samples):
        raise RuntimeError("Cached section order does not match primary_per_section.csv")

    config = AnalysisConfig(
        "Primary: mask shared edge; RCTD >=0.5; block 500 um",
        "mask_shared_edge",
        primary=True,
    )
    summary = summarise_config(
        config,
        per_sample,
        [null_counts[index] for index in range(len(null_counts))],
        args.bootstrap,
        args.seed + 900000,
    )
    summary.to_csv(args.outdir / "primary_summary.csv", index=False)
    summary.to_csv(args.plot_data_dir / "primary_summary.csv", index=False)
    source_data = primary_source_data(per_sample, summary)
    source_data.to_csv(args.outdir / "Figure6e_source_data.csv", index=False)
    source_data.to_csv(args.plot_data_dir / "Figure6e_source_data.csv", index=False)
    leave_one_section_out_summary(config, per_sample).to_csv(
        args.outdir / "primary_leave_one_section_out.csv", index=False
    )
    plot_primary(per_sample, summary, args.outdir)

    manifest_path = args.outdir / "run_manifest.json"
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        manifest["statistics_refreshed_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["script_sha256"] = sha256(Path(__file__).resolve())
        manifest["primary_results"] = summary.to_dict(orient="records")
        manifest["primary_inference"] = {
            "unadjusted": "plus-one one-sided permutation p value",
            "monte_carlo_interval": "95% Clopper-Pearson interval for the permutation-tail probability",
            "multiple_testing_family": [str(value) for value in expected_pairs],
            "adjustments": [
                "Benjamini-Hochberg",
                "Holm",
                "single-step studentized maxT",
            ],
            "maxT_statistic": "pooled edge-count excess standardized by its permutation SD",
        }
        refreshed = [
            args.outdir / "primary_summary.csv",
            args.outdir / "primary_leave_one_section_out.csv",
            args.outdir / "Figure6e_source_data.csv",
            args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency.pdf",
            args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency.png",
            args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency_grayscale.png",
            args.plot_data_dir / "primary_summary.csv",
            args.plot_data_dir / "Figure6e_source_data.csv",
        ]
        manifest.setdefault("outputs", {}).update(
            {str(path): sha256(path) for path in refreshed}
        )
        with manifest_path.open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, ensure_ascii=False)
    print("Refreshed primary inference from cached joint permutation draws")
    print(summary.to_string(index=False))


def write_manifest(
    args: argparse.Namespace,
    configs: list[AnalysisConfig],
    outputs: list[Path],
    primary_summary: pd.DataFrame,
) -> None:
    boundary_files = sorted(args.boundaries.glob("*_cell_boundaries.parquet.gz"))
    mask_files = sorted(args.masks.glob("*_cells.zarr.zip"))
    weight_files = sorted(args.rctd.glob("GSE312932_Xenium__*_weights.csv"))
    spotlight_files = sorted(
        args.spotlight.glob("GSE312932_Xenium__*_spotlight.csv")
    )
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "command_parameters": vars(args) | {
            "h5ad": str(args.h5ad),
            "boundaries": str(args.boundaries),
            "masks": str(args.masks),
            "rctd": str(args.rctd),
            "spotlight": str(args.spotlight),
            "outdir": str(args.outdir),
            "plot_data_dir": str(args.plot_data_dir),
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "anndata": package_version("anndata"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "h5py": h5py.__version__,
            "shapely": shapely.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
        },
        "geo_accession": "GSE312932",
        "geo_boundary_url_template": (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9357nnn/{GSM}/suppl/"
            "{GSM}_{SSc_sample}_cell_boundaries.parquet.gz"
        ),
        "inputs": {
            str(args.h5ad): sha256(args.h5ad),
            **{str(path): sha256(path) for path in mask_files},
            **{str(path): sha256(path) for path in boundary_files},
            **{str(path): sha256(path) for path in weight_files},
            **{str(path): sha256(path) for path in spotlight_files},
        },
        "configs": [config.__dict__ for config in configs],
        "primary_results": primary_summary.to_dict(orient="records"),
        "outputs": {str(path): sha256(path) for path in outputs if path.exists()},
    }
    with (args.outdir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)


def main() -> None:
    args = parse_args()
    if args.permutations < 999 or args.sensitivity_permutations < 999:
        raise ValueError("At least 999 permutations are required")
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.plot_data_dir.mkdir(parents=True, exist_ok=True)
    if args.refresh_primary_statistics:
        refresh_primary_statistics(args)
        return

    configs = [
        AnalysisConfig(
            "Primary: mask shared edge; RCTD >=0.5; block 500 um",
            "mask_shared_edge",
            primary=True,
        ),
        AnalysisConfig(
            "Mask corner-touch contact", "mask_corner_touch"
        ),
        AnalysisConfig(
            "Mask <=0.5 um dilation", "mask_0.5um"
        ),
        AnalysisConfig(
            "Mask <=1 um dilation", "mask_1um"
        ),
        AnalysisConfig(
            "Mask <=2 um dilation", "mask_2um"
        ),
        AnalysisConfig(
            "Exported polygon <=1 um", "boundary_1um"
        ),
        AnalysisConfig(
            "Exported polygon <=0.5 um", "boundary_0.5um"
        ),
        AnalysisConfig(
            "Exported polygon <=2 um", "boundary_2um"
        ),
        AnalysisConfig(
            "Centroid Delaunay <=20 um", "delaunay_20um"
        ),
        AnalysisConfig(
            "RCTD >=0.7", "mask_shared_edge", rctd_threshold=0.70
        ),
        AnalysisConfig(
            "RCTD-SPOTlight consensus >=0.5",
            "mask_shared_edge",
            annotation="consensus",
        ),
        AnalysisConfig(
            "Spatial block 250 um", "mask_shared_edge", block_size_um=250.0
        ),
        AnalysisConfig(
            "Spatial block 1000 um", "mask_shared_edge", block_size_um=1000.0
        ),
        AnalysisConfig(
            "Spatial block 500 um; x-origin shifted 250 um",
            "mask_shared_edge",
            block_offset_xy_um=(250.0, 0.0),
        ),
        AnalysisConfig(
            "Spatial block 500 um; y-origin shifted 250 um",
            "mask_shared_edge",
            block_offset_xy_um=(0.0, 250.0),
        ),
        AnalysisConfig(
            "Spatial block 500 um; xy-origin shifted 250 um",
            "mask_shared_edge",
            block_offset_xy_um=(250.0, 250.0),
        ),
        AnalysisConfig(
            "Whole-section type-stratified null", "mask_shared_edge", block_size_um=None
        ),
    ]

    sections, concordance = load_sections(args)
    available_graphs = set().union(*(section.graphs.keys() for section in sections))
    configs = [config for config in configs if config.graph in available_graphs]
    if not any(config.primary for config in configs):
        raise RuntimeError("Primary mask shared-edge graph was not available")

    concordance.to_csv(args.outdir / "annotation_concordance_by_section.csv", index=False)

    all_per_sample = []
    all_summary = []
    primary_null = None
    primary_per_sample = None
    primary_summary = None
    for config_index, config in enumerate(configs):
        n_permutations = (
            args.permutations if config.primary else args.sensitivity_permutations
        )
        per_sample, summary, null = run_config(
            sections,
            config,
            n_permutations,
            args.bootstrap,
            args.seed + config_index * 1000000,
        )
        all_per_sample.append(per_sample)
        all_summary.append(summary)
        if config.primary:
            primary_per_sample = per_sample
            primary_summary = summary
            primary_null = null

    if primary_per_sample is None or primary_summary is None or primary_null is None:
        raise RuntimeError("Primary configuration was not run")
    per_sample_all = pd.concat(all_per_sample, ignore_index=True)
    summary_all = pd.concat(all_summary, ignore_index=True)

    primary_per_sample.to_csv(args.outdir / "primary_per_section.csv", index=False)
    primary_summary.to_csv(args.outdir / "primary_summary.csv", index=False)
    per_sample_all.to_csv(args.outdir / "sensitivity_per_section.csv", index=False)
    summary_all.to_csv(args.outdir / "sensitivity_summary.csv", index=False)
    np.savez_compressed(
        args.outdir / "primary_null_counts.npz",
        null_counts=primary_null,
        samples=np.asarray([section.sample for section in sections]),
        pairs=np.asarray([PAIR_LABELS[pair] for pair in PAIR_SPECS]),
    )

    source_data = primary_source_data(primary_per_sample, primary_summary)
    source_data.to_csv(args.outdir / "Figure6e_source_data.csv", index=False)
    source_data.to_csv(args.plot_data_dir / "Figure6e_source_data.csv", index=False)
    primary_summary.to_csv(args.plot_data_dir / "primary_summary.csv", index=False)
    primary_leave_one_out = leave_one_section_out_summary(
        configs[0], primary_per_sample
    )
    primary_leave_one_out.to_csv(
        args.outdir / "primary_leave_one_section_out.csv", index=False
    )

    plot_primary(primary_per_sample, primary_summary, args.outdir)
    plot_sensitivity(summary_all, args.outdir)

    output_files = [
        args.outdir / "annotation_concordance_by_section.csv",
        args.outdir / "primary_per_section.csv",
        args.outdir / "primary_summary.csv",
        args.outdir / "primary_leave_one_section_out.csv",
        args.outdir / "sensitivity_per_section.csv",
        args.outdir / "sensitivity_summary.csv",
        args.outdir / "primary_null_counts.npz",
        args.outdir / "Figure6e_source_data.csv",
        args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency.pdf",
        args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency.png",
        args.outdir / "Figure6e_Xenium_mask_shared_edge_adjacency_grayscale.png",
        args.outdir / "Figure6e_Xenium_adjacency_sensitivity.pdf",
        args.outdir / "Figure6e_Xenium_adjacency_sensitivity.png",
        args.outdir / "Figure6e_Xenium_adjacency_sensitivity_grayscale.png",
        args.plot_data_dir / "Figure6e_source_data.csv",
        args.plot_data_dir / "primary_summary.csv",
    ]
    write_manifest(args, configs, output_files, primary_summary)

    print("Primary Xenium mask shared-edge adjacency results")
    print(primary_summary.to_string(index=False))
    print(f"Outputs: {args.outdir}")


if __name__ == "__main__":
    main()
