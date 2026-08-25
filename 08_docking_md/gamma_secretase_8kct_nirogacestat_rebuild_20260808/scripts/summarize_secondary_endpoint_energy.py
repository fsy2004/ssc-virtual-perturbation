#!/usr/bin/env python3
"""Summarize frozen 300-frame endpoint energies by replica and fixed blocks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import statistics
from pathlib import Path
from typing import Any


MODELS = ("PB_membrane_indi4",)
REPLICAS = ("rep01", "rep02", "rep03")
FROZEN_DECOMP_RESIDUES = (
    "R:B:VAL:261",
    "R:B:LEU:268",
    "R:B:VAL:272",
    "R:B:LEU:282",
    "R:B:ILE:287",
    "R:B:LYS:380",
    "R:B:LEU:381",
    "R:B:ALA:431",
    "R:B:LEU:432",
    "L:B:O6U:502",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_delta_title(row: list[str]) -> bool:
    text = " ".join(cell.strip().lower() for cell in row)
    return all(word in text for word in ("delta", "complex", "receptor", "ligand"))


def parse_delta_energy_csv(text: str) -> list[dict[str, float | int]]:
    rows = list(csv.reader(io.StringIO(text)))
    in_delta = False
    header: list[str] | None = None
    start = 0
    parsed: list[dict[str, float | int]] = []
    for row in rows:
        if _is_delta_title(row):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        normalized = [cell.strip() for cell in row]
        if header is None:
            frame_positions = [index for index, value in enumerate(normalized) if value.lower() in ("frame", "frame #")]
            if not frame_positions:
                continue
            start = frame_positions[0]
            header = normalized[start:]
            if len(header) < 2 or len(header) != len(set(header)):
                raise ValueError("delta CSV header is invalid")
            continue
        if not any(normalized):
            if parsed:
                break
            continue
        values = normalized[start:start + len(header)]
        if len(values) != len(header):
            if parsed:
                break
            continue
        try:
            frame = int(float(values[0]))
            numeric = [float(value) for value in values[1:]]
        except ValueError:
            if parsed:
                break
            continue
        if any(not math.isfinite(value) for value in numeric):
            raise ValueError("delta CSV contains non-finite values")
        parsed.append({"frame": frame, **dict(zip(header[1:], numeric))})
    if not parsed:
        raise ValueError("delta binding-energy section was not found")
    return parsed


def fixed_block_means(rows: list[dict[str, float | int]]) -> list[dict[str, Any]]:
    if len(rows) != 300:
        raise ValueError("formal endpoint-energy series must contain exactly 300 frames")
    if [int(row.get("frame", -1)) for row in rows] != list(range(1, 301)):
        raise ValueError("formal endpoint-energy frame numbers must be 1..300 without gaps")
    terms = [key for key in rows[0] if key != "frame"]
    if not terms or any(set(row) != set(rows[0]) for row in rows):
        raise ValueError("formal endpoint-energy term columns are empty or inconsistent")
    for row in rows:
        if any(not math.isfinite(float(row[term])) for term in terms):
            raise ValueError("formal endpoint-energy series contains non-finite values")
    blocks = []
    for block_index in range(5):
        subset = rows[block_index * 60:(block_index + 1) * 60]
        blocks.append({
            "block_index_zero_based": block_index,
            "start_ns": 200.0 + block_index * 60.0,
            "end_ns": 260.0 + block_index * 60.0,
            "frame_count": 60,
            "means": {term: statistics.fmean(float(row[term]) for row in subset) for term in terms},
        })
    return blocks


def parse_decomposition_csv(text: str) -> list[dict[str, Any]]:
    """Parse only the DELTAS per-residue TDC table from gmx_MMPBSA -deo output."""
    rows = list(csv.reader(io.StringIO(text)))
    in_delta = False
    header: list[str] | None = None
    start = 0
    parsed: list[dict[str, Any]] = []
    frozen = set(FROZEN_DECOMP_RESIDUES)
    for row in rows:
        normalized = [cell.strip() for cell in row]
        joined = " ".join(normalized).strip().lower()
        if joined.startswith("deltas"):
            in_delta = True
            header = None
            parsed = []
            continue
        if not in_delta:
            continue
        if header is None:
            frame_positions = [
                index for index, value in enumerate(normalized)
                if value.lower() in ("frame", "frame #")
            ]
            if not frame_positions:
                continue
            start = frame_positions[0]
            header = normalized[start:]
            lowered = [value.lower() for value in header]
            if "resid 2" in lowered or "residue 2" in lowered:
                raise ValueError("decomposition output must be per-residue, not pairwise")
            if not any(value in lowered for value in ("residue", "resid 1")):
                raise ValueError("per-residue decomposition header is invalid")
            if "total" not in lowered or len(header) != len(set(header)):
                raise ValueError("per-residue decomposition terms are invalid")
            continue
        if not any(normalized):
            if parsed:
                break
            continue
        values = normalized[start:start + len(header)]
        if len(values) != len(header):
            if parsed:
                break
            continue
        record = dict(zip(header, values))
        residue_key = next(key for key in header if key.lower() in ("residue", "resid 1"))
        try:
            frame = int(float(record[header[0]]))
            residue = record[residue_key]
            numeric = {
                key: float(value)
                for key, value in record.items()
                if key not in (header[0], residue_key)
            }
        except ValueError:
            if parsed:
                break
            continue
        if residue not in frozen:
            raise ValueError(f"decomposition output contains non-frozen residue: {residue}")
        if not numeric or any(not math.isfinite(value) for value in numeric.values()):
            raise ValueError("decomposition output contains missing or non-finite terms")
        parsed.append({"frame": frame, "residue": residue, **numeric})
    if not parsed:
        raise ValueError("DELTAS per-residue decomposition table was not found")
    return parsed


def fixed_decomposition_summary(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Summarize the fixed residue list descriptively in its frozen order."""
    validate_fixed_decomposition_rows(rows, frame_count=300)
    grouped = {residue: [] for residue in FROZEN_DECOMP_RESIDUES}
    for row in rows:
        residue = str(row.get("residue", ""))
        if residue not in grouped:
            raise ValueError(f"decomposition output contains non-frozen residue: {residue}")
        grouped[residue].append({key: value for key, value in row.items() if key != "residue"})
    result: dict[str, dict[str, Any]] = {}
    for residue in FROZEN_DECOMP_RESIDUES:
        residue_rows = grouped[residue]
        blocks = fixed_block_means(residue_rows)
        terms = list(blocks[0]["means"])
        result[residue] = {
            "frame_count": 300,
            "window_means": {
                term: statistics.fmean(float(row[term]) for row in residue_rows)
                for term in terms
            },
            "blocks": blocks,
        }
    return result


def validate_fixed_decomposition_rows(
    rows: list[dict[str, Any]], frame_count: int
) -> dict[str, Any]:
    if frame_count <= 0:
        raise ValueError("decomposition frame count must be positive")
    grouped = {residue: [] for residue in FROZEN_DECOMP_RESIDUES}
    expected_columns: set[str] | None = None
    for row in rows:
        residue = str(row.get("residue", ""))
        if residue not in grouped:
            raise ValueError(f"decomposition output contains non-frozen residue: {residue}")
        columns = set(row)
        if expected_columns is None:
            expected_columns = columns
        elif columns != expected_columns:
            raise ValueError("decomposition term columns drifted between rows")
        grouped[residue].append(int(row.get("frame", -1)))
    expected_frames = list(range(1, frame_count + 1))
    if any(grouped[residue] != expected_frames for residue in FROZEN_DECOMP_RESIDUES):
        raise ValueError(
            "decomposition frames must be complete and ordered for every frozen residue"
        )
    return {
        "status": "pass",
        "frame_count": frame_count,
        "residues_in_frozen_order": list(FROZEN_DECOMP_RESIDUES),
    }


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def hierarchical_block_bootstrap(
    replica_blocks: dict[str, list[float]], seed: int, draws: int = 10000
) -> dict[str, float | int]:
    if set(replica_blocks) != set(REPLICAS) or any(len(values) != 5 for values in replica_blocks.values()):
        raise ValueError("hierarchical bootstrap requires three replicas and five fixed blocks each")
    if draws <= 0 or any(not math.isfinite(value) for values in replica_blocks.values() for value in values):
        raise ValueError("bootstrap inputs must be finite and draws positive")
    estimate = statistics.fmean(value for values in replica_blocks.values() for value in values)
    generator = random.Random(seed)
    samples = []
    for _ in range(draws):
        selected_replicas = [generator.choice(REPLICAS) for _ in REPLICAS]
        values = []
        for replica in selected_replicas:
            values.extend(generator.choices(replica_blocks[replica], k=5))
        samples.append(statistics.fmean(values))
    return {
        "estimate": estimate,
        "ci95_low": _quantile(samples, 0.025),
        "ci95_high": _quantile(samples, 0.975),
        "draws": draws,
        "seed": seed,
    }


def three_replica_descriptive(values: dict[str, float]) -> dict[str, float | int]:
    if set(values) != set(REPLICAS):
        raise ValueError("descriptive summary requires exactly rep01, rep02, and rep03")
    ordered = [float(values[replica]) for replica in REPLICAS]
    if any(not math.isfinite(value) for value in ordered):
        raise ValueError("replica means must be finite")
    return {
        "n_realizations": 3,
        "mean": statistics.fmean(ordered),
        "sample_sd": statistics.stdev(ordered),
    }


def summarize(run_dir: Path, completion_report: Path, output_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    completion = json.loads(completion_report.read_text(encoding="utf-8"))
    if completion.get("status") != "pass":
        raise ValueError("formal run completion report has not passed")
    output_dir.mkdir(parents=True)
    raw: dict[str, dict[str, Any]] = {}
    long_rows = []
    decomp_long_rows = []
    for model in MODELS:
        raw[model] = {}
        expected_terms: set[str] | None = None
        for replica in REPLICAS:
            path = run_dir / model / replica / "FINAL_RESULTS_MMPBSA.csv"
            if not path.is_file():
                raise ValueError(f"missing formal frame-energy CSV: {model}/{replica}")
            rows = parse_delta_energy_csv(path.read_text(encoding="utf-8", errors="replace"))
            decomp_path = run_dir / model / replica / "FINAL_DECOMP_MMPBSA.csv"
            if not decomp_path.is_file():
                raise ValueError(f"missing formal decomposition CSV: {model}/{replica}")
            decomp_rows = parse_decomposition_csv(
                decomp_path.read_text(encoding="utf-8", errors="replace")
            )
            decomp_summary = fixed_decomposition_summary(decomp_rows)
            blocks = fixed_block_means(rows)
            terms = set(blocks[0]["means"])
            if expected_terms is None:
                expected_terms = terms
            elif terms != expected_terms:
                raise ValueError(f"{model}: frame-energy terms drifted between replicas")
            window_means = {term: statistics.fmean(float(row[term]) for row in rows) for term in sorted(terms)}
            raw[model][replica] = {
                "source": str(path),
                "source_sha256": file_sha256(path),
                "frame_count": 300,
                "window_means": window_means,
                "blocks": blocks,
                "fixed_residue_decomposition": {
                    "source": str(decomp_path),
                    "source_sha256": file_sha256(decomp_path),
                    "descriptive_only": True,
                    "data_driven_ranking": False,
                    "residues_in_frozen_order": list(FROZEN_DECOMP_RESIDUES),
                    "residues": decomp_summary,
                },
            }
            for block in blocks:
                for term, value in block["means"].items():
                    long_rows.append({
                        "model": model,
                        "replica": replica,
                        "block_index_zero_based": block["block_index_zero_based"],
                        "start_ns": block["start_ns"],
                        "end_ns": block["end_ns"],
                        "term": term,
                        "mean_kcal_per_mol": value,
                    })
            for residue in FROZEN_DECOMP_RESIDUES:
                for block in decomp_summary[residue]["blocks"]:
                    for term, value in block["means"].items():
                        decomp_long_rows.append({
                            "model": model,
                            "replica": replica,
                            "residue": residue,
                            "block_index_zero_based": block["block_index_zero_based"],
                            "start_ns": block["start_ns"],
                            "end_ns": block["end_ns"],
                            "term": term,
                            "mean_kcal_per_mol": value,
                        })
    inference = {}
    for model in MODELS:
        terms = sorted(raw[model]["rep01"]["window_means"])
        inference[model] = {}
        for term in terms:
            replica_blocks = {
                replica: [block["means"][term] for block in raw[model][replica]["blocks"]]
                for replica in REPLICAS
            }
            inference[model][term] = {
                "replica_window_means": {
                    replica: raw[model][replica]["window_means"][term] for replica in REPLICAS
                },
                "all_three_descriptive": three_replica_descriptive({
                    replica: raw[model][replica]["window_means"][term]
                    for replica in REPLICAS
                }),
                "hierarchical_block_bootstrap": hierarchical_block_bootstrap(
                    replica_blocks, seed=20260818, draws=10000
                ),
            }
    decomposition_descriptive = {}
    for model in MODELS:
        decomposition_descriptive[model] = {}
        for residue in FROZEN_DECOMP_RESIDUES:
            decomposition_descriptive[model][residue] = {}
            terms = list(
                raw[model]["rep01"]["fixed_residue_decomposition"]["residues"]
                [residue]["window_means"]
            )
            for term in terms:
                replica_means = {
                    replica: raw[model][replica]["fixed_residue_decomposition"]
                    ["residues"][residue]["window_means"][term]
                    for replica in REPLICAS
                }
                decomposition_descriptive[model][residue][term] = {
                    "replica_window_means": replica_means,
                    "all_three_descriptive": three_replica_descriptive(replica_means),
                }
    summary = {
        "schema_version": "1.0",
        "report_type": "secondary_exploratory_endpoint_energy_summary",
        "status": "pass",
        "window_ns": [200.0, 500.0],
        "replicas_are_independent_units": True,
        "frames_are_independent_units": False,
        "frame_level_p_values": False,
        "decomposition_is_descriptive_only": True,
        "decomposition_data_driven_ranking": False,
        "models": raw,
        "inference": inference,
        "fixed_residue_decomposition_descriptive": decomposition_descriptive,
        "interpretation_ceiling": (
            "Exploratory model-dependent endpoint-energy summaries only; not absolute binding free energy, "
            "affinity, potency, target occupancy, efficacy, or between-ligand ranking."
        ),
        "formal_completion_sha256": file_sha256(completion_report),
    }
    summary_path = output_dir / "SECONDARY_ENDPOINT_ENERGY_SUMMARY.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    csv_path = output_dir / "SECONDARY_ENDPOINT_ENERGY_FIXED_BLOCKS.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(long_rows[0]))
        writer.writeheader()
        writer.writerows(long_rows)
    decomp_csv_path = output_dir / "SECONDARY_ENDPOINT_ENERGY_FIXED_RESIDUE_DECOMPOSITION_BLOCKS.csv"
    with decomp_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(decomp_long_rows[0]))
        writer.writeheader()
        writer.writerows(decomp_long_rows)
    manifest = {
        "status": "pass",
        "outputs": {
            summary_path.name: {"sha256": file_sha256(summary_path), "bytes": summary_path.stat().st_size},
            csv_path.name: {"sha256": file_sha256(csv_path), "bytes": csv_path.stat().st_size},
            decomp_csv_path.name: {
                "sha256": file_sha256(decomp_csv_path),
                "bytes": decomp_csv_path.stat().st_size,
            },
        },
    }
    manifest_path = output_dir / "SUMMARY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = file_sha256(manifest_path)
    manifest_path.with_suffix(".json.sha256").write_text(f"{digest}  {manifest_path.name}\n", encoding="ascii")
    return {"status": "pass", "manifest": str(manifest_path), "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--completion-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.run_dir, args.completion_report, args.output_dir), sort_keys=True))


if __name__ == "__main__":
    main()
