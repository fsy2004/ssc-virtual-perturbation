#!/usr/bin/env python3
"""Prepare compact time-series + block-summary data for the MD main figure (Figure 7).
Reads the preserved raw_unsmoothed CSVs per realization and writes a JSON with
downsampled per-frame series and fixed-window block means for plotting.
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path

RUN = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/analysis/primary_postprocessing/20260822T170054Z_primary_qc_rep01_completion"
REPS = ["rep01", "rep02", "rep03"]
OUT = "/root/autodl-tmp/o6u_md_release_3x500ns_v4/analysis/primary_postprocessing/20260822T170054Z_primary_qc_rep01_completion/figure7_data.json"

STRUCT_METRICS = [
    "pocket_aligned_o6u_heavy_rmsd_nm",
    "pocket_aligned_o6u_com_displacement_nm",
    "native_contact_fraction",
    "tm_core_ca_rmsd_nm",
    "protein_ca_rmsd_nm",
]
MEMB_METRICS = ["phosphate_peak_thickness_nm", "protein_tilt_deg", "outside_protein_hydrophobic_core_water_count",
                "cell_lateral_area_nm2_not_apl", "box_z_vector_length_nm", "cell_volume_nm3"]
ENERGY_METRICS = ["temperature_k", "density_kg_m3", "volume_nm3", "box_x_nm", "box_y_nm", "box_z_nm", "pressure_bar"]


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def extract(path: Path, metrics: list[str], step: int = 25) -> dict:
    rows = read_rows(path)
    t = []
    for r in rows:
        key = {}
        for m in metrics:
            try:
                key[m] = float(r[m])
            except (KeyError, ValueError):
                key[m] = None
        key["t"] = float(r["time_ns"])
        key["rep"] = r["realization_id"]
        t.append(key)
    ns = len(t)
    # downsample for the raw curve
    down = t[::step] if step else t
    # block means over 0-500 per 50 ns
    blocks = []
    for lo in range(0, 500, 50):
        sub = [r for r in t if lo <= r["t"] < lo + 50]
        if sub:
            b = {"lo_ns": lo, "hi_ns": lo + 50}
            for m in metrics:
                vals = [r[m] for r in sub if r[m] is not None]
                b[m] = sum(vals) / len(vals) if vals else None
            blocks.append(b)
    return {"n_frames": ns, "downsampled": down, "blocks": blocks}


def main() -> int:
    data = {"run": RUN, "reps": REPS}
    for rep in REPS:
        d = {}
        sp = Path(RUN) / "structural_analysis" / rep / "structural_raw_unsmoothed.csv"
        mp = Path(RUN) / "membrane_qc" / rep / "membrane_raw_unsmoothed.csv"
        ep = Path(RUN) / "energy_qc" / rep / "energy_raw_unsmoothed.csv"
        if sp.exists():
            d["structural"] = extract(sp, STRUCT_METRICS)
        if mp.exists():
            d["membrane"] = extract(mp, MEMB_METRICS)
        if ep.exists():
            d["energy"] = extract(ep, ENERGY_METRICS)
        data[rep] = d
    Path(OUT).write_text(json.dumps(data, indent=1, default=str), encoding="utf-8")
    print(f"wrote {OUT} size={Path(OUT).stat().st_size}")
    print("keys per rep:", list(data["rep01"].keys()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
