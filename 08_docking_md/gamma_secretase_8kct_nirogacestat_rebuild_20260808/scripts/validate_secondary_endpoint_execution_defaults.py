#!/usr/bin/env python3
"""Validate the versioned endpoint-energy execution supplement and sources."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


DEFAULTS_SHA256 = "e8e3d19e133e89679e9c9e44b158d477ef1e6d6d7a893fb4cbe58d5d1db95381"
SUPPLEMENT_SHA256 = "3639bd54694ee11fe684d9940bec19d55a7c3c3e91b8e290bd6b7d50a2de39be"
PLIP_SHA256 = "69c57604a13cb454ee55bd98b4614eefaac0527fc0a2636351393aefb874969d"
PRINT_RES = "B/261,268,272,282,287,380-381,431-432,502"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, field: str) -> None:
    if not condition:
        raise ValueError(f"{field}: frozen execution value mismatch")


def validate(defaults_path: Path, supplement_path: Path, plip_path: Path) -> dict[str, Any]:
    _require(file_sha256(defaults_path) == DEFAULTS_SHA256, "defaults_sha256")
    _require(file_sha256(supplement_path) == SUPPLEMENT_SHA256, "supplement_sha256")
    _require(file_sha256(plip_path) == PLIP_SHA256, "plip_sha256")
    defaults = json.loads(defaults_path.read_text(encoding="utf-8"))
    _require(defaults.get("status") == "frozen_before_endpoint_energy_results", "status")
    software = defaults.get("software", {})
    _require(software == {
        "gmx_mmpbsa_version": "1.6.5",
        "gmx_mmpbsa_git_commit": "64e994c71aaff315f3c82dd0852919aecb1ab62e",
        "python": "3.11.8",
        "ambertools": "23.3",
        "gromacs": "2023.4",
        "pbradii": 7,
        "openmpi": "4.1.6",
        "parmed": "4.3.0",
    }, "software")
    _require(defaults.get("pb_numerical") == {
        "radiopt": 0, "fillratio": 1.25, "inp": 2, "sasopt": 0,
        "solvopt": 2, "ipb": 1, "bcopt": 10, "nfocus": 1,
        "linit": 1000, "eneopt": 1, "cutfd": 7.0, "cutnb": 99.0,
        "maxarcdot": 15000, "npbverb": 1,
    }, "pb_numerical")
    decomp = defaults.get("decomposition", {})
    _require(decomp.get("idecomp") == 2, "decomposition.idecomp")
    _require(decomp.get("dec_verbose") == 0, "decomposition.dec_verbose")
    _require(decomp.get("csv_format") == 1, "decomposition.csv_format")
    _require(decomp.get("print_res") == PRINT_RES, "decomposition.print_res")
    _require(decomp.get("descriptive_only") is True, "decomposition.descriptive_only")
    _require(decomp.get("data_driven_residue_ranking") is False, "decomposition.data_driven_residue_ranking")
    plip = json.loads(plip_path.read_text(encoding="utf-8"))
    observed = sorted({
        (str(row["protein"]["chain"]), str(row["protein"]["residue_name"]), int(row["protein"]["residue_number"]))
        for row in plip.get("interactions", [])
    })
    expected = sorted([
        ("B", "VAL", 261), ("B", "LEU", 268), ("B", "VAL", 272),
        ("B", "LEU", 282), ("B", "ILE", 287), ("B", "LYS", 380),
        ("B", "LEU", 381), ("B", "ALA", 431), ("B", "LEU", 432),
    ])
    _require(observed == expected, "PLIP protein residue set")
    binding = plip.get("binding_site", {})
    _require((binding.get("chain"), binding.get("hetid"), int(binding.get("position", -1))) == ("B", "O6U", 502), "PLIP ligand residue")
    return {
        "schema_version": "1.0",
        "status": "pass",
        "defaults_sha256": DEFAULTS_SHA256,
        "supplement_sha256": SUPPLEMENT_SHA256,
        "plip_sha256": PLIP_SHA256,
        "decomposition_print_res": PRINT_RES,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--defaults", type=Path, required=True)
    parser.add_argument("--supplement", type=Path, required=True)
    parser.add_argument("--plip", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate(args.defaults, args.supplement, args.plip)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report is not None:
        if args.report.exists():
            raise FileExistsError(f"refusing to overwrite {args.report}")
        args.report.write_text(rendered, encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
