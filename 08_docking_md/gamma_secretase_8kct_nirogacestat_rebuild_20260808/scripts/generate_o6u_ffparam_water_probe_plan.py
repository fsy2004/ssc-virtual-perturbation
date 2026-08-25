#!/usr/bin/env python3
"""Generate and validate the frozen O6U FFParam water-probe orientation plan."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import pathlib
from datetime import datetime, timezone

from ffparam.script_core import qmfileio, toppario
from ffparam.script_core.CoordinateWriter import crdformat
from ffparam.script_core.moleculereader import Molecule
from ffparam.script_core.rtftopsf import rtftopsf


EXPECTED_TYPE_COUNTS = {
    "A2": 2,
    "A31": 16,
    "AP": 2,
    "APL": 6,
    "D": 38,
    "DOP": 6,
}
EXPECTED_TARGETS = {
    "F1",
    "F2",
    "N1",
    "N2",
    "N3",
    "N4",
    "N5",
    "O",
    *(f"H{i}" for i in range(1, 42)),
}


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: pathlib.Path) -> dict:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rtf", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument("--output-da", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rtf = pathlib.Path(args.rtf)
    coordinates = pathlib.Path(args.coordinates)
    output = pathlib.Path(args.output_da)
    report_path = pathlib.Path(args.report)
    for path in (rtf, coordinates):
        if not path.is_file():
            raise FileNotFoundError(path)

    topology = rtftopsf(topparlist=[str(rtf)], resilist=["o6u"], psf=False).topparPsfs[0]
    mapped_indices = toppario.createmap(topology["atomnames"])["mapindices"]
    molecule = Molecule()
    if molecule.readcoor(str(coordinates), center=False) is None:
        raise RuntimeError(f"FFParam could not parse coordinates: {coordinates}")
    formatted = crdformat(molecule.mergepos(), mapped_indices)
    if len(topology["atomnames"]) != 76 or len(formatted) != 76:
        raise RuntimeError("Frozen O6U requires 76 topology and coordinate atoms")
    if len(topology["bondnames"]) != 78:
        raise RuntimeError("Frozen O6U requires 78 bonds")

    output.parent.mkdir(parents=True, exist_ok=True)
    qmfileio.createda(topology, str(output), atomchoice=[], coords=formatted)
    lines = [line.strip() for line in output.read_text(encoding="utf-8").splitlines() if line.strip()]
    parsed = [line.split() for line in lines]
    type_counts = dict(sorted(collections.Counter(fields[0] for fields in parsed).items()))
    targets = {fields[2] for fields in parsed}
    if len(lines) != 70:
        raise RuntimeError(f"Expected 70 frozen O6U probe orientations, found {len(lines)}")
    if type_counts != EXPECTED_TYPE_COUNTS:
        raise RuntimeError(f"Unexpected probe-type counts: {type_counts}")
    if targets != EXPECTED_TARGETS:
        raise RuntimeError(
            f"Probe target mismatch; missing={sorted(EXPECTED_TARGETS-targets)}, "
            f"extra={sorted(targets-EXPECTED_TARGETS)}"
        )
    report = {
        "schema_version": "1.0",
        "status": "pass",
        "production_approved": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "component_id": "O6U",
        "formal_charge_e": 0,
        "topology_atom_count": 76,
        "coordinate_atom_count": 76,
        "bond_count": 78,
        "orientation_count": len(lines),
        "orientation_type_counts": type_counts,
        "target_atoms": sorted(targets),
        "inputs": {"rtf": file_record(rtf), "coordinates": file_record(coordinates)},
        "output_da": file_record(output),
        "interpretation": (
            "This is an exhaustive FFParam orientation-generation artifact, not an accepted "
            "interaction set. Every orientation must receive a source-hashed prospective "
            "disposition. Retained applicable polar interactions require HF/6-31G(d) distance "
            "optimization; prespecified nonpolar, sterically clashing, or competing orientations "
            "may be excluded with retained visual or chemical review evidence."
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(report_path)
    print(sha256_file(report_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
