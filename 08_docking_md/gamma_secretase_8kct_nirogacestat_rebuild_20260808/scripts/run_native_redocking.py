#!/usr/bin/env python3
"""Run the preregistered 8KCT/O6U self-redocking protocol.

This is a protocol-QA calculation only. It must not be interpreted as an
affinity calculation or as evidence that the experimental pose is stable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SEEDS = (11111, 22222, 33333)
BOX = {
    "center_x": 164.105,
    "center_y": 174.6835,
    "center_z": 142.9475,
    "size_x": 19.594,
    "size_y": 18.381,
    "size_z": 18.711,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pdbqt_charges(path: Path) -> dict[str, float | int]:
    charges: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            charges.append(float(line[70:76]))
    if not charges:
        raise ValueError(f"no PDBQT atoms in {path}")
    return {
        "atom_count": len(charges),
        "rounded_partial_charge_sum_e": round(sum(charges), 6),
        "exact_zero_partial_charge_count": sum(abs(q) < 5e-7 for q in charges),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vina", type=Path, required=True)
    ap.add_argument("--receptor", type=Path, required=True)
    ap.add_argument("--ligand", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()

    for path in (args.vina, args.receptor, args.ligand):
        if not path.is_file():
            ap.error(f"missing required file: {path}")
    if args.cpu < 1:
        ap.error("--cpu must be positive")

    ligand_q = pdbqt_charges(args.ligand)
    if ligand_q["exact_zero_partial_charge_count"]:
        raise SystemExit("NO-GO: ligand PDBQT contains zero partial charges")

    version = subprocess.run(
        [str(args.vina), "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    if version != "AutoDock Vina v1.2.7":
        raise SystemExit(f"NO-GO: expected AutoDock Vina v1.2.7, observed {version!r}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    records = []
    for seed in SEEDS:
        run_dir = args.output_root / f"seed{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)
        output = run_dir / "O6U_redocked.pdbqt"
        stdout_log = run_dir / "vina.stdout.log"
        stderr_log = run_dir / "vina.stderr.log"
        command = [
            str(args.vina),
            "--receptor", str(args.receptor),
            "--ligand", str(args.ligand),
            "--center_x", str(BOX["center_x"]),
            "--center_y", str(BOX["center_y"]),
            "--center_z", str(BOX["center_z"]),
            "--size_x", str(BOX["size_x"]),
            "--size_y", str(BOX["size_y"]),
            "--size_z", str(BOX["size_z"]),
            "--scoring", "vina",
            "--exhaustiveness", "32",
            "--num_modes", "20",
            "--energy_range", "5",
            "--seed", str(seed),
            "--cpu", str(args.cpu),
            "--out", str(output),
        ]
        started = datetime.now(timezone.utc).isoformat()
        proc = subprocess.run(command, capture_output=True, text=True)
        stdout_log.write_text(proc.stdout, encoding="utf-8")
        stderr_log.write_text(proc.stderr, encoding="utf-8")
        if proc.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
            raise SystemExit(f"NO-GO: seed {seed} failed; see {run_dir}")
        records.append(
            {
                "seed": seed,
                "started_utc": started,
                "finished_utc": datetime.now(timezone.utc).isoformat(),
                "command_argv": command,
                "returncode": proc.returncode,
                "output": str(output.resolve()),
                "output_sha256": sha256(output),
                "stdout_sha256": sha256(stdout_log),
                "stderr_sha256": sha256(stderr_log),
            }
        )

    manifest = {
        "schema_version": 1,
        "purpose": "native-pose self-redocking protocol QA only",
        "interpretation_prohibited": [
            "binding affinity",
            "binding stability",
            "potency",
            "mechanism",
        ],
        "software": {"version": version, "executable_sha256": sha256(args.vina)},
        "inputs": {
            "receptor": str(args.receptor.resolve()),
            "receptor_sha256": sha256(args.receptor),
            "receptor_charge_audit": pdbqt_charges(args.receptor),
            "ligand": str(args.ligand.resolve()),
            "ligand_sha256": sha256(args.ligand),
            "ligand_charge_audit": ligand_q,
        },
        "protocol": {
            "box_source": "native O6U heavy-atom envelope plus 5 A padding per Meeko tutorial",
            "box_angstrom": BOX,
            "scoring": "vina",
            "exhaustiveness": 32,
            "num_modes": 20,
            "energy_range_kcal_mol": 5,
            "seeds": list(SEEDS),
            "cpu": args.cpu,
        },
        "runs": records,
    }
    out = args.output_root / "redocking_run_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
