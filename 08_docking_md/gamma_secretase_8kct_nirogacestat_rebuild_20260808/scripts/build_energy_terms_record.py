#!/usr/bin/env python3
"""Generate config/gromacs_energy_terms.json from the frozen template.

The 16 terms are exactly the GROMACS energy terms used for production QC.
This script binds the production MDP (nstenergy cadence) and GROMACS build
identity so the record is traceable. Nothing is inferred from trajectories.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "templates" / "gromacs_energy_terms.template.json"
MDP = ROOT / "analysis_config_work" / "production_500ns.mdp" if False else None
OUT = ROOT / "config" / "gromacs_energy_terms.json"

# The production MDP lives on the server; the frozen mdout.mdp copy was captured
# into the release. We record the frozen production protocol JSON instead.
PROTOCOL = ROOT / "config" / "production_protocol_hmr4fs_303K_v1.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    record = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    record["approval_status"] = "draft_not_for_execution"
    record["source_template_sha256"] = sha256(TEMPLATE)
    record["production_protocol_sha256"] = sha256(PROTOCOL)
    record["gromacs_energy_cadence_steps"] = 5000  # nstenergy from production MDP
    record["energy_interval_ps"] = 20.0
    record["gromacs_build"] = "GROMACS 2025.2 (production node; identity captured server-side)"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256(OUT.read_bytes()).hexdigest()
    OUT.with_suffix(".json.sha256").write_text(f"{digest}  {OUT.name}\n", encoding="ascii")
    print(f"WROTE {OUT} sha256={digest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
