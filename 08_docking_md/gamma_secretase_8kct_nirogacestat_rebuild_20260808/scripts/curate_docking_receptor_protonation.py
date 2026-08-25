#!/usr/bin/env python3
"""Reproduce the retired double-deprotonated docking-QA receptor only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_identity(line: str) -> tuple[str, str, str, int]:
    return line[12:16].strip(), line[17:21].strip(), line[21:22], int(line[22:26])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    source_lines = args.input.read_text(encoding="ascii", errors="strict").splitlines()
    output_lines: list[str] = []
    removed: list[str] = []
    normalized = 0
    for line in source_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            output_lines.append(line)
            continue
        atom, resname, chain, resid = atom_identity(line)
        four_character_resname = line[16:20].strip()
        if chain == "B" and resid == 385 and atom == "HD2":
            removed.append(line)
            continue
        if chain == "B" and resid == 385 and four_character_resname == "ASPP":
            line = line[:16] + " ASP" + line[20:]
            normalized += 1
        output_lines.append(line)

    if len(removed) != 1:
        raise SystemExit(f"Expected exactly one ASPP B385 side-chain proton, removed {len(removed)}")
    if normalized and normalized < 4:
        raise SystemExit(f"Too few ASPP B385 records were normalized: {normalized}")
    dyad: dict[str, list[str]] = {"B:257": [], "B:385": []}
    for line in output_lines:
        if not line.startswith("ATOM  "):
            continue
        atom, resname, chain, resid = atom_identity(line)
        key = f"{chain}:{resid}"
        if key in dyad:
            if resname != "ASP":
                raise SystemExit(f"Unexpected dyad residue name after curation: {key} {resname}")
            dyad[key].append(atom)
    for key, atoms in dyad.items():
        if not {"CG", "OD1", "OD2"}.issubset(atoms):
            raise SystemExit(f"Incomplete dyad side chain: {key} {atoms}")
        if any(atom in {"HD1", "HD2"} for atom in atoms):
            raise SystemExit(f"Dyad side-chain proton remains: {key} {atoms}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(output_lines) + "\n", encoding="ascii", newline="\n")
    report = {
        "schema_version": "1.0",
        "purpose": "docking receptor preparation only",
        "model": "PSEN1 Asp257(-1) / Asp385(-1), retired docking-QA-only dyad",
        "input": {"path": str(args.input), "sha256": sha256(args.input)},
        "output": {"path": str(args.output), "sha256": sha256(args.output)},
        "removed_sidechain_protons": [{"chain": "B", "resid": 385, "atom": "HD2"}],
        "normalized_aspp_records": normalized,
        "dyad_atom_names": dyad,
        "warning": "PROPKA predicted Asp385 pKa 7.77. This exact override is retained only to reproduce the failed docking QA and must not enter the production MD build.",
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "pass", "report": str(args.report), "report_sha256": sha256(args.report)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
