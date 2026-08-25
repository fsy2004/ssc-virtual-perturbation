#!/usr/bin/env python3
"""Prescreen frozen O6U water orientations by chemical role before visual review."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem


EXPECTED_SMILES = "CCC[C@H](N[C@H]1CCc2cc(F)cc(F)c2C1)C(=O)Nc1cn(C(C)(C)CNCC(C)(C)C)cn1"
EXPECTED_ORIENTATION_COUNT = 70
EXPECTED_RETAINED_COUNT = 20
EXPECTED_EXCLUDED_COUNT = 50

ACCEPTOR_TARGETS = {
    "N1": "secondary tetrahydroisoquinoline amine acceptor",
    "N3": "secondary aliphatic amine acceptor",
    "N5": "pyridine-like imidazole nitrogen acceptor",
    "O": "amide carbonyl oxygen acceptor",
}
DONOR_H_TARGETS = {
    "H39": "secondary tetrahydroisoquinoline N-H donor",
    "H40": "secondary aliphatic N-H donor",
    "H41": "amide N-H donor adjacent to the high-penalty charge region",
}
NONACCEPTOR_N_TARGETS = {
    "N2": "N-substituted pyrrole-like imidazole nitrogen; not an acceptor or donor",
    "N4": "amide nitrogen; donor through H41 but not an acceptor at nitrogen",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orientation-da", required=True, type=Path)
    parser.add_argument("--sdf", required=True, type=Path)
    parser.add_argument("--correspondence-tsv", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    paths = {
        "orientation_da": args.orientation_da.resolve(),
        "sdf": args.sdf.resolve(),
        "correspondence_tsv": args.correspondence_tsv.resolve(),
        "policy": args.policy.resolve(),
    }
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty {label}: {path}")

    policy = json.loads(paths["policy"].read_text(encoding="utf-8"))
    if (
        policy.get("report_type") != "o6u_water_probe_disposition_policy"
        or policy.get("status") != "pass"
        or policy.get("production_approved") is not False
    ):
        raise SystemExit("Disposition policy is invalid")

    supplier = Chem.SDMolSupplier(str(paths["sdf"]), removeHs=False, sanitize=True)
    molecule = supplier[0] if supplier and len(supplier) else None
    if molecule is None:
        raise SystemExit("RDKit could not read the frozen O6U SDF")
    smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
    if smiles != EXPECTED_SMILES or molecule.GetNumAtoms() != 76:
        raise SystemExit(f"Frozen O6U identity differs: atoms={molecule.GetNumAtoms()} smiles={smiles}")

    with paths["correspondence_tsv"].open("r", encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    if len(mapping) != 76:
        raise SystemExit(f"Expected 76 atom-correspondence rows, found {len(mapping)}")
    atom_by_name = {row["cgenff_atom_name"]: row for row in mapping}
    if len(atom_by_name) != 76:
        raise SystemExit("CGenFF atom names are not unique")

    definitions = [
        line.strip()
        for line in paths["orientation_da"].read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(definitions) != EXPECTED_ORIENTATION_COUNT:
        raise SystemExit(f"Expected 70 orientation rows, found {len(definitions)}")

    rows: list[dict[str, object]] = []
    for index, definition in enumerate(definitions, start=1):
        fields = definition.split()
        probe_type, target = fields[0], fields[2]
        if target not in atom_by_name:
            raise SystemExit(f"Unknown target atom on row {index}: {target}")
        if target in ACCEPTOR_TARGETS:
            suggestion = "retain_for_visual_review"
            rationale = ACCEPTOR_TARGETS[target]
            interaction_role = "acceptor"
        elif target in DONOR_H_TARGETS:
            suggestion = "retain_for_visual_review"
            rationale = DONOR_H_TARGETS[target]
            interaction_role = "donor"
        elif target in NONACCEPTOR_N_TARGETS:
            suggestion = "exclude_chemical_role"
            rationale = NONACCEPTOR_N_TARGETS[target]
            interaction_role = "nonacceptor_nitrogen"
        elif target in {"F1", "F2"}:
            suggestion = "exclude_unselected_zero_penalty_site"
            rationale = "organofluorine site outside the prespecified high-penalty charge region"
            interaction_role = "unselected_organofluorine"
        elif target.startswith("H"):
            suggestion = "exclude_nonpolar_site"
            rationale = "carbon-bound hydrogen; not a polar hydrogen-bond donor target"
            interaction_role = "nonpolar_carbon_bound_hydrogen"
        else:
            raise SystemExit(f"No prespecified chemical-role rule for row {index}: {definition}")
        rows.append(
            {
                "orientation_id": f"O6U_WP_{index:03d}",
                "source_line_number": index,
                "source_definition": definition,
                "probe_type": probe_type,
                "target_atom": target,
                "target_element": atom_by_name[target]["element"],
                "target_charge_penalty": float(atom_by_name[target]["charge_penalty"]),
                "interaction_role": interaction_role,
                "prescreen_suggestion": suggestion,
                "prescreen_rationale": rationale,
                "final_disposition": "pending_visual_review" if suggestion == "retain_for_visual_review" else "pending_signed_exclusion",
            }
        )

    counts = dict(sorted(collections.Counter(row["prescreen_suggestion"] for row in rows).items()))
    retained = counts.get("retain_for_visual_review", 0)
    excluded = len(rows) - retained
    if retained != EXPECTED_RETAINED_COUNT or excluded != EXPECTED_EXCLUDED_COUNT:
        raise SystemExit(f"Unexpected prescreen counts: retained={retained} excluded={excluded}")

    report = {
        "schema_version": "1.0",
        "report_type": "o6u_water_probe_chemical_role_prescreen",
        "status": "pass_chemical_role_prescreen_visual_review_required",
        "production_approved": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "component_id": "O6U",
        "frozen_identity": {"atom_count": 76, "isomeric_smiles": smiles},
        "inputs": {label: artifact(path) for label, path in paths.items()},
        "orientation_count": len(rows),
        "prescreen_counts": counts,
        "retained_for_visual_review_count": retained,
        "chemically_excludable_count": excluded,
        "decision_boundary": (
            "This report is a prospective chemical-role prescreen, not a final disposition table. "
            "The 20 retained rows require geometry-specific visual review after the MP2 geometry is "
            "frozen; all exclusions require signed transfer into the final 70-row disposition report."
        ),
        "orientations": rows,
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "output": str(output), "sha256": sha256(output), "counts": counts}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
