#!/usr/bin/env python3
"""Audit deposited 8KCT facts before the internal CHARMM-GUI PDB Reader review.

The output is deliberately non-approving: it verifies deposited components,
resolved polymer segments, and covalent-link provenance, then preserves the
remaining primary-plus-independent technical decisions as hard blockers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from collections import Counter
from pathlib import Path

import gemmi


EXPECTED_SEGMENTS = {
    "A": [(34, 700)],
    "B": [(76, 291), (377, 467)],
    "C": [(2, 244)],
    "D": [(6, 101)],
}
EXPECTED_COMPONENTS = {"O6U": 1, "CLR": 3, "PC1": 2, "NAG": 18, "BMA": 3}
EXPECTED_DISULFIDES = {
    ("A", "CYS", "50", "SG", "A", "CYS", "62", "SG"),
    ("A", "CYS", "140", "SG", "A", "CYS", "159", "SG"),
    ("A", "CYS", "230", "SG", "A", "CYS", "248", "SG"),
    ("A", "CYS", "586", "SG", "A", "CYS", "620", "SG"),
}
EXPECTED_N_LINKED_SITES = {45, 55, 187, 264, 387, 435, 464, 506, 530, 562, 573, 580}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contiguous_segments(numbers: list[int]) -> list[tuple[int, int]]:
    require(numbers, "Cannot derive segments from an empty residue list")
    require(numbers == sorted(set(numbers)), "Polymer residue numbers are duplicated or unordered")
    segments: list[tuple[int, int]] = []
    start = previous = numbers[0]
    for number in numbers[1:]:
        if number != previous + 1:
            segments.append((start, previous))
            start = number
        previous = number
    segments.append((start, previous))
    return segments


def struct_conn_records(block: gemmi.cif.Block) -> list[dict[str, str]]:
    tags = [
        "_struct_conn.id",
        "_struct_conn.conn_type_id",
        "_struct_conn.ptnr1_label_asym_id",
        "_struct_conn.ptnr1_label_comp_id",
        "_struct_conn.ptnr1_label_seq_id",
        "_struct_conn.ptnr1_label_atom_id",
        "_struct_conn.ptnr2_label_asym_id",
        "_struct_conn.ptnr2_label_comp_id",
        "_struct_conn.ptnr2_label_seq_id",
        "_struct_conn.ptnr2_label_atom_id",
    ]
    keys = [tag.split(".")[-1] for tag in tags]
    return [{key: str(value) for key, value in zip(keys, row)} for row in block.find(tags)]


def audit_structure(structure_path: Path, validation_report: Path) -> dict[str, object]:
    document = gemmi.cif.read_file(str(structure_path))
    block = document.sole_block()
    structure = gemmi.make_structure_from_block(block)
    require(len(structure) == 1, f"Expected one 8KCT model, observed {len(structure)}")
    model = structure[0]

    polymer_segments: dict[str, list[tuple[int, int]]] = {}
    for chain in model:
        numbers = [residue.seqid.num for residue in chain if residue.het_flag == "A"]
        if numbers:
            polymer_segments[chain.name] = contiguous_segments(numbers)
    require(polymer_segments == EXPECTED_SEGMENTS, f"Resolved polymer segments differ: {polymer_segments}")

    components = Counter(
        residue.name
        for chain in model
        for residue in chain
        if residue.name in EXPECTED_COMPONENTS
    )
    require(dict(components) == EXPECTED_COMPONENTS, f"Native component counts differ: {dict(components)}")
    waters = [residue for chain in model for residue in chain if residue.name in {"HOH", "WAT", "TIP3"}]
    require(not waters, f"Unexpected resolved waters: {len(waters)}")

    all_atoms = [atom for chain in model for residue in chain for atom in residue]
    require(all(float(atom.occ) >= 0.99 for atom in all_atoms), "At least one deposited atom has occupancy below 0.99")
    alternate_atoms = [
        {"chain": chain.name, "residue": residue.seqid.num, "resname": residue.name, "atom": atom.name, "altloc": atom.altloc}
        for chain in model
        for residue in chain
        for atom in residue
        if atom.altloc not in {"\x00", " "}
    ]
    require(not alternate_atoms, f"Deposited model contains alternate locations: {alternate_atoms[:5]}")

    connections = struct_conn_records(block)
    disulfides = [record for record in connections if record["conn_type_id"] == "disulf"]
    observed_disulfides = {
        (
            record["ptnr1_label_asym_id"], record["ptnr1_label_comp_id"], record["ptnr1_label_seq_id"], record["ptnr1_label_atom_id"],
            record["ptnr2_label_asym_id"], record["ptnr2_label_comp_id"], record["ptnr2_label_seq_id"], record["ptnr2_label_atom_id"],
        )
        for record in disulfides
    }
    require(observed_disulfides == EXPECTED_DISULFIDES, f"Disulfide records differ: {observed_disulfides}")

    covalent = [record for record in connections if record["conn_type_id"] == "covale"]
    n_linked = [
        record for record in covalent
        if record["ptnr1_label_comp_id"] == "ASN"
        and record["ptnr1_label_atom_id"] == "ND2"
        and record["ptnr2_label_comp_id"] == "NAG"
        and record["ptnr2_label_atom_id"] == "C1"
    ]
    n_linked_sites = {int(record["ptnr1_label_seq_id"]) for record in n_linked}
    require(n_linked_sites == EXPECTED_N_LINKED_SITES, f"N-linked glycosylation sites differ: {sorted(n_linked_sites)}")
    glycan_internal = [record for record in covalent if record not in n_linked]
    require(len(covalent) == 21, f"Expected 21 glycan covalent links, observed {len(covalent)}")
    require(len(n_linked) == 12 and len(glycan_internal) == 9, "Expected 12 protein-NAG plus 9 glycan-internal links")

    o6u = [(chain.name, residue) for chain in model for residue in chain if residue.name == "O6U"]
    require(len(o6u) == 1 and o6u[0][0] == "B" and o6u[0][1].seqid.num == 502, "Native O6U site differs from B:502")
    require(len(o6u[0][1]) == 35, "Native O6U does not have 35 deposited heavy atoms")

    return {
        "schema_version": "1.0",
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "overall_status": "deposited_structure_pass_manual_pdb_reader_blocked",
        "pdb_reader_approved": False,
        "quick_bilayer_submission_allowed": False,
        "source": {
            "pdb_id": "8KCT",
            "path": str(structure_path.resolve()),
            "bytes": structure_path.stat().st_size,
            "sha256": sha256(structure_path),
            "wwpdb_validation_report": {
                "path": str(validation_report.resolve()),
                "bytes": validation_report.stat().st_size,
                "sha256": sha256(validation_report),
            },
        },
        "deposited_model": {
            "model_count": len(structure),
            "atom_count": len(all_atoms),
            "minimum_occupancy": min(float(atom.occ) for atom in all_atoms),
            "alternate_location_atom_count": len(alternate_atoms),
            "resolved_polymer_segments": {chain: [[start, end] for start, end in segments] for chain, segments in polymer_segments.items()},
            "unresolved_psen1_interval": [292, 376],
            "native_component_counts": dict(components),
            "resolved_water_count": len(waters),
            "native_o6u": {"chain": "B", "residue_number": 502, "heavy_atom_count": 35},
        },
        "deposited_covalent_connections": {
            "total_struct_conn_records": len(connections),
            "disulfide_count": len(disulfides),
            "disulfides": disulfides,
            "glycan_covalent_count": len(covalent),
            "protein_to_nag_count": len(n_linked),
            "protein_to_nag_sites": sorted(n_linked_sites),
            "glycan_internal_count": len(glycan_internal),
            "protein_to_nag_records": n_linked,
            "glycan_internal_records": glycan_internal,
        },
        "manual_blockers": [
            "PDB Reader segment and terminal patches have not been reviewed and signed",
            "PSEN1 NTF and CTF must remain separate topology segments across the unresolved 292-376 interval",
            "all four deposited disulfides and all 21 glycan covalent links must be reproduced in the generated topology",
            "18 NAG, 3 BMA, 3 CLR, 2 PC1/DSPC, and native O6U dispositions must be checked in coordinates and topology",
            "Asp257 deprotonated and Asp385 protonated atom/patch mapping must be recorded",
            "neutral O6U parameters must pass the independent ligand parameter gate",
            "orthogonal visual inspection and an independent reviewer signature are absent",
            "no PDB Reader job ID, sanitized evidence bundle, or curated-coordinate hash exists",
        ],
        "prohibited_actions": [
            "Do not submit Quick Bilayer from the deposited mmCIF alone.",
            "Do not silently strip native glycans, CLR, PC1/DSPC, or O6U.",
            "Do not bridge or model PSEN1 residues 292-376 for this narrow native-pose study.",
            "Do not approve the structure record from this automated audit alone.",
        ],
    }


def checklist_text(audit_path: Path, payload: dict[str, object]) -> str:
    model = payload["deposited_model"]
    links = payload["deposited_covalent_connections"]
    return f"""# PDB Reader manual review checklist\n\nThis checklist is a mandatory human checkpoint. Passing the deposited-structure audit does not approve a CHARMM-GUI model.\n\n## Immutable deposited facts\n\n- Source: 8KCT biological assembly 1, one cryo-EM model.\n- Polymer segments: NCSTN A 34-700; PSEN1 B 76-291 and 377-467; APH1A C 2-244; PEN2 D 6-101.\n- The PSEN1 292-376 interval is unresolved and must not be bridged or newly modelled.\n- Retain exactly: 1 O6U, 3 CLR, 2 PC1/DSPC, 18 NAG, and 3 BMA.\n- Reproduce exactly: {links['disulfide_count']} disulfides, {links['protein_to_nag_count']} protein-to-NAG links, and {links['glycan_internal_count']} glycan-internal links.\n- Native O6U is chain B residue 502 with {model['native_o6u']['heavy_atom_count']} deposited heavy atoms.\n\n## Required PDB Reader review\n\n- [ ] Record the PDB Reader job ID, UTC time, official host, sanitized request/response evidence, and all downloaded file hashes.\n- [ ] Confirm chain-to-segment mapping and keep PSEN1 NTF and CTF as separate topology segments.\n- [ ] Review every N- and C-terminal patch; do not accept automatic terminal choices without recording them.\n- [ ] Confirm all four NCSTN disulfides: 50-62, 140-159, 230-248, and 586-620.\n- [ ] Confirm all 12 protein-to-NAG sites: 45, 55, 187, 264, 387, 435, 464, 506, 530, 562, 573, and 580.\n- [ ] Confirm the nine deposited glycan-internal covalent links; retain 18 NAG and 3 BMA in the final coordinates and topology.\n- [ ] Retain the three resolved CLR and two resolved PC1/DSPC molecules exactly once; record any CHARMM residue-name conversion.\n- [ ] Retain native neutral O6U exactly once and bind its 76-row atom correspondence plus accepted parameter-record hash.\n- [ ] Apply the predeclared dyad state: Asp257 deprotonated and Asp385 protonated; record the exact CHARMM residue/patch and proton atom.\n- [ ] Inspect missing heavy atoms, clashes, peptide geometry, chirality, and generated hydrogens; record every repair or confirm none.\n- [ ] Save orthogonal whole-complex, membrane-domain, glycan, structural-lipid, catalytic-pocket, and O6U close-up images.\n- [ ] Independently compare final coordinate and topology counts against the immutable deposited facts above.\n- [ ] Curator and independent reviewer sign the structure record before Quick Bilayer submission.\n\n## Release rule\n\nQuick Bilayer remains blocked until this checklist, the approved structure record, and the independently approved O6U parameter record are all complete. The source audit is `{audit_path.name}`; bind its SHA-256 in the signed review record.\n"""


def write_outputs(structure_path: Path, validation_report: Path, output_dir: Path) -> tuple[Path, Path]:
    payload = audit_structure(structure_path, validation_report)
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = output_dir / "STRUCTURE_PRECURATION_AUDIT.json"
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    checklist_path = output_dir / "PDB_READER_REVIEW_CHECKLIST.md"
    checklist = checklist_text(audit_path, payload)
    checklist = checklist.replace(
        "# PDB Reader manual review checklist\n\nThis checklist is a mandatory human checkpoint.",
        "# PDB Reader internal technical review checklist\n\n"
        "This checklist is an internal model-release safeguard, not manuscript evidence.",
    ).replace(
        "Curator and independent reviewer sign the structure record",
        "The primary technical curator and an independently tasked reviewer sign the structure record",
    ).replace(
        "Polymer segments: NCSTN A 34-700; PSEN1 B 76-291 and 377-467; APH1A C 2-244; PEN2 D 6-101.",
        "Polymer segments by author-chain ID: NCSTN A 34-700; PSEN1 B 76-291 and 377-467; "
        "APH1A C 2-244; PEN2 D 6-101. Preserve the explicit author-chain-to-label-asym mapping.",
    ).replace(
        "Native O6U is chain B residue 502",
        "Native O6U is author chain B residue 502 (label asym V)",
    ).replace(
        "- [ ] Confirm chain-to-segment mapping and keep PSEN1 NTF and CTF as separate topology segments.",
        "- [ ] Record and verify the complete author-chain-to-label-asym-to-CHARMM-segment mapping; keep PSEN1 NTF and CTF as separate topology segments.\n"
        "- [ ] Create five protein topology segments in total. Explicitly prove that PSEN1 author-chain B residues 291 and 377 are not peptide-bonded; the deposited PDB has no `TER` at this unresolved 292-376 interval.",
    ).replace(
        "- [ ] Inspect missing heavy atoms, clashes, peptide geometry, chirality, and generated hydrogens; record every repair or confirm none.",
        "- [ ] Restore and label as generated the missing PSEN1 Lys76 side-chain atoms CG/CD/CE/NZ and the missing tails of both PC1/DSPC molecules; never describe the deposited model as having no missing heavy atoms.\n"
        "- [ ] Inspect clashes, peptide geometry, chirality, generated hydrogens, and all other missing atoms. Retain before/after reports and quantify deposited-heavy-atom displacement for every repair.\n"
        "- [ ] Review the deposited O6U geometry outliers and the PSEN1 Leu432 O--O6U N04 short contact without moving or re-fitting the experimental pose silently.\n"
        "- [ ] Confirm that PDB Reader introduces no unrecorded waters or ions; later bulk solvent and ions must remain distinguishable as generated components.",
    )
    checklist_path.write_text(checklist, encoding="utf-8", newline="\n")
    return audit_path, checklist_path


def self_test(structure_path: Path, validation_report: Path) -> None:
    payload = audit_structure(structure_path, validation_report)
    require(payload["overall_status"] == "deposited_structure_pass_manual_pdb_reader_blocked", "Valid source did not retain manual block")
    require(payload["quick_bilayer_submission_allowed"] is False, "Automated source audit incorrectly released Quick Bilayer")
    rejected = False
    try:
        audit_structure(validation_report, validation_report)
    except (ValueError, RuntimeError):
        rejected = True
    require(rejected, "A non-mmCIF structure input was not rejected")
    with tempfile.TemporaryDirectory(prefix="8kct_precuration_") as tmp:
        audit_path, checklist_path = write_outputs(structure_path, validation_report, Path(tmp))
        require(audit_path.is_file() and checklist_path.is_file(), "Self-test outputs are missing")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--validation-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    require(args.structure.is_file(), f"Missing structure: {args.structure}")
    require(args.validation_report.is_file(), f"Missing validation report: {args.validation_report}")
    if args.self_test:
        self_test(args.structure, args.validation_report)
        print(json.dumps({"status": "pass", "test": "deposited facts plus invalid-input rejection"}))
        return 0
    require(args.output_dir is not None, "--output-dir is required unless --self-test is used")
    audit_path, checklist_path = write_outputs(args.structure, args.validation_report, args.output_dir)
    print(json.dumps({"status": "pass", "audit": str(audit_path), "checklist": str(checklist_path), "audit_sha256": sha256(audit_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
