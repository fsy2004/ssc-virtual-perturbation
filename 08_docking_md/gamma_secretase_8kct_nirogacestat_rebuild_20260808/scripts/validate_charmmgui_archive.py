#!/usr/bin/env python3
"""Safely audit one downloaded CHARMM-GUI archive without extracting it."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import sys
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

from md_contract import build_contract_sha256, load_json, report_payload_sha256


TEXT_LIMIT = 64 * 1024 * 1024
GRO_RE = re.compile(r"(?:^|/)gromacs/(?:step6\.6_equilibration|step5_input|step4_lipid|system)\.gro$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_members(archive: tarfile.TarFile) -> list[tarfile.TarInfo]:
    members = archive.getmembers()
    if not members:
        raise ValueError("archive is empty")
    for member in members:
        path = PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"unsafe member path: {member.name}")
        if member.issym() or member.islnk():
            raise ValueError(f"archive link is not accepted: {member.name}")
        if member.size < 0:
            raise ValueError(f"invalid member size: {member.name}")
    return members


def read_member_text(archive: tarfile.TarFile, member: tarfile.TarInfo) -> str:
    if member.size > TEXT_LIMIT:
        raise ValueError(f"refusing to read unexpectedly large text member: {member.name}")
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"cannot read archive member: {member.name}")
    return handle.read().decode("utf-8", errors="replace")


def sha256_member(archive: tarfile.TarFile, member: tarfile.TarInfo) -> str:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"cannot read archive member: {member.name}")
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()


def select_gro(members: list[tarfile.TarInfo]) -> tarfile.TarInfo:
    candidates = [member for member in members if member.isfile() and GRO_RE.search(member.name)]
    if not candidates:
        candidates = [
            member
            for member in members
            if member.isfile() and "/gromacs/" in f"/{member.name}" and member.name.endswith(".gro")
        ]
    if not candidates:
        raise ValueError("no GROMACS coordinate file found")
    priority = ("step6.6_equilibration.gro", "step5_input.gro", "system.gro")
    for suffix in priority:
        matches = [member for member in candidates if member.name.endswith(suffix)]
        if matches:
            return sorted(matches, key=lambda item: item.name)[0]
    return sorted(candidates, key=lambda item: item.name)[-1]


def parse_gro(text: str) -> dict[str, Any]:
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError("GRO file is too short")
    try:
        expected_atoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError("GRO atom-count line is invalid") from exc
    if expected_atoms <= 0 or len(lines) < expected_atoms + 3:
        raise ValueError("GRO atom count does not match file length")
    residues: Counter[str] = Counter()
    residue_molecules: Counter[str] = Counter()
    ligand_atom_names: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    previous_residue: tuple[int, str] | None = None
    for index, line in enumerate(lines[2 : 2 + expected_atoms], start=1):
        if len(line) < 44:
            raise ValueError(f"GRO atom line {index} is truncated")
        try:
            resid = int(line[0:5])
        except ValueError as exc:
            raise ValueError(f"GRO residue number on atom line {index} is invalid") from exc
        resname = line[5:10].strip().upper()
        atomname = line[10:15].strip()
        residues[resname] += 1
        residue_key = (resid, resname)
        if residue_key != previous_residue:
            residue_molecules[resname] += 1
            previous_residue = residue_key
        if resname == "O6U":
            ligand_atom_names.append(atomname)
        try:
            coords = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
        except ValueError as exc:
            raise ValueError(f"GRO coordinate line {index} is invalid") from exc
        if not all(math.isfinite(value) for value in coords):
            raise ValueError(f"GRO contains non-finite coordinates at atom line {index}")
        coordinates.append(coords)
    try:
        box = [float(value) for value in lines[2 + expected_atoms].split()]
    except ValueError as exc:
        raise ValueError("GRO box line is invalid") from exc
    if len(box) not in {3, 9} or not all(math.isfinite(value) and value > 0 for value in box[:3]):
        raise ValueError("GRO box vectors are non-positive or malformed")
    return {
        "atom_count": expected_atoms,
        "residue_atom_counts": dict(sorted(residues.items())),
        "residue_molecule_counts": dict(sorted(residue_molecules.items())),
        "o6u_atom_count": residues.get("O6U", 0),
        "o6u_atom_names": ligand_atom_names,
        "box_nm": box,
    }


def parse_itp_charges(text: str, residue: str) -> list[float]:
    section = ""
    charges: list[float] = []
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[] ").lower()
            continue
        if section != "atoms" or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 8:
            continue
        if fields[3].upper() != residue.upper():
            continue
        try:
            value = float(fields[6])
        except ValueError:
            continue
        if not math.isfinite(value):
            raise ValueError("ligand topology contains a non-finite partial charge")
        charges.append(value)
    return charges


def parse_hydrogen_masses(text: str) -> list[float]:
    """Read explicit [ atoms ] masses whose atom name begins with H."""
    section = ""
    masses: list[float] = []
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[] ").lower()
            continue
        if section != "atoms" or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 8 or not fields[4].upper().startswith("H"):
            continue
        try:
            mass = float(fields[7])
        except ValueError:
            continue
        if math.isfinite(mass):
            masses.append(mass)
    return masses


def parse_mdp_text(text: str, label: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        normalized = key.strip().lower().replace("-", "_")
        if normalized in values:
            raise ValueError(f"{label}: duplicate active MDP key {key.strip()!r}")
        values[normalized] = value.strip()
    return values


def mdp_bool(value: str) -> bool:
    return re.sub(r"[-_\s]", "", value.lower()) in {"yes", "true"}


def audit_charmm_gui_stage_mdp(name: str, text: str, temperature_k: float) -> dict[str, Any]:
    mdp = parse_mdp_text(text, name)
    required = {"integrator", "nsteps", "gen_vel"}
    missing = sorted(required.difference(mdp))
    if missing:
        raise ValueError(f"{name}: missing staged MDP fields {missing}")
    integrator = mdp["integrator"].strip().lower()
    nsteps = int(mdp["nsteps"])
    if nsteps <= 0:
        raise ValueError(f"{name}: nsteps must be positive")
    dt = float(mdp.get("dt", "0")) if integrator in {"md", "md-vv", "md_vv"} else None
    if dt is not None and (not math.isfinite(dt) or dt <= 0.0 or dt > 0.002):
        raise ValueError(f"{name}: dynamic stage dt must be in (0, 0.002] ps")
    if name == "step6.0_minimization.mdp":
        if integrator not in {"steep", "cg", "l-bfgs"} or mdp_bool(mdp["gen_vel"]):
            raise ValueError("step6.0 minimization integrator/velocity contract failed")
    else:
        match = re.fullmatch(r"step6\.([1-6])_equilibration\.mdp", name)
        if match is None or integrator not in {"md", "md-vv", "md_vv"}:
            raise ValueError(f"{name}: unexpected CHARMM-GUI dynamic stage")
        index = int(match.group(1))
        if index == 1:
            if not mdp_bool(mdp["gen_vel"]) or mdp_bool(mdp.get("continuation", "no")):
                raise ValueError("step6.1 must generate velocities without continuation")
            if re.sub(r"[-_\s]", "", mdp.get("pcoupl", "no").lower()) != "no":
                raise ValueError("step6.1 must be NVT before pressure coupling")
        else:
            if mdp_bool(mdp["gen_vel"]) or not mdp_bool(mdp.get("continuation", "no")):
                raise ValueError(f"{name}: later equilibration must continue without regenerating velocities")
        ref_t_values = [float(value) for value in mdp.get("ref_t", "").split()]
        if not ref_t_values or any(abs(value - temperature_k) > 1e-9 for value in ref_t_values):
            raise ValueError(f"{name}: ref-t differs from frozen {temperature_k} K")
        if re.sub(r"[-_\s]", "", mdp.get("tcoupl", "no").lower()) == "no":
            raise ValueError(f"{name}: temperature coupling is disabled")
    forbidden = [
        key for key in mdp
        if key == "mass_repartition_factor"
        or key in {"pull", "awh", "free_energy", "simulated_tempering", "annealing", "mts", "deform", "qmmm"}
        or key.startswith(("pull_", "awh_", "fep_", "sim_temp", "electric_field_"))
    ]
    if forbidden:
        raise ValueError(f"{name}: prohibited staged controls present: {sorted(forbidden)}")
    return {
        "name": name,
        "integrator": integrator,
        "dt_ps": dt,
        "nsteps": nsteps,
        "duration_ps": dt * nsteps if dt is not None else None,
        "continuation": mdp_bool(mdp.get("continuation", "no")),
        "gen_vel": mdp_bool(mdp["gen_vel"]),
        "gen_seed": mdp.get("gen_seed"),
        "tcoupl": mdp.get("tcoupl"),
        "pcoupl": mdp.get("pcoupl"),
        "define": mdp.get("define"),
        "constraints": mdp.get("constraints"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Frozen study manifest that owns this build")
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--system-id", default="8kct_nirogacestat_native")
    parser.add_argument("--ligand-resname", default="O6U")
    parser.add_argument("--expected-ligand-charge", type=float, default=0.0)
    parser.add_argument("--expected-ligand-atoms", type=int, default=76)
    parser.add_argument(
        "--expected-artifact",
        action="append",
        default=[],
        metavar="ARCHIVE_MEMBER=SHA256",
        help="require an exact approved ligand/topology artifact and hash; repeat as needed",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--strict", action="store_true", help="treat all warnings as failures")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    archive_path = args.archive.resolve()
    manifest_path = args.manifest.resolve()
    package_root = manifest_path.parent.parent.resolve()
    manifest = load_json(manifest_path)
    systems = manifest.get("systems", [])
    if not isinstance(systems, list) or len(systems) != 1:
        raise SystemExit("manifest must contain exactly one system")
    system = systems[0]
    construction = system.get("construction", {})
    report: dict[str, Any] = {
        "schema_version": "2.0",
        "report_type": "charmm_gui_build_validation",
        "strict": args.strict,
        "study_id": manifest.get("study_id"),
        "system_id": args.system_id,
        "construction_id": construction.get("id"),
        "pdb_reader_jobid": construction.get("pdb_reader_jobid"),
        "quick_bilayer_jobid": construction.get("quick_bilayer_jobid"),
        "build_contract_sha256": build_contract_sha256(manifest),
        "manifest_path": str(manifest_path),
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "archive": str(archive_path),
        "expect_ligand": True,
        "ligand_resname": args.ligand_resname.upper(),
    }
    try:
        if args.system_id != "8kct_nirogacestat_native" or system.get("id") != args.system_id:
            raise ValueError("this frozen validator accepts only 8kct_nirogacestat_native")
        if construction.get("id") != "build01":
            raise ValueError("manifest construction must be build01")
        expected_archive = construction.get("charmm_gui_archive", {})
        expected_archive_path = (package_root / str(expected_archive.get("path", ""))).resolve()
        if archive_path != expected_archive_path:
            raise ValueError("--archive differs from manifest build01 archive path")
        if args.expected_ligand_charge != 0.0 or args.expected_ligand_atoms != 76:
            raise ValueError("the frozen O6U model requires charge 0 and 76 atoms including hydrogens")
        if not archive_path.is_file() or archive_path.stat().st_size == 0:
            raise ValueError(f"archive is missing or empty: {archive_path}")
        report["archive_bytes"] = archive_path.stat().st_size
        report["archive_sha256"] = sha256(archive_path)
        if report["archive_sha256"] != expected_archive.get("sha256"):
            raise ValueError("archive SHA-256 differs from manifest")
        report["membrane_orientation_record_sha256"] = construction.get("membrane_orientation_record", {}).get("sha256")
        report["gromacs_input_tree_manifest_sha256"] = construction.get("gromacs_input_tree_manifest", {}).get("sha256")
        report["starting_coordinates_sha256"] = construction.get("starting_coordinates", {}).get("sha256")
        report["topology_sha256"] = construction.get("topology", {}).get("sha256")
        report["index_sha256"] = construction.get("index", {}).get("sha256")
        report["analysis_index_sha256"] = construction.get("analysis_index", {}).get("sha256")
        report["production_mdp_sha256"] = construction.get("production_mdp", {}).get("sha256")
        report["minimization_mdp_sha256"] = construction.get("minimization_mdp", {}).get("sha256")
        report["equilibration_mdp_sha256"] = construction.get("equilibration_mdp_sha256")
        with tarfile.open(archive_path, mode="r:*") as archive:
            members = safe_members(archive)
            names = [member.name for member in members]
            report["member_count"] = len(members)
            report["safe_member_paths"] = True

            topols = [member for member in members if member.isfile() and member.name.endswith("/gromacs/topol.top")]
            if not topols:
                topols = [member for member in members if member.isfile() and member.name == "gromacs/topol.top"]
            if len(topols) != 1:
                errors.append(f"expected one gromacs/topol.top, found {len(topols)}")
            else:
                report["topology_member"] = topols[0].name

            if len(topols) == 1:
                gromacs_prefix = PurePosixPath(topols[0].name).parent
                extracted_root = (package_root / str(construction.get("gromacs_input_dir", ""))).resolve()
                if extracted_root != package_root and package_root not in extracted_root.parents:
                    errors.append("gromacs_input_dir escapes the package root")
                elif not extracted_root.is_dir():
                    errors.append(f"extracted gromacs_input_dir is missing: {extracted_root}")
                else:
                    tree_binding: list[dict[str, Any]] = []
                    for member in sorted(
                        (value for value in members if value.isfile() and PurePosixPath(value.name).is_relative_to(gromacs_prefix)),
                        key=lambda value: value.name,
                    ):
                        relative = PurePosixPath(member.name).relative_to(gromacs_prefix)
                        if str(relative) == ".":
                            continue
                        extracted = extracted_root.joinpath(*relative.parts).resolve()
                        if extracted_root not in extracted.parents:
                            errors.append(f"archive-to-extracted path escapes gromacs root: {member.name}")
                            continue
                        member_hash = sha256_member(archive, member)
                        extracted_hash = sha256(extracted) if extracted.is_file() else None
                        tree_binding.append({
                            "archive_member": member.name,
                            "relative_path": relative.as_posix(),
                            "archive_bytes": member.size,
                            "archive_sha256": member_hash,
                            "extracted_path": str(extracted),
                            "extracted_bytes": extracted.stat().st_size if extracted.is_file() else None,
                            "extracted_sha256": extracted_hash,
                            "match": extracted.is_file() and extracted.stat().st_size == member.size and extracted_hash == member_hash,
                        })
                    report["archive_extracted_gromacs_tree_binding"] = tree_binding
                    if not tree_binding or not all(item["match"] for item in tree_binding):
                        errors.append("the extracted GROMACS tree is incomplete or differs byte-for-byte from the archive")

            mdps = [member.name for member in members if member.isfile() and "/gromacs/" in f"/{member.name}" and member.name.endswith(".mdp")]
            equil_mdps = [name for name in mdps if "equilibration" in name.lower()]
            production_mdps = [name for name in mdps if "production" in name.lower()]
            report["mdp_members"] = sorted(mdps)
            required_stage_names = ["step6.0_minimization.mdp"] + [
                f"step6.{index}_equilibration.mdp" for index in range(1, 7)
            ]
            stage_binding: list[dict[str, Any]] = []
            for required_name in required_stage_names:
                matches = [
                    member for member in members
                    if member.isfile() and (
                        member.name == "gromacs/" + required_name
                        or member.name.endswith("/gromacs/" + required_name)
                    )
                ]
                if len(matches) != 1:
                    errors.append(f"expected one archive member for {required_name}, found {len(matches)}")
                    continue
                member = matches[0]
                observed_hash = sha256_member(archive, member)
                if required_name == "step6.0_minimization.mdp":
                    expected_hash = construction.get("minimization_mdp", {}).get("sha256")
                else:
                    expected_hash = construction.get("equilibration_mdp_sha256", {}).get(required_name)
                if observed_hash != expected_hash:
                    errors.append(f"archive {required_name} differs from the manifest hash")
                try:
                    physics = audit_charmm_gui_stage_mdp(
                        required_name, read_member_text(archive, member), float(manifest["global_model"]["temperature_k"])
                    )
                except (ValueError, KeyError) as exc:
                    errors.append(str(exc))
                    physics = None
                stage_binding.append({
                    "name": required_name,
                    "archive_member": member.name,
                    "archive_sha256": observed_hash,
                    "manifest_sha256": expected_hash,
                    "physics": physics,
                })
            report["staged_mdp_archive_binding_and_physics"] = stage_binding
            if len(equil_mdps) != 6:
                errors.append(f"expected exactly six CHARMM-GUI equilibration MDPs, found {len(equil_mdps)}")
            if len(production_mdps) < 1:
                warnings.append("no production MDP was found; it must be reconciled before production")

            gro_member = select_gro(members)
            gro = parse_gro(read_member_text(archive, gro_member))
            report["coordinate_member"] = gro_member.name
            report["coordinate_summary"] = gro
            ligand_count = gro["residue_atom_counts"].get(args.ligand_resname.upper(), 0)
            ligand_molecules = gro["residue_molecule_counts"].get(args.ligand_resname.upper(), 0)
            if ligand_count != args.expected_ligand_atoms or ligand_molecules != 1:
                errors.append(
                    f"expected one {args.ligand_resname} with {args.expected_ligand_atoms} atoms, "
                    f"found {ligand_molecules} molecule(s) and {ligand_count} atoms"
                )

            text_candidates = [
                member
                for member in members
                if member.isfile()
                and "/gromacs/" in f"/{member.name}"
                and member.name.lower().endswith((".itp", ".top"))
                and member.size <= TEXT_LIMIT
            ]
            charge_candidates: list[dict[str, Any]] = []
            hydrogen_masses: list[float] = []
            for member in text_candidates:
                member_text = read_member_text(archive, member)
                charges = parse_itp_charges(member_text, args.ligand_resname)
                hydrogen_masses.extend(parse_hydrogen_masses(member_text))
                if charges:
                    charge_candidates.append(
                        {
                            "member": member.name,
                            "atom_count": len(charges),
                            "charge_sum_e": sum(charges),
                            "all_zero": not any(abs(value) > 1e-8 for value in charges),
                        }
                    )
            report["ligand_charge_candidates"] = charge_candidates
            report["explicit_hydrogen_mass_count"] = len(hydrogen_masses)
            report["maximum_explicit_hydrogen_mass_u"] = max(hydrogen_masses) if hydrogen_masses else None
            report["hydrogen_mass_repartitioning_detected"] = any(value > 1.5 for value in hydrogen_masses)
            if not hydrogen_masses:
                errors.append("no explicit hydrogen masses were found; HMR could not be excluded")
            elif report["hydrogen_mass_repartitioning_detected"]:
                errors.append("hydrogen mass >1.5 u detected; no-HMR contract failed")
            if not charge_candidates:
                errors.append("no O6U [ atoms ] topology section was found")
            else:
                chosen = max(charge_candidates, key=lambda item: item["atom_count"])
                report["ligand_charge_topology"] = chosen
                if chosen["atom_count"] != args.expected_ligand_atoms:
                    errors.append(
                        f"ligand topology contains {chosen['atom_count']} atoms; expected {args.expected_ligand_atoms}"
                    )
                if chosen["all_zero"]:
                    errors.append("all ligand topology partial charges are zero")
                if abs(chosen["charge_sum_e"] - args.expected_ligand_charge) > 0.0001:
                    errors.append(
                        f"ligand topology charge {chosen['charge_sum_e']:.6f} does not match expected "
                        f"{args.expected_ligand_charge:.6f} within 0.0001 e"
                    )

            member_by_name = {member.name: member for member in members if member.isfile()}
            artifact_audit: list[dict[str, str]] = []
            for item in args.expected_artifact:
                if "=" not in item:
                    errors.append(f"invalid --expected-artifact value: {item}")
                    continue
                member_name, expected_hash = item.rsplit("=", 1)
                expected_hash = expected_hash.lower()
                member = member_by_name.get(member_name)
                if member is None:
                    errors.append(f"approved artifact is absent from archive: {member_name}")
                    continue
                observed_hash = sha256_member(archive, member)
                artifact_audit.append(
                    {"member": member_name, "expected_sha256": expected_hash, "observed_sha256": observed_hash}
                )
                if not re.fullmatch(r"[0-9a-f]{64}", expected_hash) or observed_hash != expected_hash:
                    errors.append(f"approved artifact hash mismatch: {member_name}")
            report["approved_artifact_hashes"] = artifact_audit
            if not args.expected_artifact:
                warnings.append(
                    "no --expected-artifact hashes were supplied; approved O6U topology/parameter identity is unaudited"
                )

            counts = gro["residue_atom_counts"]
            molecule_counts = gro["residue_molecule_counts"]
            report["composition_atom_counts"] = {
                key: counts.get(key, 0)
                for key in ("POPC", "PC1", "DSPC", "CHL1", "CHOL", "CLR", "NAG", "BMA", "SOD", "NA", "POT", "K", "CLA", "CL", "TIP3", "SOL")
            }
            report["composition_molecule_counts"] = {
                key: molecule_counts.get(key, 0)
                for key in ("POPC", "PC1", "DSPC", "CHL1", "CHOL", "CLR", "NAG", "BMA", "SOD", "NA", "POT", "K", "CLA", "CL", "TIP3", "SOL")
            }
            cholesterol_count = sum(molecule_counts.get(key, 0) for key in ("CHL1", "CHOL", "CLR"))
            pc1_count = sum(molecule_counts.get(key, 0) for key in ("PC1", "DSPC"))
            if molecule_counts.get("POPC", 0) <= 0:
                errors.append("bulk POPC is absent")
            if cholesterol_count != 3:
                errors.append(f"expected exactly 3 retained cholesterol molecules and no bulk cholesterol; found {cholesterol_count}")
            if pc1_count != 2:
                errors.append(f"expected exactly 2 retained PC1/DSPC molecules; found {pc1_count}")
            if molecule_counts.get("NAG", 0) != 18 or molecule_counts.get("BMA", 0) != 3:
                errors.append(
                    f"resolved glycan count mismatch: NAG={molecule_counts.get('NAG', 0)}, "
                    f"BMA={molecule_counts.get('BMA', 0)}"
                )
            sodium_count = molecule_counts.get("SOD", 0) + molecule_counts.get("NA", 0)
            chloride_count = molecule_counts.get("CLA", 0) + molecule_counts.get("CL", 0)
            potassium_count = molecule_counts.get("POT", 0) + molecule_counts.get("K", 0)
            if sodium_count <= 0 or chloride_count <= 0:
                errors.append("expected NaCl ions were not detected")
            if potassium_count != 0:
                errors.append(f"unexpected potassium ions detected: {potassium_count}")
            forcefield_markers = [
                name for name in names if any(marker in name.lower() for marker in ("forcefield.itp", "charmm36", "cgenff"))
            ]
            report["forcefield_marker_members"] = sorted(forcefield_markers)[:100]
            if not forcefield_markers:
                warnings.append("force-field release could not be identified from archive member names")
            if not any("cgenff" in name.lower() for name in forcefield_markers):
                warnings.append("CGenFF release marker could not be identified from archive member names")
    except (OSError, ValueError, tarfile.TarError) as exc:
        errors.append(str(exc))

    report["warnings"] = warnings
    if args.strict and warnings:
        errors.extend(f"strict warning: {warning}" for warning in warnings)
    report["errors"] = errors
    report["status"] = "pass" if not errors else "fail"
    report["integrity"] = {"payload_sha256": "UNSEALED"}
    report["integrity"]["payload_sha256"] = report_payload_sha256(report, ("integrity", "payload_sha256"))
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": str(output), "errors": len(errors), "warnings": len(warnings)}, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
