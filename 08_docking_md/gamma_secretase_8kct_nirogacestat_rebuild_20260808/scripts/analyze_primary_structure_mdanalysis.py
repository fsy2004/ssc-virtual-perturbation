#!/usr/bin/env python3
"""Primary unsmoothed 8KCT/O6U structural analysis with frozen atom maps.

The input is the full 0-500 ns centered-system trajectory for each of rep01,
rep02, and rep03.  The only inferentially relevant descriptive window is
200-500 ns.  No frame is smoothed, deleted, interpolated, or pooled across
realizations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import MDAnalysis as mda
import numpy as np

from primary_postprocessing_common import (
    PRIMARY_WINDOW_NS,
    REALIZATION_IDS,
    STRUCTURAL_STATIONARITY_METRICS,
    ContractError,
    apply_transform,
    atomic_write_csv,
    atomic_write_json,
    block_diagnostics,
    check_mdanalysis_version,
    continuous_true_events,
    has_placeholder,
    kabsch_transform,
    load_json,
    mapped_indices,
    primary_window_mask,
    require,
    resolve_record,
    robust_first_difference,
    rmsd_nm,
    sha256_file,
    stationarity_diagnostics,
    validate_primary_manifest,
    validate_time_axis,
    verify_atom,
)


RAW_BASE_FIELDS = [
    "system_id",
    "realization_id",
    "frame_index_zero_based",
    "time_ns",
    "in_primary_window_200_500_ns",
    "pocket_aligned_o6u_heavy_rmsd_nm",
    "pocket_aligned_o6u_com_displacement_nm",
    "tm_core_ca_rmsd_nm",
    "protein_ca_rmsd_nm",
    "native_contact_fraction",
]


def topology_identity_sha256(universe: mda.Universe) -> str:
    records = []
    for atom in universe.atoms:
        record = {"index": int(atom.index), "name": str(atom.name), "resname": str(atom.resname), "resid": int(atom.resid)}
        for field in ("segid", "chainID"):
            try:
                value = getattr(atom, field)
            except Exception:
                continue
            if value not in (None, ""):
                record[field] = str(value)
        records.append(record)
    payload = json.dumps(records, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def atom_is_hydrogen(atom: Any) -> bool:
    """Identify hydrogens by chemical identity, independent of HMR mass."""
    try:
        element = str(atom.element).strip().upper()
    except Exception:
        element = ""
    normalized_name = str(atom.name).strip().upper().lstrip("0123456789")
    return element == "H" or normalized_name.startswith("H")


def _endpoint_indices(reference: mda.Universe, trajectory: mda.Universe, endpoint: Mapping[str, Any], label: str) -> tuple[int, int]:
    reference_indices, trajectory_indices = mapped_indices(reference, trajectory, [endpoint], label)
    return int(reference_indices[0]), int(trajectory_indices[0])


def _validate_mapping_record(mapping: Mapping[str, Any], reference: mda.Universe, trajectory: mda.Universe) -> dict[str, Any]:
    require(mapping.get("schema_version") == "1.0", "Structural mapping schema_version must be 1.0")
    require(mapping.get("approval_status") in {"approved", "synthetic_self_test"}, "Structural mapping is not approved")
    require(not has_placeholder(mapping), "Structural mapping contains TODO/REPLACE_ME placeholders")
    require(mapping.get("system_id") == "8kct_nirogacestat_native", "Structural mapping system_id differs")
    require(mapping.get("ligand_resname") == "O6U", "Structural mapping ligand must be O6U")
    require(float(mapping.get("native_contact_cutoff_nm")) == 0.45, "Native-contact cutoff must be 0.45 nm")
    require(float(mapping.get("hydrogen_bond_distance_cutoff_nm")) == 0.35, "Hydrogen-bond D-A cutoff must be 0.35 nm")
    require(float(mapping.get("hydrogen_bond_angular_deviation_cutoff_deg")) == 30.0, "Hydrogen-bond angular-deviation cutoff must be 30 degrees")
    require(mapping.get("coordinate_space") == "pocket_aligned_euclidean_after_validated_whole_cluster_nojump_center", "Structural coordinate-space contract differs")
    observed_fingerprint = topology_identity_sha256(trajectory)
    require(mapping.get("trajectory_atom_identity_sha256") == observed_fingerprint, "Trajectory atom-identity fingerprint differs from the frozen mapping")

    atom_mappings = mapping.get("atom_mappings", {})
    pocket_ref, pocket_traj = mapped_indices(reference, trajectory, atom_mappings.get("pocket_alignment", []), "pocket_alignment")
    tm_ref, tm_traj = mapped_indices(reference, trajectory, atom_mappings.get("tm_core_ca", []), "tm_core_ca")
    ligand_ref, ligand_traj = mapped_indices(reference, trajectory, atom_mappings.get("o6u_heavy", []), "o6u_heavy")
    protein_ref, protein_traj = mapped_indices(reference, trajectory, atom_mappings.get("protein_ca", []), "protein_ca")
    require(len(pocket_ref) >= 3 and len(tm_ref) >= 3, "Alignment mappings require at least three atoms")
    require(len(ligand_ref) >= 2, "O6U heavy-atom mapping is unexpectedly small")
    require(len(protein_ref) >= len(tm_ref), "protein_ca mapping is smaller than tm_core_ca mapping")
    reference_o6u = np.flatnonzero(np.asarray(reference.atoms.resnames) == "O6U")
    trajectory_o6u = np.flatnonzero(np.asarray(trajectory.atoms.resnames) == "O6U")
    require(len(reference_o6u) > 0 and len(trajectory_o6u) > 0, "O6U is absent from the reference or trajectory topology")
    reference_o6u_masses = np.asarray(reference.atoms[reference_o6u].masses, dtype=np.float64)
    trajectory_o6u_masses = np.asarray(trajectory.atoms[trajectory_o6u].masses, dtype=np.float64)
    require(np.all(np.isfinite(reference_o6u_masses)) and np.all(reference_o6u_masses > 0.0), "Reference O6U atoms lack positive finite masses")
    require(np.all(np.isfinite(trajectory_o6u_masses)) and np.all(trajectory_o6u_masses > 0.0), "Trajectory O6U atoms lack positive finite masses")
    reference_o6u_heavy = {int(index) for index in reference_o6u if not atom_is_hydrogen(reference.atoms[index])}
    trajectory_o6u_heavy = {int(index) for index in trajectory_o6u if not atom_is_hydrogen(trajectory.atoms[index])}
    require(set(int(index) for index in ligand_ref) == reference_o6u_heavy, "o6u_heavy mapping is not exhaustive for the reference O6U heavy atoms")
    require(set(int(index) for index in ligand_traj) == trajectory_o6u_heavy, "o6u_heavy mapping is not exhaustive for the trajectory O6U heavy atoms")
    reference_protein_ca = set(int(index) for index in reference.select_atoms("protein and name CA").indices)
    require(reference_protein_ca and set(int(index) for index in protein_ref) == reference_protein_ca, "protein_ca mapping must cover every resolved reference protein CA exactly once")
    ligand_masses = np.asarray(trajectory.atoms[ligand_traj].masses, dtype=np.float64)
    require(np.all(np.isfinite(ligand_masses)) and np.all(ligand_masses > 0.0), "O6U heavy-atom mapping lacks a valid mass")
    require(all(not atom_is_hydrogen(trajectory.atoms[index]) for index in ligand_traj), "O6U heavy-atom mapping contains a hydrogen")

    native_contacts = []
    reference_tolerance = float(mapping.get("reference_distance_tolerance_nm"))
    seen_contact_ids: set[str] = set()
    for number, contact in enumerate(mapping.get("native_contacts", [])):
        contact_id = str(contact.get("contact_id", ""))
        require(contact_id and contact_id not in seen_contact_ids, f"Duplicate/empty native contact ID at row {number}")
        seen_contact_ids.add(contact_id)
        protein_ref_index, protein_traj_index = _endpoint_indices(reference, trajectory, contact["protein_atom"], f"native_contacts[{number}].protein")
        ligand_ref_index, ligand_traj_index = _endpoint_indices(reference, trajectory, contact["ligand_atom"], f"native_contacts[{number}].ligand")
        observed_reference_nm = float(np.linalg.norm(reference.atoms[protein_ref_index].position - reference.atoms[ligand_ref_index].position) / 10.0)
        frozen_reference_nm = float(contact["reference_distance_nm"])
        cutoff_nm = float(contact["cutoff_nm"])
        require(cutoff_nm == 0.45, f"{contact_id} cutoff must be 0.45 nm")
        require(abs(observed_reference_nm - frozen_reference_nm) <= reference_tolerance, f"{contact_id} frozen native distance does not match 8KCT")
        require(observed_reference_nm <= cutoff_nm + reference_tolerance, f"{contact_id} is not a native <=0.45 nm contact")
        native_contacts.append({**contact, "contact_id": contact_id, "protein_ref_index": protein_ref_index, "ligand_ref_index": ligand_ref_index, "protein_traj_index": protein_traj_index, "ligand_traj_index": ligand_traj_index, "observed_reference_distance_nm": observed_reference_nm})
    require(native_contacts, "At least one explicit native contact is required")
    protein_atoms = reference.select_atoms("protein")
    protein_masses = np.asarray(protein_atoms.masses, dtype=np.float64)
    require(np.all(np.isfinite(protein_masses)) and np.all(protein_masses > 0.0), "Reference protein atoms lack positive finite masses")
    protein_heavy_indices = np.asarray(
        [int(index) for index in protein_atoms.indices if not atom_is_hydrogen(reference.atoms[index])],
        dtype=np.int64,
    )
    expected_native_pairs: set[tuple[int, int]] = set()
    cutoff_angstrom = float(mapping["native_contact_cutoff_nm"]) * 10.0
    for ligand_index in sorted(reference_o6u_heavy):
        deltas = reference.atoms[protein_heavy_indices].positions - reference.atoms[ligand_index].position
        for local_index in np.flatnonzero(np.linalg.norm(deltas, axis=1) <= cutoff_angstrom + reference_tolerance * 10.0):
            expected_native_pairs.add((int(protein_heavy_indices[local_index]), int(ligand_index)))
    observed_native_pairs = {(int(item["protein_ref_index"]), int(item["ligand_ref_index"])) for item in native_contacts}
    require(observed_native_pairs == expected_native_pairs, "native_contacts must exhaustively enumerate every reference protein-heavy/O6U-heavy pair within 0.45 nm")

    distances = []
    seen_distance_ids: set[str] = set()
    for number, distance in enumerate(mapping.get("prespecified_distances", [])):
        metric_id = str(distance.get("metric_id", ""))
        require(metric_id and metric_id not in seen_distance_ids, f"Duplicate/empty distance metric ID at row {number}")
        seen_distance_ids.add(metric_id)
        _, first_index = _endpoint_indices(reference, trajectory, distance["atom1"], f"prespecified_distances[{number}].atom1")
        _, second_index = _endpoint_indices(reference, trajectory, distance["atom2"], f"prespecified_distances[{number}].atom2")
        distances.append({**distance, "metric_id": metric_id, "atom1_traj_index": first_index, "atom2_traj_index": second_index})

    hydrogen_bonds = []
    seen_hbond_ids: set[str] = set()
    for number, bond in enumerate(mapping.get("hydrogen_bonds", [])):
        metric_id = str(bond.get("metric_id", ""))
        require(metric_id and metric_id not in seen_hbond_ids, f"Duplicate/empty hydrogen-bond ID at row {number}")
        seen_hbond_ids.add(metric_id)
        donor_ref_index, donor_index = _endpoint_indices(reference, trajectory, bond["donor"], f"hydrogen_bonds[{number}].donor")
        hydrogen_ref_index, hydrogen_index = _endpoint_indices(reference, trajectory, bond["hydrogen"], f"hydrogen_bonds[{number}].hydrogen")
        acceptor_ref_index, acceptor_index = _endpoint_indices(reference, trajectory, bond["acceptor"], f"hydrogen_bonds[{number}].acceptor")
        require(len({donor_index, hydrogen_index, acceptor_index}) == 3, f"{metric_id} donor/hydrogen/acceptor indices must differ")
        endpoint_masses = np.asarray(trajectory.atoms[[donor_index, hydrogen_index, acceptor_index]].masses, dtype=np.float64)
        require(np.all(np.isfinite(endpoint_masses)) and np.all(endpoint_masses > 0.0), f"{metric_id} donor/H/acceptor masses are invalid")
        require(
            not atom_is_hydrogen(trajectory.atoms[donor_index])
            and atom_is_hydrogen(trajectory.atoms[hydrogen_index])
            and not atom_is_hydrogen(trajectory.atoms[acceptor_index]),
            f"{metric_id} donor/H/acceptor identities do not define an explicit hydrogen bond",
        )
        reference_triplet = reference.atoms[[donor_ref_index, hydrogen_ref_index, acceptor_ref_index]].positions
        donor_hydrogen_nm = float(np.linalg.norm(reference_triplet[0] - reference_triplet[1]) / 10.0)
        donor_acceptor_nm = float(np.linalg.norm(reference_triplet[0] - reference_triplet[2]) / 10.0)
        reference_deviation_deg = float(180.0 - _angle_deg(reference_triplet[0], reference_triplet[1], reference_triplet[2]))
        require(donor_hydrogen_nm <= 0.13, f"{metric_id} donor and explicit H are not covalently proximate in the reference")
        require(donor_acceptor_nm <= float(mapping["hydrogen_bond_distance_cutoff_nm"]) and reference_deviation_deg <= float(mapping["hydrogen_bond_angular_deviation_cutoff_deg"]), f"{metric_id} is not a native reference hydrogen-bond geometry")
        hydrogen_bonds.append({**bond, "metric_id": metric_id, "donor_traj_index": donor_index, "hydrogen_traj_index": hydrogen_index, "acceptor_traj_index": acceptor_index})
    require(hydrogen_bonds, "At least one endpoint-resolved hydrogen-bond geometry must be frozen")
    return {
        "pocket_ref": pocket_ref,
        "pocket_traj": pocket_traj,
        "tm_ref": tm_ref,
        "tm_traj": tm_traj,
        "ligand_ref": ligand_ref,
        "ligand_traj": ligand_traj,
        "protein_ref": protein_ref,
        "protein_traj": protein_traj,
        "ligand_masses": ligand_masses,
        "native_contacts": native_contacts,
        "distances": distances,
        "hydrogen_bonds": hydrogen_bonds,
        "topology_identity_sha256": observed_fingerprint,
    }


def _angle_deg(first: np.ndarray, vertex: np.ndarray, third: np.ndarray) -> float:
    vector1 = np.asarray(first, dtype=np.float64) - vertex
    vector2 = np.asarray(third, dtype=np.float64) - vertex
    denominator = float(np.linalg.norm(vector1) * np.linalg.norm(vector2))
    require(denominator > 0.0, "Hydrogen-bond angle has a zero-length vector")
    cosine = float(np.clip(np.dot(vector1, vector2) / denominator, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def _process_realization(
    manifest: Mapping[str, Any],
    manifest_base: Path,
    realization: Mapping[str, Any],
    reference: mda.Universe,
    mapping_record: Mapping[str, Any],
    output_directory: Path,
) -> tuple[dict[str, Any], np.ndarray]:
    realization_id = str(realization["realization_id"])
    topology_path = resolve_record(manifest_base, realization["topology"], f"{realization_id}.topology")
    trajectory_path = resolve_record(manifest_base, realization["centered_system_trajectory"], f"{realization_id}.centered_system_trajectory")
    try:
        universe = mda.Universe(str(topology_path), str(trajectory_path))
    except Exception as exc:
        raise ContractError(f"MDAnalysis could not open {realization_id}: {exc}") from exc
    mapping = _validate_mapping_record(mapping_record, reference, universe)
    times_ns = np.asarray([float(ts.time) / 1000.0 for ts in universe.trajectory], dtype=np.float64)
    tolerance_ns = float(manifest["time_contract"]["endpoint_tolerance_ns"])
    saved_step_ns = validate_time_axis(times_ns, tolerance_ns)
    window_mask = primary_window_mask(times_ns, tolerance_ns)
    reference_positions = np.asarray(reference.atoms.positions, dtype=np.float64).copy()
    require(np.all(np.isfinite(reference_positions)), "Reference coordinates contain NaN or infinity")

    pocket_reference = reference_positions[mapping["pocket_ref"]]
    tm_reference = reference_positions[mapping["tm_ref"]]
    ligand_reference = reference_positions[mapping["ligand_ref"]]
    protein_reference = reference_positions[mapping["protein_ref"]]
    ligand_reference_com = np.average(ligand_reference, axis=0, weights=mapping["ligand_masses"])

    contact_fields: list[str] = []
    for contact in mapping["native_contacts"]:
        contact_fields.extend([f"contact__{contact['contact_id']}__distance_nm", f"contact__{contact['contact_id']}__present"])
    distance_fields = [f"distance__{item['metric_id']}__nm" for item in mapping["distances"]]
    hbond_fields: list[str] = []
    for bond in mapping["hydrogen_bonds"]:
        hbond_fields.extend([
            f"hbond__{bond['metric_id']}__donor_acceptor_nm",
            f"hbond__{bond['metric_id']}__dha_angle_deg",
            f"hbond__{bond['metric_id']}__angular_deviation_deg",
            f"hbond__{bond['metric_id']}__present",
        ])
    raw_fields = RAW_BASE_FIELDS + contact_fields + distance_fields + hbond_fields
    raw_rows: list[dict[str, Any]] = []

    protein_mean = np.zeros((len(mapping["protein_traj"]), 3), dtype=np.float64)
    protein_m2 = np.zeros_like(protein_mean)
    protein_rmsf_count = 0
    for frame_index, ts in enumerate(universe.trajectory):
        coordinates = np.asarray(universe.atoms.positions, dtype=np.float64)
        require(np.all(np.isfinite(coordinates)), f"{realization_id} frame {frame_index} contains nonfinite coordinates")
        dimensions = np.asarray(ts.dimensions, dtype=np.float64)
        require(dimensions.shape == (6,) and np.all(np.isfinite(dimensions)) and np.all(dimensions[:3] > 0.0), f"{realization_id} frame {frame_index} has an invalid box")

        pocket_rotation, pocket_mobile_center, pocket_reference_center = kabsch_transform(coordinates[mapping["pocket_traj"]], pocket_reference)
        tm_rotation, tm_mobile_center, tm_reference_center = kabsch_transform(coordinates[mapping["tm_traj"]], tm_reference)
        aligned_ligand = apply_transform(coordinates[mapping["ligand_traj"]], pocket_rotation, pocket_mobile_center, pocket_reference_center)
        aligned_tm = apply_transform(coordinates[mapping["tm_traj"]], tm_rotation, tm_mobile_center, tm_reference_center)
        aligned_protein = apply_transform(coordinates[mapping["protein_traj"]], tm_rotation, tm_mobile_center, tm_reference_center)
        ligand_com = np.average(aligned_ligand, axis=0, weights=mapping["ligand_masses"])

        row: dict[str, Any] = {
            "system_id": manifest["system_id"],
            "realization_id": realization_id,
            "frame_index_zero_based": frame_index,
            "time_ns": float(times_ns[frame_index]),
            "in_primary_window_200_500_ns": int(bool(window_mask[frame_index])),
            "pocket_aligned_o6u_heavy_rmsd_nm": rmsd_nm(aligned_ligand, ligand_reference),
            "pocket_aligned_o6u_com_displacement_nm": float(np.linalg.norm(ligand_com - ligand_reference_com) / 10.0),
            "tm_core_ca_rmsd_nm": rmsd_nm(aligned_tm, tm_reference),
            "protein_ca_rmsd_nm": rmsd_nm(aligned_protein, protein_reference),
        }

        contact_presence = []
        for contact in mapping["native_contacts"]:
            pair = coordinates[[contact["protein_traj_index"], contact["ligand_traj_index"]]]
            aligned_pair = apply_transform(pair, pocket_rotation, pocket_mobile_center, pocket_reference_center)
            distance_nm = float(np.linalg.norm(aligned_pair[0] - aligned_pair[1]) / 10.0)
            present = int(distance_nm <= float(contact["cutoff_nm"]))
            row[f"contact__{contact['contact_id']}__distance_nm"] = distance_nm
            row[f"contact__{contact['contact_id']}__present"] = present
            contact_presence.append(present)
        row["native_contact_fraction"] = float(np.mean(contact_presence))

        for distance in mapping["distances"]:
            pair = coordinates[[distance["atom1_traj_index"], distance["atom2_traj_index"]]]
            aligned_pair = apply_transform(pair, pocket_rotation, pocket_mobile_center, pocket_reference_center)
            row[f"distance__{distance['metric_id']}__nm"] = float(np.linalg.norm(aligned_pair[0] - aligned_pair[1]) / 10.0)

        for bond in mapping["hydrogen_bonds"]:
            triplet = coordinates[[bond["donor_traj_index"], bond["hydrogen_traj_index"], bond["acceptor_traj_index"]]]
            aligned_triplet = apply_transform(triplet, pocket_rotation, pocket_mobile_center, pocket_reference_center)
            donor_acceptor_nm = float(np.linalg.norm(aligned_triplet[0] - aligned_triplet[2]) / 10.0)
            angle = _angle_deg(aligned_triplet[0], aligned_triplet[1], aligned_triplet[2])
            deviation = float(180.0 - angle)
            present = int(donor_acceptor_nm <= float(mapping_record["hydrogen_bond_distance_cutoff_nm"]) and deviation <= float(mapping_record["hydrogen_bond_angular_deviation_cutoff_deg"]))
            row[f"hbond__{bond['metric_id']}__donor_acceptor_nm"] = donor_acceptor_nm
            row[f"hbond__{bond['metric_id']}__dha_angle_deg"] = angle
            row[f"hbond__{bond['metric_id']}__angular_deviation_deg"] = deviation
            row[f"hbond__{bond['metric_id']}__present"] = present

        require(all(math.isfinite(float(value)) for key, value in row.items() if key not in {"system_id", "realization_id"}), f"{realization_id} frame {frame_index} produced a nonfinite metric")
        raw_rows.append(row)
        if window_mask[frame_index]:
            protein_rmsf_count += 1
            delta = aligned_protein - protein_mean
            protein_mean += delta / protein_rmsf_count
            protein_m2 += delta * (aligned_protein - protein_mean)

    require(len(raw_rows) == len(times_ns), f"{realization_id} lost frames during analysis")
    require(protein_rmsf_count == int(np.count_nonzero(window_mask)), f"{realization_id} RMSF frame accounting differs")
    rmsf_nm = np.sqrt(np.sum(protein_m2 / protein_rmsf_count, axis=1)) / 10.0
    require(np.all(np.isfinite(rmsf_nm)), f"{realization_id} RMSF contains nonfinite values")

    realization_directory = output_directory / realization_id
    raw_path = realization_directory / "structural_raw_unsmoothed.csv"
    raw_count = atomic_write_csv(raw_path, raw_fields, raw_rows)
    rmsf_rows = []
    for mapping_entry, value in zip(mapping_record["atom_mappings"]["protein_ca"], rmsf_nm, strict=True):
        atom = mapping_entry["trajectory"]
        rmsf_rows.append({
            "system_id": manifest["system_id"],
            "realization_id": realization_id,
            "window_start_ns": PRIMARY_WINDOW_NS[0],
            "window_end_ns": PRIMARY_WINDOW_NS[1],
            "trajectory_atom_index_zero_based": int(atom["index"]),
            "segid": atom.get("segid", ""),
            "resid": int(atom["resid"]),
            "resname": atom["resname"],
            "atom_name": atom["name"],
            "rmsf_nm": float(value),
        })
    rmsf_path = realization_directory / "protein_ca_rmsf_200_500ns.csv"
    atomic_write_csv(rmsf_path, list(rmsf_rows[0]), rmsf_rows)

    scalar_metrics = [field for field in raw_fields if field not in {"system_id", "realization_id", "frame_index_zero_based", "time_ns", "in_primary_window_200_500_ns"}]
    diagnostics_payload: dict[str, Any] = {}
    block_rows: list[dict[str, Any]] = []
    review_rows: list[dict[str, Any]] = []
    for metric in scalar_metrics:
        values = np.asarray([float(row[metric]) for row in raw_rows], dtype=np.float64)
        diagnostics_payload[metric] = {}
        for window_name, mask in (("full_0_500_ns", np.ones(len(times_ns), dtype=bool)), ("primary_200_500_ns", window_mask)):
            diagnostic, current_blocks = block_diagnostics(values[mask], times_ns[mask], manifest["diagnostics"])
            if window_name == "primary_200_500_ns" and metric in STRUCTURAL_STATIONARITY_METRICS:
                diagnostic["stationarity"] = stationarity_diagnostics(
                    values[mask],
                    times_ns[mask],
                    manifest["diagnostics"]["stationarity"],
                    float(manifest["acceptance_gates"]["stationarity_scale_floors"]["structural"][metric]),
                )
            diagnostics_payload[metric][window_name] = diagnostic
            for block in current_blocks:
                block_rows.append({"system_id": manifest["system_id"], "realization_id": realization_id, "metric": metric, "window": window_name, **block})
        if metric in STRUCTURAL_STATIONARITY_METRICS:
            first_difference_summary, flags = robust_first_difference(
                values,
                times_ns,
                float(manifest["diagnostics"]["robust_first_difference_z_threshold"]),
            )
            diagnostics_payload[metric]["first_difference_review"] = first_difference_summary
            for flag in flags:
                review_rows.append({
                    "system_id": manifest["system_id"],
                    "realization_id": realization_id,
                    "metric": metric,
                    **flag,
                })
    block_path = realization_directory / "structural_block_summaries.csv"
    block_fields = ["system_id", "realization_id", "metric", "window", "block_index_zero_based", "start_time_ns", "end_time_ns", "frame_count", "mean", "median", "minimum", "maximum"]
    atomic_write_csv(block_path, block_fields, block_rows)
    review_path = realization_directory / "structural_first_difference_review_flags.csv"
    review_fields = [
        "system_id", "realization_id", "metric", "row_index_zero_based",
        "time_before_ns", "time_after_ns", "value_before", "value_after",
        "first_difference", "median_first_difference", "mad_first_difference",
        "robust_z", "method", "review_required", "point_retained",
    ]
    atomic_write_csv(review_path, review_fields, review_rows)

    occupancy_rows = []
    for contact in mapping["native_contacts"]:
        field = f"contact__{contact['contact_id']}__present"
        values = np.asarray([int(row[field]) for row in raw_rows], dtype=np.int64)
        for window_name, mask in (("full_0_500_ns", np.ones(len(times_ns), dtype=bool)), ("primary_200_500_ns", window_mask)):
            occupancy_rows.append({
                "system_id": manifest["system_id"],
                "realization_id": realization_id,
                "window": window_name,
                "contact_id": contact["contact_id"],
                "frames": int(np.count_nonzero(mask)),
                "occupancy": float(np.mean(values[mask])),
                "cutoff_nm": float(contact["cutoff_nm"]),
                "reference_distance_nm": float(contact["observed_reference_distance_nm"]),
            })
    occupancy_path = realization_directory / "native_contact_occupancy.csv"
    atomic_write_csv(occupancy_path, list(occupancy_rows[0]), occupancy_rows)

    mandatory_primary = [
        "pocket_aligned_o6u_heavy_rmsd_nm",
        "pocket_aligned_o6u_com_displacement_nm",
        "tm_core_ca_rmsd_nm",
        "protein_ca_rmsd_nm",
        "native_contact_fraction",
    ]
    sampling_failures = []
    stationarity_failures = []
    for metric in mandatory_primary:
        primary_diagnostic = diagnostics_payload[metric]["primary_200_500_ns"]
        if primary_diagnostic["status"] != "pass":
            sampling_failures.append(metric)
        if primary_diagnostic["stationarity"]["status"] != "pass":
            stationarity_failures.append(metric)
    pose_gates = manifest["acceptance_gates"]["native_pose"]
    primary_rows = [row for row, selected in zip(raw_rows, window_mask, strict=True) if selected]
    joint_pose_pass = np.asarray([
        float(row["pocket_aligned_o6u_heavy_rmsd_nm"]) <= float(pose_gates["maximum_pocket_aligned_o6u_heavy_rmsd_nm"])
        and float(row["pocket_aligned_o6u_com_displacement_nm"]) <= float(pose_gates["maximum_o6u_com_displacement_nm"])
        and float(row["native_contact_fraction"]) >= float(pose_gates["minimum_native_contact_fraction"])
        for row in primary_rows
    ], dtype=bool)
    joint_pose_fraction = float(np.mean(joint_pose_pass))
    scientific_failures = []
    if joint_pose_fraction < float(pose_gates["minimum_fraction_of_primary_frames_meeting_all_pose_gates"]):
        scientific_failures.append("primary_native_pose_fraction_gate_failed")
    full_ligand_rmsd = np.asarray([float(row["pocket_aligned_o6u_heavy_rmsd_nm"]) for row in raw_rows], dtype=np.float64)
    full_ligand_com = np.asarray([float(row["pocket_aligned_o6u_com_displacement_nm"]) for row in raw_rows], dtype=np.float64)
    full_native_fraction = np.asarray([float(row["native_contact_fraction"]) for row in raw_rows], dtype=np.float64)
    minimum_event_duration_ns = float(pose_gates["minimum_continuous_event_duration_ns"])
    event_definitions = (
        ("o6u_heavy_rmsd_egress", full_ligand_rmsd > float(pose_gates["maximum_pocket_aligned_o6u_heavy_rmsd_nm"])),
        ("o6u_com_displacement_egress", full_ligand_com > float(pose_gates["maximum_o6u_com_displacement_nm"])),
        ("native_contact_loss", full_native_fraction < float(pose_gates["minimum_native_contact_fraction"])),
    )
    ligand_events: list[dict[str, Any]] = []
    for event_type, condition in event_definitions:
        ligand_events.extend(continuous_true_events(condition, times_ns, minimum_event_duration_ns, event_type))
    ligand_events.sort(key=lambda item: (float(item["start_time_ns"]), str(item["event_type"])))
    if ligand_events:
        scientific_failures.append("prespecified_continuous_ligand_egress_or_contact_loss_event")
    event_path = realization_directory / "native_pose_continuous_events.csv"
    event_fields = [
        "system_id", "realization_id", "event_type", "start_frame_index_zero_based", "end_frame_index_zero_based",
        "start_time_ns", "end_time_ns", "continuous_duration_ns", "frame_count", "minimum_duration_ns",
        "gap_bridging", "all_frames_retained",
    ]
    atomic_write_csv(
        event_path,
        event_fields,
        ({"system_id": manifest["system_id"], "realization_id": realization_id, **event} for event in ligand_events),
    )
    pose_gate_summary = {
        "window_ns": list(PRIMARY_WINDOW_NS),
        "criteria": dict(pose_gates),
        "primary_frame_count": len(primary_rows),
        "frames_meeting_all_pose_gates": int(np.count_nonzero(joint_pose_pass)),
        "fraction_of_primary_frames_meeting_all_pose_gates": joint_pose_fraction,
        "continuous_event_search_window_ns": [0.0, 500.0],
        "continuous_event_rule": pose_gates["continuous_event_rule"],
        "minimum_continuous_event_duration_ns": minimum_event_duration_ns,
        "qualifying_continuous_event_count": len(ligand_events),
        "qualifying_continuous_events": ligand_events,
        "status": "pass" if not scientific_failures else "scientific_fail",
        "failure_triggers_rerun_or_extension": False,
    }
    summary = {
        "schema_version": "1.0",
        "analysis": "primary_structural_mdanalysis",
        "analysis_role": "single_system_native_pose_geometric_compatibility_only",
        "system_id": manifest["system_id"],
        "realization_id": realization_id,
        "construction_count": 1,
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "input_frame_count": len(times_ns),
        "primary_window_frame_count": int(np.count_nonzero(window_mask)),
        "saved_step_ns": saved_step_ns,
        "first_time_ns": float(times_ns[0]),
        "last_time_ns": float(times_ns[-1]),
        "technical_status": "pass",
        "sampling_status": "pass" if not sampling_failures and not stationarity_failures else "inconclusive",
        "sampling_failures": sampling_failures,
        "stationarity_status": "pass" if not stationarity_failures else "fail",
        "stationarity_failures": stationarity_failures,
        "scientific_status": "pass" if not scientific_failures else "fail",
        "scientific_failures": scientific_failures,
        "first_difference_review_flag_count": len(review_rows),
        "first_difference_flags_remove_points": False,
        "native_pose_acceptance": pose_gate_summary,
        "scientific_failure_triggers_rerun_or_extension": False,
        "mdanalysis_version": mda.__version__,
        "trajectory_atom_identity_sha256": mapping["topology_identity_sha256"],
        "input_files": {
            "topology": {"path": str(topology_path), "sha256": sha256_file(topology_path)},
            "centered_system_trajectory": {"path": str(trajectory_path), "sha256": sha256_file(trajectory_path)},
        },
        "data_handling": {"raw_frames_retained": True, "smoothing": False, "frame_deletion": False, "interpolation": False, "realization_pooling": False},
        "diagnostics": diagnostics_payload,
        "outputs": {},
    }
    summary_path = realization_directory / "structural_summary.json"
    for path in (raw_path, rmsf_path, block_path, review_path, occupancy_path, event_path):
        summary["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    atomic_write_json(summary_path, summary)
    return summary, times_ns


def run(manifest_path: Path, output_root: Path, allow_synthetic: bool = False) -> dict[str, Any]:
    check_mdanalysis_version(mda.__version__)
    manifest, manifest_base = validate_primary_manifest(manifest_path, allow_synthetic=allow_synthetic)
    mapping_path = resolve_record(manifest_base, manifest["mapping_records"]["structural"], "mapping_records.structural")
    reference_topology = resolve_record(manifest_base, manifest["reference"]["topology"], "reference.topology")
    reference_coordinates = resolve_record(manifest_base, manifest["reference"]["coordinates"], "reference.coordinates")
    mapping_record = load_json(mapping_path)
    try:
        reference = mda.Universe(str(reference_topology), str(reference_coordinates))
    except Exception as exc:
        raise ContractError(f"MDAnalysis could not open the frozen 8KCT reference: {exc}") from exc
    output_directory = output_root.resolve() / "structural_analysis"
    require(not output_directory.exists(), f"Refusing to overwrite an existing structural output directory: {output_directory}")
    output_directory.mkdir(parents=True)
    summaries = []
    shared_times: np.ndarray | None = None
    for realization in manifest["realizations"]:
        summary, times = _process_realization(manifest, manifest_base, realization, reference, mapping_record, output_directory)
        if shared_times is None:
            shared_times = times
        else:
            require(np.allclose(times, shared_times, rtol=0.0, atol=float(manifest["time_contract"]["endpoint_tolerance_ns"])), "Saved frame times differ among realizations")
        summaries.append(summary)
    require([summary["realization_id"] for summary in summaries] == list(REALIZATION_IDS), "Structural analysis lost or reordered a realization")
    overall_status = "pass" if all(summary["sampling_status"] == "pass" and summary["scientific_status"] == "pass" for summary in summaries) else "inconclusive"
    complete = {
        "schema_version": "1.0",
        "status": overall_status,
        "technical_status": "pass",
        "system_id": manifest["system_id"],
        "construction_count": 1,
        "realization_ids": list(REALIZATION_IDS),
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": list(PRIMARY_WINDOW_NS),
        "extension_or_recovery_window": False,
        "mdanalysis_version": mda.__version__,
        "manifest": {"path": str(manifest_path.resolve()), "sha256": sha256_file(manifest_path.resolve())},
        "mapping_record": {"path": str(mapping_path), "sha256": sha256_file(mapping_path)},
        "realization_summaries": [
            {"realization_id": item["realization_id"], "technical_status": item["technical_status"], "sampling_status": item["sampling_status"], "stationarity_status": item["stationarity_status"], "scientific_status": item["scientific_status"], "scientific_failures": item["scientific_failures"]}
            for item in summaries
        ],
        "prohibited_outputs": {"mmgbsa": False, "mmpbsa": False, "smoothed_trace": False, "deleted_frames": False, "interpolated_frames": False},
    }
    atomic_write_json(output_directory / "COMPLETE.json", complete)
    return complete


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True, help="Approved primary postprocessing manifest")
    parser.add_argument("--output-root", type=Path, required=True, help="New output root; existing structural_analysis is never overwritten")
    parser.add_argument("--allow-synthetic", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = run(args.manifest, args.output_root, allow_synthetic=args.allow_synthetic)
    except (ContractError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 3


if __name__ == "__main__":
    raise SystemExit(main())
