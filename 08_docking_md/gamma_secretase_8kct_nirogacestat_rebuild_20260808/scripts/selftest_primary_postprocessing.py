#!/usr/bin/env python3
"""End-to-end synthetic self-test for the primary postprocessing bundle."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np

import analyze_membrane_qc_mdanalysis
import analyze_primary_structure_mdanalysis
import gmx_energy_qc
import validate_primary_postprocessing
from primary_postprocessing_common import (
    ContractError,
    atom_identity,
    atomic_write_csv,
    atomic_write_json,
    continuous_true_events,
    require,
    sha256_file,
    stationarity_diagnostics,
)


NAMES = ["CA", "CA", "CA", "CA", "NZ", "HZ1", "O", "N1", "C1", "O1", "P", "P", "P", "P", "O", "O"]
ELEMENTS = ["C", "C", "C", "C", "N", "H", "O", "N", "C", "O", "P", "P", "P", "P", "O", "O"]
RESNAMES = ["ALA", "LEU", "VAL", "ILE", "LYS", "LYS", "LEU", "O6U", "O6U", "O6U", "POPC", "POPC", "POPC", "POPC", "TIP3", "TIP3"]
RESIDS = [10, 20, 30, 40, 380, 380, 432, 900, 900, 900, 1001, 1002, 1003, 1004, 2001, 2002]
SEGMENT_FOR_ATOM = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3]
SEGIDS = ["PROA", "LIGA", "MEMB", "SOLV"]


def create_universe() -> mda.Universe:
    atom_resindex = np.arange(len(NAMES), dtype=np.int64)
    residue_segindex = np.asarray(SEGMENT_FOR_ATOM, dtype=np.int64)
    universe = mda.Universe.empty(
        len(NAMES),
        n_residues=len(NAMES),
        n_segments=len(SEGIDS),
        atom_resindex=atom_resindex,
        residue_segindex=residue_segindex,
        trajectory=True,
    )
    universe.add_TopologyAttr("names", NAMES)
    universe.add_TopologyAttr("types", ELEMENTS)
    universe.add_TopologyAttr("elements", ELEMENTS)
    universe.add_TopologyAttr("resnames", RESNAMES)
    universe.add_TopologyAttr("resids", RESIDS)
    universe.add_TopologyAttr("segids", SEGIDS)
    universe.add_TopologyAttr("chainIDs", ["A"] * 7 + ["L"] * 3 + ["M"] * 4 + ["W"] * 2)
    universe.add_TopologyAttr("masses", [12.011, 12.011, 12.011, 12.011, 14.007, 1.008, 15.999, 14.007, 12.011, 15.999, 30.974, 30.974, 30.974, 30.974, 15.999, 15.999])
    return universe


def reference_coordinates() -> np.ndarray:
    return np.asarray(
        [
            [-5.0, 0.0, 10.0],
            [5.0, 0.0, 10.0],
            [-5.0, 0.0, -10.0],
            [5.0, 0.0, -10.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 0.0, 1.0],
            [-15.0, -10.0, 20.0],
            [15.0, 10.0, 20.0],
            [-15.0, 10.0, -20.0],
            [15.0, -10.0, -20.0],
            [20.0, 20.0, 30.0],
            [1.0, 1.0, 0.0],
        ],
        dtype=np.float32,
    )


def write_topology_and_trajectories(root: Path) -> tuple[Path, list[Path]]:
    universe = create_universe()
    coordinates = reference_coordinates()
    universe.atoms.positions = coordinates
    universe.dimensions = np.asarray([60.0, 60.0, 100.0, 90.0, 90.0, 90.0], dtype=np.float32)
    topology = root / "synthetic_reference_and_topology.pdb"
    with mda.Writer(str(topology), n_atoms=len(universe.atoms)) as writer:
        writer.write(universe.atoms)
    trajectories = []
    for realization_index, realization_id in enumerate(("rep01", "rep02", "rep03"), start=1):
        trajectory = root / f"{realization_id}.dcd"
        rng = np.random.default_rng(1000 + realization_index)
        with mda.Writer(str(trajectory), n_atoms=len(universe.atoms), dt=1000.0) as writer:
            for frame in range(501):
                current = coordinates.astype(np.float64).copy()
                angle = float(rng.normal(0.0, 0.002))
                rotation = np.asarray([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
                current[:10] = current[:10] @ rotation.T
                current[:10] += np.asarray([rng.normal(0.0, 0.02), rng.normal(0.0, 0.02), 0.0])
                current[7:10] += rng.normal(0.0, 0.01, size=(3, 3))
                current[10:14, 2] += rng.normal(0.0, 0.03, size=4)
                universe.atoms.positions = current.astype(np.float32)
                universe.dimensions = np.asarray([60.0, 60.0, 100.0, 90.0, 90.0, 90.0], dtype=np.float32)
                universe.trajectory.ts.time = frame * 1000.0
                writer.write(universe.atoms)
        trajectories.append(trajectory)
    return topology, trajectories


def identities(topology: Path) -> tuple[mda.Universe, list[dict[str, Any]], str]:
    universe = mda.Universe(str(topology))
    records = [atom_identity(atom) for atom in universe.atoms]
    payload = json.dumps(records, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    import hashlib

    fingerprint = hashlib.sha256(payload).hexdigest()
    return universe, records, fingerprint


def pair(records: list[dict[str, Any]], index: int) -> dict[str, Any]:
    return {"reference": records[index], "trajectory": records[index]}


def write_mapping_records(root: Path, topology: Path) -> tuple[Path, Path, Path]:
    universe, records, fingerprint = identities(topology)
    synthetic_native_contacts = []
    for protein_index in (0, 1, 2, 3, 4, 6):
        for ligand_index in (7, 8, 9):
            distance_nm = float(np.linalg.norm(universe.atoms[protein_index].position - universe.atoms[ligand_index].position) / 10.0)
            if distance_nm <= 0.45:
                synthetic_native_contacts.append({
                    "contact_id": f"ref{protein_index}__o6u{ligand_index}",
                    "protein_atom": pair(records, protein_index),
                    "ligand_atom": pair(records, ligand_index),
                    "reference_distance_nm": distance_nm,
                    "cutoff_nm": 0.45,
                })
    structural = {
        "schema_version": "1.0",
        "approval_status": "synthetic_self_test",
        "system_id": "8kct_nirogacestat_native",
        "reference_structure_id": "SYNTHETIC_8KCT",
        "ligand_resname": "O6U",
        "trajectory_atom_identity_sha256": fingerprint,
        "coordinate_space": "pocket_aligned_euclidean_after_validated_whole_cluster_nojump_center",
        "native_contact_cutoff_nm": 0.45,
        "reference_distance_tolerance_nm": 0.0001,
        "hydrogen_bond_distance_cutoff_nm": 0.35,
        "hydrogen_bond_angular_deviation_cutoff_deg": 30.0,
        "atom_mappings": {
            "pocket_alignment": [pair(records, index) for index in (0, 1, 2, 3)],
            "tm_core_ca": [pair(records, index) for index in (0, 1, 2, 3)],
            "o6u_heavy": [pair(records, index) for index in (7, 8, 9)],
            "protein_ca": [pair(records, index) for index in (0, 1, 2, 3)],
        },
        "native_contacts": synthetic_native_contacts,
        "prespecified_distances": [
            {"metric_id": "leu432_o__o6u_c1", "atom1": pair(records, 6), "atom2": pair(records, 8)}
        ],
        "hydrogen_bonds": [
            {"metric_id": "lys380_nz_hz1__o6u_n1", "literature_residue": "synthetic Lys380", "donor": pair(records, 4), "hydrogen": pair(records, 5), "acceptor": pair(records, 7)}
        ],
    }
    structural_path = root / "structural_mapping.json"
    atomic_write_json(structural_path, structural)
    membrane = {
        "schema_version": "1.0",
        "approval_status": "synthetic_self_test",
        "system_id": "8kct_nirogacestat_native",
        "trajectory_atom_identity_sha256": fingerprint,
        "box_geometry_required": "orthorhombic",
        "frozen_atom_groups": {
            "upper_leaflet_phosphate_atoms": [records[10], records[11]],
            "lower_leaflet_phosphate_atoms": [records[12], records[13]],
            "protein_tilt_upper_anchor_atoms": [records[0], records[1]],
            "protein_tilt_lower_anchor_atoms": [records[2], records[3]],
            "protein_heavy_atoms": [records[index] for index in (0, 1, 2, 3, 4, 6)],
            "water_oxygen_atoms": [records[14], records[15]],
        },
        "metric_settings": {
            "phosphate_density_bandwidth_nm": 0.1,
            "phosphate_density_grid_nm": [-3.0, 3.0, 301],
            "leaflet_hysteresis_nm": 0.5,
            "hydrophobic_core_half_thickness_nm": 1.2,
            "protein_xy_exclusion_nm": 0.5,
            "water_cluster_cutoff_nm": 0.35,
            "orthorhombic_angle_tolerance_deg": 0.01,
        },
        "qc_gates": {
            "maximum_cumulative_leaflet_flip_events": 0,
            "water_defect_largest_cluster_threshold": 3,
            "water_defect_persistence_frames": 5,
            "maximum_absolute_scd_adjacent_block_change": 0.1,
            "maximum_absolute_scd_first_last_change": 0.1,
        },
        "external_metrics": {
            "protein_aware_area_per_lipid": {"status": "not_available", "reason": "Synthetic test confirms that cell lateral area is not relabeled as protein-aware area per lipid."},
            "popc_deuterium_order_parameters": {"status": "not_available", "reason": "Synthetic topology deliberately lacks a validated CHARMM36 carbon-hydrogen chain mapping."},
        },
    }
    membrane_path = root / "membrane_mapping.json"
    atomic_write_json(membrane_path, membrane)
    energy = {
        "schema_version": "1.0",
        "approval_status": "synthetic_self_test",
        "gromacs_time_unit": "ps",
        "selection_mode": "one_exact_named_term_per_gmx_energy_invocation",
        "terms": [
            {"key": key, "gmx_name": name, "unit": unit, "required": True}
            for key, name, unit in zip(
                gmx_energy_qc.REQUIRED_TERM_KEYS,
                ("Temperature", "Pressure", "Pres-XX", "Pres-YY", "Pres-ZZ", "Pres-XY", "Pres-XZ", "Pres-YZ", "Potential", "Kinetic-En.", "Total-Energy", "Density", "Volume", "Box-X", "Box-Y", "Box-Z"),
                ("K", "bar", "bar", "bar", "bar", "bar", "bar", "bar", "kJ/mol", "kJ/mol", "kJ/mol", "kg/m^3", "nm^3", "nm", "nm", "nm"),
                strict=True,
            )
        ],
    }
    energy_path = root / "energy_terms.json"
    atomic_write_json(energy_path, energy)
    return structural_path, membrane_path, energy_path


def write_existing_xvg(root: Path, energy_terms_path: Path) -> None:
    terms = json.loads(energy_terms_path.read_text(encoding="utf-8"))["terms"]
    synthetic_values = {
        "temperature_k": 310.0,
        "pressure_bar": 1.0,
        "pressure_xx_bar": 0.9,
        "pressure_yy_bar": 1.0,
        "pressure_zz_bar": 1.1,
        "pressure_xy_bar": 0.0,
        "pressure_xz_bar": 0.0,
        "pressure_yz_bar": 0.0,
        "potential_energy_kj_mol": -1000.0,
        "kinetic_energy_kj_mol": 300.0,
        "total_energy_kj_mol": -700.0,
        "density_kg_m3": 1000.0,
        "volume_nm3": 216.0,
        "box_x_nm": 6.0,
        "box_y_nm": 6.0,
        "box_z_nm": 6.0,
    }
    for realization_index, realization_id in enumerate(("rep01", "rep02", "rep03"), start=1):
        directory = root / realization_id
        directory.mkdir(parents=True)
        for term_index, term in enumerate(terms, start=1):
            selection = directory / f"selection__{term['key']}.txt"
            selection.write_text(f"{term['gmx_name']}\n0\n", encoding="utf-8", newline="\n")
            command = directory / f"command__{term['key']}.json"
            atomic_write_json(command, {"argv": ["gmx", "energy", "-f", "synthetic.edr", "-o", f"raw__{term['key']}.xvg"], "stdin_file": selection.name})
            log = directory / f"gmx_energy__{term['key']}.log"
            log.write_text("returncode=0\nsynthetic parser input; no real EDR execution\n", encoding="utf-8", newline="\n")
            xvg = directory / f"raw__{term['key']}.xvg"
            base_value = synthetic_values[term["key"]]
            lines = ["# synthetic exact-term XVG", f'@ s0 legend "{term["gmx_name"]}"']
            for frame in range(501):
                lines.append(f"{frame * 1000.0:.1f} {base_value:.8f}")
            xvg.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def write_validated_external_membrane_mapping(root: Path, base_mapping: Path, trajectories: list[Path]) -> tuple[Path, Path]:
    mapping = json.loads(base_mapping.read_text(encoding="utf-8"))

    def record(path: Path) -> dict[str, str]:
        return {"path": str(path), "sha256": sha256_file(path)}

    common_records: dict[str, dict[str, dict[str, str]]] = {}
    for metric, tool_version in (("apl", "3.3"), ("scd", "synthetic-gorder-1.0")):
        files = {}
        for label in ("source_code", "version_capture", "command", "validation_report"):
            path = root / f"external_{metric}_{label}.txt"
            path.write_text(f"Synthetic {metric} {label}; deterministic evidence only; tool version {tool_version}.\n", encoding="utf-8", newline="\n")
            files[f"{label}_record"] = record(path)
        atom_mapping_path = root / f"external_{metric}_atom_mapping.json"
        if metric == "apl":
            atom_mapping_payload = {
                "schema_version": "1.0",
                "system_id": "8kct_nirogacestat_native",
                "metric": "protein_aware_area_per_lipid",
                "trajectory_atom_identity_sha256": mapping["trajectory_atom_identity_sha256"],
                "protein_footprint_included": True,
                "lipid_resnames": ["POPC"],
                "protein_atom_indices_sha256": "a" * 64,
                "popc_atom_indices_sha256": "b" * 64,
            }
        else:
            atom_mapping_payload = {
                "schema_version": "1.0",
                "system_id": "8kct_nirogacestat_native",
                "metric": "popc_deuterium_order_parameters",
                "trajectory_atom_identity_sha256": mapping["trajectory_atom_identity_sha256"],
                "force_field_family": "CHARMM36",
                "lipid_resname": "POPC",
                "unsaturated_chain_geometry_explicit": True,
                "carbon_hydrogen_mappings": [
                    {"chain_id": "sn1", "carbon_id": "C2", "carbon_atom_name": "C22", "hydrogen_atom_names": ["H2R", "H2S"]},
                    {"chain_id": "sn2", "carbon_id": "C2", "carbon_atom_name": "C32", "hydrogen_atom_names": ["H2X", "H2Y"]},
                ],
            }
        atomic_write_json(atom_mapping_path, atom_mapping_payload)
        files["atom_mapping_record"] = record(atom_mapping_path)
        common_records[metric] = files

    apl_outputs = []
    scd_outputs = []
    apl_mapping_sha = common_records["apl"]["atom_mapping_record"]["sha256"]
    scd_mapping_sha = common_records["scd"]["atom_mapping_record"]["sha256"]
    first_apl_path: Path | None = None
    for realization_id, trajectory in zip(("rep01", "rep02", "rep03"), trajectories, strict=True):
        trajectory_sha = sha256_file(trajectory)
        apl_path = root / f"{realization_id}_external_apl.csv"
        atomic_write_csv(
            apl_path,
            ["realization_id", "time_ns", "protein_aware_popc_area_per_lipid_nm2", "source_trajectory_sha256", "atom_mapping_sha256", "tool_version"],
            ({
                "realization_id": realization_id,
                "time_ns": float(frame),
                "protein_aware_popc_area_per_lipid_nm2": 0.65,
                "source_trajectory_sha256": trajectory_sha,
                "atom_mapping_sha256": apl_mapping_sha,
                "tool_version": "3.3",
            } for frame in range(501)),
        )
        first_apl_path = first_apl_path or apl_path
        apl_outputs.append({"realization_id": realization_id, "output": record(apl_path)})
        scd_path = root / f"{realization_id}_external_scd.csv"
        scd_rows = []
        for block in range(5):
            for chain_id, carbon_id in (("sn1", "C2"), ("sn2", "C2")):
                scd_rows.append({
                    "realization_id": realization_id,
                    "block_index_zero_based": block,
                    "block_start_ns": 200.0 + block * 60.0,
                    "block_end_ns": 260.0 + block * 60.0,
                    "chain_id": chain_id,
                    "carbon_id": carbon_id,
                    "s_cd": -0.2,
                    "source_trajectory_sha256": trajectory_sha,
                    "atom_mapping_sha256": scd_mapping_sha,
                    "tool_version": "synthetic-gorder-1.0",
                })
        atomic_write_csv(
            scd_path,
            ["realization_id", "block_index_zero_based", "block_start_ns", "block_end_ns", "chain_id", "carbon_id", "s_cd", "source_trajectory_sha256", "atom_mapping_sha256", "tool_version"],
            scd_rows,
        )
        scd_outputs.append({"realization_id": realization_id, "output": record(scd_path)})
    mapping["external_metrics"] = {
        "protein_aware_area_per_lipid": {
            "status": "validated",
            "frozen_tool_route_before_production": True,
            "tool": {"name": "APL@Voro", "version": "3.3"},
            "output_schema": "per_saved_frame_protein_aware_popc_apl_v1",
            **common_records["apl"],
            "per_realization_outputs": apl_outputs,
        },
        "popc_deuterium_order_parameters": {
            "status": "validated",
            "frozen_tool_route_before_production": True,
            "tool": {"name": "gorder", "version": "synthetic-gorder-1.0"},
            "output_schema": "five_fixed_60ns_blocks_charmm36_popc_scd_v1",
            **common_records["scd"],
            "per_realization_outputs": scd_outputs,
        },
    }
    path = root / "membrane_mapping_external_validated.json"
    atomic_write_json(path, mapping)
    assert first_apl_path is not None
    return path, first_apl_path


def write_manifest(root: Path, topology: Path, trajectories: list[Path], structural: Path, membrane: Path, energy: Path) -> Path:
    dummy_edr = root / "synthetic.edr"
    dummy_log = root / "synthetic.log"
    dummy_edr.write_bytes(b"synthetic EDR placeholder for parse-existing self-test only\n")
    dummy_log.write_text("synthetic production log\n", encoding="utf-8", newline="\n")
    acceptance_source = root / "synthetic_acceptance_gate_source.txt"
    acceptance_source.write_text("Synthetic-only gate provenance for deterministic end-to-end validation; never a scientific cutoff source.\n", encoding="utf-8", newline="\n")

    def file_record(path: Path) -> dict[str, str]:
        return {"path": str(path), "sha256": sha256_file(path)}

    manifest = {
        "schema_version": "1.0",
        "approval_status": "synthetic_self_test",
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "system_id": "8kct_nirogacestat_native",
        "construction_count": 1,
        "required_realization_ids": ["rep01", "rep02", "rep03"],
        "production_duration_ns": 500.0,
        "primary_analysis_window_ns": [200.0, 500.0],
        "extension_or_recovery_window": False,
        "required_mdanalysis_version": "2.10.x",
        "reference": {"topology": file_record(topology), "coordinates": file_record(topology)},
        "mapping_records": {"structural": file_record(structural), "membrane": file_record(membrane), "gromacs_energy_terms": file_record(energy)},
        "time_contract": {
            "mdanalysis_trajectory_time_unit": "ps",
            "reported_time_unit": "ns",
            "required_start_ns": 0.0,
            "required_end_ns": 500.0,
            "endpoint_tolerance_ns": 0.0001,
            "strictly_increasing": True,
            "uniform_saved_step": True,
            "identical_saved_times_across_realizations": True,
        },
        "data_handling": {"frame_stride": 1, "smooth_frames": False, "delete_frames": False, "interpolate_frames": False, "drop_realizations": False, "retain_full_0_500_ns_traces": True},
        "inference_contract": {"frame_level_hypothesis_tests": False, "frames_are_independent_replicates": False, "realizations_are_biological_replicates": False, "mmgbsa_mmpbsa": False, "binding_affinity_inference": False, "between_system_comparison": False},
        "diagnostics": {
            "autocorrelation_method": "geyer_initial_positive_sequence",
            "max_lag_fraction": 0.25,
            "block_tau_multiplier": 10.0,
            "minimum_complete_blocks": 5,
            "robust_first_difference_z_threshold": 12.0,
            "stationarity": {
                "method": "fixed_time_block_median_drift_and_change",
                "window_ns": [200.0, 500.0],
                "fixed_time_blocks": 5,
                "minimum_frames_per_block": 5,
                "scale_estimator": "max_1.4826_mad_and_prespecified_metric_floor",
                "maximum_abs_normalized_linear_change": 3.0,
                "maximum_abs_normalized_first_last_shift": 3.0,
                "maximum_abs_normalized_adjacent_shift": 3.0,
                "maximum_abs_normalized_change_point_shift": 3.0,
            },
        },
        "acceptance_gates": {
            "approval_status": "synthetic_self_test",
            "frozen_before_production_review": True,
            "rationale": "Synthetic-only deterministic criteria exercise scientific fail-closed branches without supplying production cutoffs.",
            "source_record": file_record(acceptance_source),
            "stationarity_scale_floors": {
                "structural": {
                    "pocket_aligned_o6u_heavy_rmsd_nm": 1.0,
                    "pocket_aligned_o6u_com_displacement_nm": 1.0,
                    "tm_core_ca_rmsd_nm": 1.0,
                    "protein_ca_rmsd_nm": 1.0,
                    "native_contact_fraction": 1.0,
                },
                "membrane": {
                    "phosphate_peak_thickness_nm": 10.0,
                    "protein_aware_area_per_lipid_nm2": 10.0,
                    "cell_lateral_area_nm2_not_apl": 100.0,
                    "box_z_vector_length_nm": 10.0,
                    "cell_volume_nm3": 1000.0,
                    "protein_tilt_deg": 90.0,
                },
                "energy": {key: 1e6 for key in gmx_energy_qc.REQUIRED_TERM_KEYS},
            },
            "native_pose": {
                "units": {"distance": "nm", "fraction": "unitless"},
                "window_ns": [200.0, 500.0],
                "event_search_window_ns": [0.0, 500.0],
                "continuous_event_rule": "each_geometry_evaluated_separately_without_gap_bridging",
                "minimum_continuous_event_duration_ns": 1.0,
                "maximum_pocket_aligned_o6u_heavy_rmsd_nm": 0.2,
                "maximum_o6u_com_displacement_nm": 0.2,
                "minimum_native_contact_fraction": 0.5,
                "minimum_fraction_of_primary_frames_meeting_all_pose_gates": 0.9,
                "ligand_egress_or_contact_loss_is_scientific_failure": True,
                "failure_triggers_rerun_or_extension": False,
            },
            "thermodynamic_cell_qc": {
                "units": {"temperature": "K", "pressure": "bar", "density": "kg/m^3", "relative_closure": "unitless"},
                "window_ns": [200.0, 500.0],
                "target_temperature_k": 310.0,
                "maximum_absolute_primary_mean_temperature_deviation_k": 1.0,
                "approved_primary_mean_pressure_range_bar": [-5.0, 5.0],
                "approved_primary_mean_density_range_kg_m3": [990.0, 1010.0],
                "maximum_relative_total_energy_closure_error": 1e-10,
                "maximum_absolute_pressure_trace_closure_bar": 1e-10,
                "maximum_relative_orthorhombic_volume_closure_error": 1e-10,
                "failure_triggers_rerun_or_extension": False,
            },
        },
        "energy_execution": {"server_execution_authorized": False, "gmx_executable": "gmx", "required_gromacs_version": "synthetic"},
        "realizations": [
            {
                "realization_id": realization_id,
                "velocity_seed": seed,
                "topology": file_record(topology),
                "centered_system_trajectory": file_record(trajectory),
                "energy_edr": file_record(dummy_edr),
                "production_log": file_record(dummy_log),
            }
            for realization_id, seed, trajectory in zip(("rep01", "rep02", "rep03"), (11111, 22222, 33333), trajectories, strict=True)
        ],
    }
    path = root / "manifest.json"
    atomic_write_json(path, manifest)
    return path


def run_self_test() -> None:
    require(mda.__version__.startswith("2.10."), f"Self-test requires MDAnalysis 2.10.x; observed {mda.__version__}")
    gmx_energy_qc.synthetic_xvg_self_test()
    event_test = continuous_true_events([False, True, True, False], [0.0, 1.0, 2.0, 3.0], 1.0, "short_egress")
    require(len(event_test) == 1 and event_test[0]["continuous_duration_ns"] == 1.0, "A short prespecified continuous ligand event was missed")
    require(not continuous_true_events([True, False, True], [0.0, 1.0, 2.0], 1.0, "gapped"), "Ligand-event detector improperly bridged a gap")
    drift_contract = {
        "method": "fixed_time_block_median_drift_and_change",
        "fixed_time_blocks": 5,
        "minimum_frames_per_block": 5,
        "maximum_abs_normalized_linear_change": 0.5,
        "maximum_abs_normalized_first_last_shift": 0.5,
        "maximum_abs_normalized_adjacent_shift": 0.5,
        "maximum_abs_normalized_change_point_shift": 0.5,
    }
    drift_times = np.linspace(200.0, 500.0, 301)
    drift_result = stationarity_diagnostics(np.linspace(0.0, 10.0, 301), drift_times, drift_contract, 0.1)
    require(drift_result["status"] == "nonstationary" and drift_result["failed_gates"], "A deterministic drift bypassed stationarity gates")
    step_values = np.zeros(301, dtype=np.float64)
    step_values[180:] = 2.0
    change_result = stationarity_diagnostics(step_values, drift_times, drift_contract, 0.1)
    require(change_result["status"] == "nonstationary" and "maximum_abs_normalized_change_point_shift" in change_result["failed_gates"], "A deterministic change point bypassed stationarity gates")
    far_count, _ = analyze_membrane_qc_mdanalysis._outside_protein_water(
        np.asarray([[1.0, 1.0, 0.0]]), np.asarray([0.0]), np.asarray([[1.0, 1.0, 50.0]]), np.asarray([5.0]),
        6.0, 6.0, 10.0, 1.2, 0.5, 0.35,
    )
    near_count, _ = analyze_membrane_qc_mdanalysis._outside_protein_water(
        np.asarray([[1.0, 1.0, 0.0]]), np.asarray([0.0]), np.asarray([[1.0, 1.0, 1.0]]), np.asarray([0.1]),
        6.0, 6.0, 10.0, 1.2, 0.5, 0.35,
    )
    require(far_count == 1 and near_count == 0, "Ectodomain XY projection still masks membrane-core water or true 3D protein exclusion failed")
    with tempfile.TemporaryDirectory(prefix="primary_postprocessing_selftest_", ignore_cleanup_errors=True) as temporary:
        root = Path(temporary)
        topology, trajectories = write_topology_and_trajectories(root)
        structural_mapping, membrane_mapping, energy_terms = write_mapping_records(root, topology)
        existing_xvg = root / "existing_xvg"
        write_existing_xvg(existing_xvg, energy_terms)
        manifest = write_manifest(root, topology, trajectories, structural_mapping, membrane_mapping, energy_terms)
        output = root / "output"
        structural_report = analyze_primary_structure_mdanalysis.run(manifest, output, allow_synthetic=True)
        membrane_report = analyze_membrane_qc_mdanalysis.run(manifest, output, allow_synthetic=True)
        energy_report = gmx_energy_qc.run(manifest, output, "parse-existing", existing_xvg, allow_synthetic=True)
        require(structural_report["status"] == "pass", f"Synthetic structural analysis did not pass: {structural_report}")
        require(membrane_report["status"] == "inconclusive", f"Synthetic membrane analysis did not preserve the external-metric NO-GO: {membrane_report}")
        require(energy_report["status"] == "pass", f"Synthetic energy analysis did not pass: {energy_report}")
        validation = validate_primary_postprocessing.validate(manifest, output, allow_synthetic=True)
        require(validation["status"] == "inconclusive" and validation["claim_gate"] == "blocked_inconclusive", f"Synthetic full validation bypassed the external-metric NO-GO: {validation}")
        require(any("protein_aware_area_per_lipid" in reason for reason in validation["reasons"]), "Protein-aware APL NO-GO is absent")
        require(any("popc_deuterium_order_parameters" in reason for reason in validation["reasons"]), "POPC order-parameter NO-GO is absent")

        adjudication_policy = root / "prospective_spike_policy.txt"
        adjudication_evidence = root / "matched_spike_evidence.txt"
        adjudication_policy.write_text("Synthetic prospective policy frozen before production.\n", encoding="utf-8", newline="\n")
        adjudication_evidence.write_text("Synthetic matched log/frame evidence.\n", encoding="utf-8", newline="\n")
        synthetic_structural_flag = {"realization_id": "rep03", "metric": "protein_ca_rmsd_nm", "row_index_zero_based": "30", "time_after_ns": "30.0"}
        synthetic_energy_flag = {"realization_id": "rep01", "term_key": "pressure_bar", "row_index_zero_based": "10", "time_after_ns": "10.0"}
        synthetic_membrane_flag = {"realization_id": "rep02", "metric": "protein_tilt_deg", "row_index_zero_based": "20", "time_after_ns": "20.0"}
        disposition_record = {
            "schema_version": "3.0",
            "approval_status": "approved",
            "policy_frozen_before_production": True,
            "policy_source_record": {"path": str(adjudication_policy), "sha256": sha256_file(adjudication_policy)},
            "source_complete_sha256": {
                "structural": sha256_file(output / "structural_analysis" / "COMPLETE.json"),
                "energy": sha256_file(output / "energy_qc" / "COMPLETE.json"),
                "membrane": sha256_file(output / "membrane_qc" / "COMPLETE.json"),
            },
            "rules": {
                "all_flags_covered_exactly_once": True,
                "source_rows_retained": True,
                "unresolved_or_blocking_disposition_blocks_claim": True,
                "every_adjudication_requires_source_hashed_evidence": True,
            },
            "allowed_dispositions": {
                "structural": ["finite_coordinate_fluctuation_retained", "pbc_or_coordinate_artifact_raw_point_retained_and_analysis_blocked", "structural_disruption_or_corruption_analysis_blocked"],
                "energy": ["finite_physical_or_barostat_fluctuation_retained", "pbc_or_output_artifact_raw_point_retained_and_analysis_blocked", "constraint_or_corruption_failure_analysis_blocked"],
                "membrane": ["finite_physical_fluctuation_retained", "pbc_or_coordinate_artifact_raw_point_retained_and_analysis_blocked", "structural_disruption_or_corruption_analysis_blocked"],
            },
            "structural_flags": [{**synthetic_structural_flag, "disposition": "finite_coordinate_fluctuation_retained", "evidence_records": [{"path": str(adjudication_evidence), "sha256": sha256_file(adjudication_evidence)}], "reviewer": "synthetic", "reviewed_at_utc": "2026-08-08T00:00:00Z"}],
            "energy_flags": [{**synthetic_energy_flag, "disposition": "finite_physical_or_barostat_fluctuation_retained", "evidence_records": [{"path": str(adjudication_evidence), "sha256": sha256_file(adjudication_evidence)}], "reviewer": "synthetic", "reviewed_at_utc": "2026-08-08T00:00:00Z"}],
            "membrane_flags": [{**synthetic_membrane_flag, "disposition": "finite_physical_fluctuation_retained", "evidence_records": [{"path": str(adjudication_evidence), "sha256": sha256_file(adjudication_evidence)}], "reviewer": "synthetic", "reviewed_at_utc": "2026-08-08T00:00:00Z"}],
        }
        disposition_path = root / "approved_spike_dispositions.json"
        atomic_write_json(disposition_path, disposition_record)
        disposition_reasons = validate_primary_postprocessing._load_review_dispositions(
            disposition_path,
            output / "structural_analysis" / "COMPLETE.json",
            output / "energy_qc" / "COMPLETE.json",
            output / "membrane_qc" / "COMPLETE.json",
            [synthetic_structural_flag],
            [synthetic_energy_flag],
            [synthetic_membrane_flag],
        )
        require(not disposition_reasons, "Valid source-hashed spike adjudication was rejected")
        with adjudication_evidence.open("a", encoding="utf-8", newline="") as handle:
            handle.write("tamper\n")
        adjudication_tamper_rejected = False
        try:
            validate_primary_postprocessing._load_review_dispositions(
                disposition_path,
                output / "structural_analysis" / "COMPLETE.json",
                output / "energy_qc" / "COMPLETE.json",
                output / "membrane_qc" / "COMPLETE.json",
                [synthetic_structural_flag],
                [synthetic_energy_flag],
                [synthetic_membrane_flag],
            )
        except ContractError:
            adjudication_tamper_rejected = True
        require(adjudication_tamper_rejected, "Tampered spike evidence bypassed disposition validation")

        validated_membrane_mapping, first_apl_output = write_validated_external_membrane_mapping(root, membrane_mapping, trajectories)
        validated_manifest_record = json.loads(manifest.read_text(encoding="utf-8"))
        validated_manifest_record["mapping_records"]["membrane"] = {"path": str(validated_membrane_mapping), "sha256": sha256_file(validated_membrane_mapping)}
        validated_manifest = root / "manifest_external_membrane_validated.json"
        atomic_write_json(validated_manifest, validated_manifest_record)
        validated_external_output = root / "validated_external_output"
        validated_membrane_report = analyze_membrane_qc_mdanalysis.run(validated_manifest, validated_external_output, allow_synthetic=True)
        require(validated_membrane_report["status"] == "pass" and validated_membrane_report["preproduction_status"] == "pass", "Source-hashed APL/gorder route was not reachable")
        with first_apl_output.open("a", encoding="utf-8", newline="") as handle:
            handle.write("tamper\n")
        external_tamper_rejected = False
        try:
            analyze_membrane_qc_mdanalysis.run(validated_manifest, root / "tampered_external_output", allow_synthetic=True)
        except ContractError:
            external_tamper_rejected = True
        require(external_tamper_rejected, "A tampered external APL output bypassed source-hash validation")

        failure_manifest_record = json.loads(manifest.read_text(encoding="utf-8"))
        failure_manifest_record["acceptance_gates"]["native_pose"]["maximum_pocket_aligned_o6u_heavy_rmsd_nm"] = 1e-12
        failure_manifest_record["acceptance_gates"]["thermodynamic_cell_qc"]["target_temperature_k"] = 100.0
        failure_manifest = root / "manifest_scientific_failure.json"
        atomic_write_json(failure_manifest, failure_manifest_record)
        failure_output = root / "scientific_failure_output"
        structural_failure = analyze_primary_structure_mdanalysis.run(failure_manifest, failure_output, allow_synthetic=True)
        energy_failure = gmx_energy_qc.run(failure_manifest, failure_output, "parse-existing", existing_xvg, allow_synthetic=True)
        require(structural_failure["status"] == "inconclusive" and all(item["scientific_status"] == "fail" for item in structural_failure["realization_summaries"]), "Ligand pose failure did not fail all realizations")
        require(energy_failure["status"] == "inconclusive" and all(item["scientific_status"] == "fail" for item in energy_failure["realization_summaries"]), "Thermodynamic failure did not fail all realizations")

        summary_path = output / "structural_analysis" / "rep01" / "structural_summary.json"
        original_summary = summary_path.read_bytes()
        tampered_summary = json.loads(original_summary.decode("utf-8"))
        tampered_summary["diagnostics"]["protein_ca_rmsd_nm"]["primary_200_500_ns"]["stationarity"]["observed"]["abs_normalized_linear_change"] += 1.0
        atomic_write_json(summary_path, tampered_summary)
        stationarity_tamper_rejected = False
        try:
            validate_primary_postprocessing.validate(manifest, output, allow_synthetic=True)
        except ContractError:
            stationarity_tamper_rejected = True
        require(stationarity_tamper_rejected, "Validator accepted a tampered stationarity report")
        summary_path.write_bytes(original_summary)

        raw_path = output / "structural_analysis" / "rep02" / "structural_raw_unsmoothed.csv"
        with raw_path.open("a", encoding="utf-8", newline="") as handle:
            handle.write("tamper\n")
        rejected = False
        try:
            validate_primary_postprocessing.validate(manifest, output, allow_synthetic=True)
        except ContractError:
            rejected = True
        require(rejected, "Validator did not reject a hash/row tamper")
    print("END-TO-END SELF-TEST PASS: fixed-window drift/change rejection, short continuous ligand-event failure, 3D membrane-core water geometry, source-hashed spike adjudication, reachable APL/gorder ingestion, external-output tamper rejection, all raw rows retained, and no rescue extension.")


if __name__ == "__main__":
    run_self_test()
