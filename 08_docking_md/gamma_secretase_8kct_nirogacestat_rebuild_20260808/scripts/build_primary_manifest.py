#!/usr/bin/env python3
"""Build config/primary_postprocessing_manifest.json from the frozen template.

Static sections (mapping records, diagnostics, acceptance gates, energy
execution) are filled now, BEFORE any production result is inspected. The
realization hashes are placeholders until rep01-rep03 artifacts exist; the
manifest remains approval_status=draft until all hashes are bound and the
record is formally approved.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "templates" / "primary_postprocessing_manifest.template.json"
STRUCTURAL = ROOT / "config" / "primary_atom_mapping_contacts.json"
MEMBRANE = ROOT / "config" / "membrane_qc_mapping.json"
ENERGY = ROOT / "config" / "gromacs_energy_terms.json"
REFERENCE = ROOT / "docking_native_redock" / "plip_native" / "8KCT_protonated.pdb"
FREEZE = ROOT / "ANALYSIS_CONFIG_FREEZE_20260821.md"
OUT = ROOT / "config" / "primary_postprocessing_manifest.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    record = json.loads(TEMPLATE.read_text(encoding="utf-8"))
    record["approval_status"] = "draft_not_for_execution"
    reference_path = str(REFERENCE.relative_to(ROOT)).replace("\\", "/")
    record["reference"] = {
        "topology": {"path": reference_path, "sha256": sha256(REFERENCE)},
        "coordinates": {"path": reference_path, "sha256": sha256(REFERENCE)},
    }
    record["mapping_records"] = {
        "structural": {"path": "config/primary_atom_mapping_contacts.json", "sha256": sha256(STRUCTURAL)},
        "membrane": {"path": "config/membrane_qc_mapping.json", "sha256": sha256(MEMBRANE)},
        "gromacs_energy_terms": {"path": "config/gromacs_energy_terms.json", "sha256": sha256(ENERGY)},
    }
    record["diagnostics"]["stationarity"].update({
        "maximum_abs_normalized_linear_change": 3.0,
        "maximum_abs_normalized_first_last_shift": 3.0,
        "maximum_abs_normalized_adjacent_shift": 3.0,
        "maximum_abs_normalized_change_point_shift": 3.0,
    })
    record["acceptance_gates"].update({
        "approval_status": "draft_not_for_execution",
        "frozen_before_production_review": True,
        "rationale": "ANALYSIS_CONFIG_FREEZE_20260821.md: values fixed from the fail-closed gate "
                     "specification, the accepted 8KCT structure and Guo 2025 Supplementary Table 3, "
                     "the constructed membrane geometry, the frozen production protocol, and "
                     "standard practice; no trajectory-derived number used.",
        "source_record": {"path": "ANALYSIS_CONFIG_FREEZE_20260821.md",
                          "sha256": sha256(FREEZE)},
    })
    floors = record["acceptance_gates"]["stationarity_scale_floors"]
    floors["structural"].update({
        "pocket_aligned_o6u_heavy_rmsd_nm": 0.05,
        "pocket_aligned_o6u_com_displacement_nm": 0.05,
        "tm_core_ca_rmsd_nm": 0.05,
        "protein_ca_rmsd_nm": 0.05,
        "native_contact_fraction": 0.05,
    })
    floors["membrane"].update({
        "phosphate_peak_thickness_nm": 0.05,
        "protein_aware_area_per_lipid_nm2": 0.05,
        "cell_lateral_area_nm2_not_apl": 0.10,
        "box_z_vector_length_nm": 0.05,
        "cell_volume_nm3": 0.50,
        "protein_tilt_deg": 1.0,
    })
    floors["energy"].update({
        "temperature_k": 1.0,
        "pressure_bar": 5.0,
        "pressure_xx_bar": 5.0, "pressure_yy_bar": 5.0, "pressure_zz_bar": 5.0,
        "pressure_xy_bar": 5.0, "pressure_xz_bar": 5.0, "pressure_yz_bar": 5.0,
        "potential_energy_kj_mol": 1000.0,
        "kinetic_energy_kj_mol": 1000.0,
        "total_energy_kj_mol": 1000.0,
        "density_kg_m3": 1.0,
        "volume_nm3": 5.0,
        "box_x_nm": 0.05, "box_y_nm": 0.05, "box_z_nm": 0.05,
    })
    record["acceptance_gates"]["native_pose"].update({
        "minimum_continuous_event_duration_ns": 5.0,
        "maximum_pocket_aligned_o6u_heavy_rmsd_nm": 0.50,
        "maximum_o6u_com_displacement_nm": 0.50,
        "minimum_native_contact_fraction": 0.50,
        "minimum_fraction_of_primary_frames_meeting_all_pose_gates": 0.80,
        "ligand_egress_or_contact_loss_is_scientific_failure": True,
        "failure_triggers_rerun_or_extension": False,
    })
    record["acceptance_gates"]["thermodynamic_cell_qc"].update({
        "target_temperature_k": 303.15,
        "maximum_absolute_primary_mean_temperature_deviation_k": 5.0,
        "approved_primary_mean_pressure_range_bar": [-100.0, 100.0],
        "approved_primary_mean_density_range_kg_m3": [950.0, 1100.0],
        "maximum_relative_total_energy_closure_error": 1.0e-3,
        "maximum_absolute_pressure_trace_closure_bar": 10.0,
        "maximum_relative_orthorhombic_volume_closure_error": 1.0e-3,
        "failure_triggers_rerun_or_extension": False,
    })
    record["energy_execution"].update({
        "server_execution_authorized": True,
        "gmx_executable": "/root/GROMACS-2025.2/bin/gmx",
        "required_gromacs_version": "2025.2",
    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    digest = hashlib.sha256(OUT.read_bytes()).hexdigest()
    OUT.with_suffix(".json.sha256").write_text(f"{digest}  {OUT.name}\n", encoding="ascii")
    print(f"WROTE {OUT} sha256={digest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
