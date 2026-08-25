#!/usr/bin/env python3
"""Fail-closed validation for the frozen O6U endpoint-energy contract."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any


FROZEN_PLAN_SHA256 = "63afdee0ae871200297576b02b3c806194039992dbe0e46e06dc191f1dc9a958"
ARCHIVE_SHA256 = "5a421f28afee664b5a8919db5f415f1205f35200950117bb3a67fceaba544a98"
PRODUCTION_RELEASE_SHA256 = "a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa"
REPLICAS = ["rep01", "rep02", "rep03"]
MODELS = ["PB_membrane_indi4"]
WITHDRAWN_MODELS = ["PB_membrane_indi1", "GB_OBC2", "GB_Neck2"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, field: str, detail: str = "frozen value mismatch") -> None:
    if not condition:
        raise ValueError(f"{field}: {detail}")


def exact_keys(payload: dict[str, Any], expected: set[str], field: str) -> None:
    actual = set(payload)
    require(actual == expected, field, f"keys differ; missing={sorted(expected-actual)}, extra={sorted(actual-expected)}")


def validate_plan(path: Path, *, bind_frozen_hash: bool = True) -> dict[str, Any]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    plan_hash = sha256(path)
    if bind_frozen_hash:
        require(plan_hash == FROZEN_PLAN_SHA256, "plan_sha256")

    exact_keys(
        payload,
        {
            "schema_version", "analysis_id", "status", "role", "amendment", "supersedes_amendment",
            "release_archive_sha256", "production_tpr_release_sha256",
            "realizations", "interim_rep01_allowed", "formal_all_three_required",
            "sampling", "groups", "trajectory", "software", "pb_primary",
            "models_to_run", "withdrawn_sensitivity_models_not_run", "sensitivity_withdrawal_reason",
            "entropy", "decomposition",
            "inference", "execution_gate",
        },
        "plan",
    )
    require(payload["schema_version"] == "1.0", "schema_version")
    require(payload["analysis_id"] == "o6u_secondary_endpoint_energy_v2_20260822", "analysis_id")
    require(payload["status"] == "frozen_before_endpoint_energy_results_sensitivity_withdrawn", "status")
    require(payload["role"] == "secondary_exploratory_endpoint_energy", "role")
    require(payload["amendment"] == "MMGBSA_PBSA_SECONDARY_ANALYSIS_SENSITIVITY_WITHDRAWAL_20260822.md", "amendment")
    require(payload["supersedes_amendment"] == "MMGBSA_PBSA_SECONDARY_ANALYSIS_AMENDMENT_20260818.md", "supersedes_amendment")
    require(payload["release_archive_sha256"] == ARCHIVE_SHA256, "release_archive_sha256")
    require(payload["production_tpr_release_sha256"] == PRODUCTION_RELEASE_SHA256, "production_tpr_release_sha256")
    require(payload["realizations"] == REPLICAS, "realizations")
    require(payload["interim_rep01_allowed"] is False, "interim_rep01_allowed")
    require(payload["formal_all_three_required"] is True, "formal_all_three_required")

    sampling = payload["sampling"]
    require(sampling["window_ns"] == [200.0, 500.0], "sampling.window_ns")
    require(sampling["stratum_width_ns"] == 1.0, "sampling.stratum_width_ns")
    require(sampling["selection"] == "nearest_to_stratum_midpoint_ties_earlier", "sampling.selection")
    require(sampling["first_midpoint_ns"] == 200.5, "sampling.first_midpoint_ns")
    require(sampling["last_midpoint_ns"] == 499.5, "sampling.last_midpoint_ns")
    require(sampling["frames_per_realization"] == 300, "sampling.frames_per_realization")
    require(
        sampling["fixed_blocks_ns"]
        == [[200.0, 260.0], [260.0, 320.0], [320.0, 380.0], [380.0, 440.0], [440.0, 500.0]],
        "sampling.fixed_blocks_ns",
    )

    require(payload["groups"] == {
        "ligand": "neutral_76_atom_O6U",
        "receptor": "gamma_secretase_protein_plus_covalent_glycans_plus_3_CLR_plus_2_PC1_DSPC",
        "strip": "bulk_POPC_water_Na_Cl",
    }, "groups")
    trajectory = payload["trajectory"]
    require(trajectory["single_trajectory_protocol"] is True, "trajectory.single_trajectory_protocol")
    require(trajectory["membrane_normal_axis"] == "z", "trajectory.membrane_normal_axis")
    require(trajectory["membrane_midplane_z_angstrom"] == 0.0, "trajectory.membrane_midplane_z_angstrom")
    require(trajectory["pbc_distance_invariance_tolerance_nm"] == 0.01, "trajectory.pbc_distance_invariance_tolerance_nm")
    require(trajectory["raw_trajectory_immutable"] is True, "trajectory.raw_trajectory_immutable")

    software = payload["software"]
    require(software["target_gmx_mmpbsa_version"] == "1.6.5", "software.target_gmx_mmpbsa_version")
    require(software["isolated_environment_required"] is True, "software.isolated_environment_required")
    require(software["three_frame_canary_required"] is True, "software.three_frame_canary_required")
    require(software["capture_versions_environment_commands_logs_and_hashes"] is True, "software.capture_versions_environment_commands_logs_and_hashes")

    pb = payload["pb_primary"]
    require(pb == {
        "memopt": 1, "emem": 7.0, "indi": 4.0, "exdi": 80.0,
        "istrng_molar": 0.15, "poretype": 1, "mctrdz_angstrom": 0.0,
        "mthick_rule": "across_realization_median_of_fixed_window_median_leaflet_phosphate_plane_separation_angstrom",
    }, "pb_primary")
    require(payload["models_to_run"] == MODELS, "models_to_run")
    require(payload["withdrawn_sensitivity_models_not_run"] == WITHDRAWN_MODELS, "withdrawn_sensitivity_models_not_run")
    require(
        payload["sensitivity_withdrawal_reason"]
        == "withdrawn_before_formal_endpoint_energy_results_because_sensitivity_branches_add_low_interpretive_value_and_nonessential_compute_multiplicity",
        "sensitivity_withdrawal_reason",
    )

    entropy = payload["entropy"]
    require(entropy and all(value is False for value in entropy.values()), "entropy")
    decomp = payload["decomposition"]
    require(decomp["enabled"] is True, "decomposition.enabled")
    require(decomp["selection"] == "hash_pinned_native_PLIP_contact_set_only", "decomposition.selection")
    require(decomp["descriptive_only"] is True, "decomposition.descriptive_only")
    require(decomp["data_driven_residue_ranking"] is False, "decomposition.data_driven_residue_ranking")

    inference = payload["inference"]
    require(inference["frames_are_independent_units"] is False, "inference.frames_are_independent_units")
    require(inference["frame_level_p_values"] is False, "inference.frame_level_p_values")
    require(inference["hierarchical_bootstrap_seed"] == 20260818, "inference.hierarchical_bootstrap_seed")
    require(inference["absolute_binding_free_energy_claim"] is False, "inference.absolute_binding_free_energy_claim")
    require(inference["affinity_potency_efficacy_claim"] is False, "inference.affinity_potency_efficacy_claim")
    require(inference["between_ligand_ranking"] is False, "inference.between_ligand_ranking")

    gate = payload["execution_gate"]
    require(gate["all_three_500ns_complete"] is True, "execution_gate.all_three_500ns_complete")
    require(gate["all_required_integrity_pbc_membrane_and_energy_gates_passed"] is True, "execution_gate.all_required_integrity_pbc_membrane_and_energy_gates_passed")
    require(gate["concurrent_with_production"] is False, "execution_gate.concurrent_with_production")

    return {
        "schema_version": "1.0",
        "report_type": "secondary_endpoint_energy_plan_validation",
        "status": "pass",
        "analysis_id": payload["analysis_id"],
        "plan_sha256": plan_hash,
        "release_archive_sha256": ARCHIVE_SHA256,
        "production_tpr_release_sha256": PRODUCTION_RELEASE_SHA256,
        "replicas": REPLICAS,
        "models": MODELS,
        "withdrawn_sensitivity_models_not_run": WITHDRAWN_MODELS,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate_plan(args.plan)
    report["validated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        if args.report.exists():
            raise FileExistsError(f"refusing to overwrite {args.report}")
        args.report.write_text(rendered, encoding="utf-8")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
