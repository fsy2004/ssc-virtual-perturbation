#!/usr/bin/env python3
"""Adversarial regression tests for the immutable MD release gates."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from charmmgui_api import DEFAULT_BASE, LOCAL_TEST_TOKEN, validate_base_url, validate_request_destination
from md_contract import build_contract_sha256, load_json, report_payload_sha256, validate_production_mdp
from validate_preflight import Audit, validate_build_report_binding


ROOT = Path(__file__).resolve().parent.parent


def production_manifest() -> dict:
    return {
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "global_model": {"temperature_k": 303.15, "pressure_bar": 1.0},
        "simulation": {
            "production_ns": 500.0, "time_step_ps": 0.004, "constraints": "h-bonds",
            "hydrogen_mass_repartitioning": True,
            "pressure_coupling": "semiisotropic",
            "production_mdp_contract": {
                "integrator": "md", "thermostat": "v-rescale", "thermostat_groups": ["SOLU", "MEMB", "SOLV"],
                "tau_t_ps": 1.0, "barostat": "C-rescale", "barostat_tau_p_ps": 5.0,
                "compressibility_bar_inverse": [4.5e-5, 4.5e-5], "cutoff_scheme": "verlet",
                "neighbor_list_update_steps": 20, "rlist_nm": 1.2, "rcoulomb_nm": 1.2,
                "vdw_type": "cut-off", "vdw_modifier": "force-switch", "rvdw_switch_nm": 1.0,
                "rvdw_nm": 1.2, "dispersion_correction": "no", "pme_order": 4,
                "fourier_spacing_nm": 0.12, "constraint_algorithm": "lincs",
                "com_removal_mode": "linear", "com_removal_groups": ["SOLU_MEMB", "SOLV"],
                "com_removal_interval_steps": 100,
                "output_cadence_steps": {"nstxout": 0, "nstvout": 0, "nstfout": 0,
                    "nstxout_compressed": 5000, "nstcalcenergy": 100, "nstenergy": 5000,
                    "nstlog": 5000, "compressed_x_precision": 1000},
            },
        },
        "systems": [{"id": "8kct_nirogacestat_native", "pdb_id": "8KCT", "ligand_component_id": "O6U",
                     "construction": {"id": "build01", "pdb_reader_jobid": "pdb1", "quick_bilayer_jobid": "qb1",
                                      "charmm_gui_archive": {"path": "archive.tgz", "sha256": "a" * 64},
                                      "gromacs_input_tree_manifest": {"path": "tree.json", "sha256": "b" * 64},
                                      "minimization_mdp": {"path": "step6.0_minimization.mdp", "sha256": "c" * 64},
                                      "equilibration_mdp_sha256": {
                                          f"step6.{index}_equilibration.mdp": str(index) * 64 for index in range(1, 7)
                                      }}}],
    }


def valid_mdp() -> str:
    return """integrator=md
dt=0.004
nsteps=125000000
continuation=yes
gen_vel=no
pbc=xyz
periodic-molecules=no
cutoff-scheme=Verlet
nstlist=20
rlist=1.2
coulombtype=PME
rcoulomb=1.2
pme-order=4
fourierspacing=0.12
vdwtype=Cut-off
vdw-modifier=Force-switch
rvdw-switch=1.0
rvdw=1.2
DispCorr=no
constraints=h-bonds
constraint-algorithm=lincs
tcoupl=v-rescale
tc-grps=SOLU MEMB SOLV
tau-t=1.0 1.0 1.0
ref-t=303.15 303.15 303.15
pcoupl=C-rescale
pcoupltype=semiisotropic
tau-p=5.0
ref-p=1.0 1.0
compressibility=4.5e-5 4.5e-5
comm-mode=linear
comm-grps=SOLU_MEMB SOLV
nstcomm=100
nstxout=0
nstvout=0
nstfout=0
nstxout-compressed=5000
nstcalcenergy=100
nstenergy=5000
nstlog=5000
compressed-x-precision=1000
"""


def main() -> int:
    if validate_base_url(DEFAULT_BASE, False) != DEFAULT_BASE:
        raise RuntimeError("official CHARMM-GUI endpoint was not accepted")
    for hostile in ("https://example.invalid/api", "https://charmm-gui.org.evil.invalid/api", "https://charmm-gui.org/api"):
        try:
            validate_base_url(hostile, False)
        except ValueError:
            pass
        else:
            raise RuntimeError(f"untrusted bearer-token destination passed: {hostile}")
    if validate_base_url("http://127.0.0.1:8765/api", True) != "http://127.0.0.1:8765/api":
        raise RuntimeError("explicit localhost test endpoint failed")
    validate_request_destination(f"{DEFAULT_BASE}/check_status?jobid=test", "synthetic-real-token")
    validate_request_destination("http://127.0.0.1:8765/api/download?jobid=test", LOCAL_TEST_TOKEN)
    for url, token in (
        ("https://example.invalid/api/check_status?jobid=test", "synthetic-real-token"),
        ("https://www.charmm-gui.org.evil.invalid/api/download?jobid=test", "synthetic-real-token"),
        ("http://127.0.0.1:8765/api/download?jobid=test", "synthetic-real-token"),
        (f"{DEFAULT_BASE}/unapproved", "synthetic-real-token"),
    ):
        try:
            validate_request_destination(url, token)
        except ValueError:
            pass
        else:
            raise RuntimeError(f"Authorization header could reach an unpinned destination: {url}")

    manifest = production_manifest()
    changed = copy.deepcopy(manifest)
    changed["systems"][0]["construction"]["quick_bilayer_jobid"] = "wrong-job"
    if build_contract_sha256(manifest) == build_contract_sha256(changed):
        raise RuntimeError("mismatched build/job report contract was not detected")
    changed = copy.deepcopy(manifest)
    changed["systems"][0]["construction"]["charmm_gui_archive"]["sha256"] = "b" * 64
    if build_contract_sha256(manifest) == build_contract_sha256(changed):
        raise RuntimeError("mismatched build archive was not detected")
    construction = manifest["systems"][0]["construction"]
    build_report = {
        "schema_version": "2.0", "report_type": "charmm_gui_build_validation", "status": "pass", "strict": True,
        "study_id": manifest["study_id"], "system_id": "8kct_nirogacestat_native", "construction_id": "build01",
        "pdb_reader_jobid": "pdb1", "quick_bilayer_jobid": "qb1", "archive_sha256": "a" * 64,
        "build_contract_sha256": build_contract_sha256(manifest), "hydrogen_mass_repartitioning_detected": True,
        "equilibration_mdp_sha256": construction["equilibration_mdp_sha256"],
        "gromacs_input_tree_manifest_sha256": construction["gromacs_input_tree_manifest"]["sha256"],
        "membrane_orientation_record_sha256": None, "starting_coordinates_sha256": None,
        "topology_sha256": None, "index_sha256": None, "analysis_index_sha256": None,
        "production_mdp_sha256": None, "minimization_mdp_sha256": construction["minimization_mdp"]["sha256"],
        "archive_extracted_gromacs_tree_binding": [{
            "relative_path": "topol.top", "match": True, "archive_bytes": 10, "extracted_bytes": 10,
            "archive_sha256": "d" * 64, "extracted_sha256": "d" * 64
        }],
        "staged_mdp_archive_binding_and_physics": ([{
            "name": "step6.0_minimization.mdp", "archive_sha256": construction["minimization_mdp"]["sha256"],
            "manifest_sha256": construction["minimization_mdp"]["sha256"], "physics": {"integrator": "steep"}
        }] + [{
            "name": f"step6.{index}_equilibration.mdp",
            "archive_sha256": construction["equilibration_mdp_sha256"][f"step6.{index}_equilibration.mdp"],
            "manifest_sha256": construction["equilibration_mdp_sha256"][f"step6.{index}_equilibration.mdp"],
            "physics": {"integrator": "md"}
        } for index in range(1, 7)]),
        "integrity": {"payload_sha256": "UNSEALED"},
    }
    build_report["integrity"]["payload_sha256"] = report_payload_sha256(build_report, ("integrity", "payload_sha256"))
    valid_audit = Audit(); validate_build_report_binding(build_report, construction, manifest, valid_audit)
    if valid_audit.errors:
        raise RuntimeError(f"valid synthetic build binding failed: {valid_audit.errors}")
    for key, bad_value in (("archive_sha256", "b" * 64), ("quick_bilayer_jobid", "wrong-job")):
        mismatched = copy.deepcopy(build_report); mismatched[key] = bad_value
        mismatched["integrity"]["payload_sha256"] = report_payload_sha256(mismatched, ("integrity", "payload_sha256"))
        mismatch_audit = Audit(); validate_build_report_binding(mismatched, construction, manifest, mismatch_audit)
        if not mismatch_audit.errors:
            raise RuntimeError(f"mismatched build report passed: {key}")
    tampered = copy.deepcopy(build_report); tampered["strict"] = False
    tamper_audit = Audit(); validate_build_report_binding(tampered, construction, manifest, tamper_audit)
    if not tamper_audit.errors:
        raise RuntimeError("manually tampered build report passed")

    with tempfile.TemporaryDirectory(prefix="core_gate_adversarial_") as temporary:
        temp = Path(temporary)
        mdp = temp / "production.mdp"; mdp.write_text(valid_mdp(), encoding="utf-8")
        validate_production_mdp(mdp, manifest)
        for old, new in (("ref-t=303.15 303.15 303.15", "ref-t=100 100 100"),
                         ("pcoupl=C-rescale", "pcoupl=no"),
                         ("DispCorr=no", "DispCorr=EnerPres"),
                         ("nstxout-compressed=5000", "nstxout-compressed=50000"),
                         ("periodic-molecules=no", "periodic-molecules=no\nmass-repartition-factor=1.0")):
            wrong = temp / f"wrong_{old.split('=')[0]}.mdp"
            wrong.write_text(valid_mdp().replace(old, new), encoding="utf-8")
            try:
                validate_production_mdp(wrong, manifest)
            except ValueError:
                pass
            else:
                raise RuntimeError(f"wrong frozen production setting passed: {new}")
        biased = temp / "wrong_pull.mdp"
        biased.write_text(valid_mdp() + "pull=yes\n", encoding="utf-8")
        try:
            validate_production_mdp(biased, manifest)
        except ValueError:
            pass
        else:
            raise RuntimeError("biased production MDP passed")

        staged = load_json(ROOT / "config" / "study_manifest.template.json")
        staged["manifest_status"] = "design_frozen"
        staged["global_model"]["temperature_k"] = 303.15
        staged["simulation"]["time_step_ps"] = 0.004
        staged["simulation"]["hydrogen_mass_repartitioning"] = True
        (temp / "templates").mkdir()
        shutil.copy2(ROOT / "templates" / "realization_record.template.json", temp / "templates" / "realization_record.template.json")
        shutil.copy2(ROOT / "PROTONATION_MODEL.md", temp / "PROTONATION_MODEL.md")
        from md_contract import sha256
        staged["global_model"]["dyad_rationale_record"]["sha256"] = sha256(temp / "PROTONATION_MODEL.md")
        manifest_path = temp / "config" / "manifest.json"; manifest_path.parent.mkdir()
        manifest_path.write_text(json.dumps(staged, indent=2) + "\n", encoding="utf-8")
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "validate_preflight.py"), "--manifest", str(manifest_path),
             "--stage", "design", "--strict"], capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(f"stage-specific strict design preflight incorrectly required future artifacts:\n{result.stderr}")
        chemistry = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "validate_preflight.py"), "--manifest", str(manifest_path),
             "--stage", "chemistry", "--strict"], capture_output=True, text=True, check=False,
        )
        if chemistry.returncode == 0:
            raise RuntimeError("chemistry stage passed without chemistry artifacts/status")
        staged["manifest_status"] = "build_approved"
        manifest_path.write_text(json.dumps(staged, indent=2) + "\n", encoding="utf-8")
        builds = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "validate_preflight.py"), "--manifest", str(manifest_path),
             "--stage", "builds", "--strict"], capture_output=True, text=True, check=False,
        )
        if builds.returncode == 0 or "server environment validation report" not in (builds.stdout + builds.stderr):
            raise RuntimeError("build/minimization release did not require the upstream server-environment report")
        plan = load_json(ROOT / "config" / "analysis_plan.template.json")
        if any(key in plan for key in ("pca", "grid", "support")):
            raise RuntimeError("analysis plan still contains an executable PCA/FEL configuration")
        prohibited = plan.get("prohibited_analyses", {})
        if prohibited.get("policy") != "hard_prohibited_no_supplementary_exception" or not all(
            prohibited.get(key) is True for key in (
                "pca", "occupancy_derived_minus_kbt_ln_p", "free_energy_landscape",
                "population_or_free_energy_surface_3d",
            )
        ):
            raise RuntimeError("analysis plan does not hard-prohibit PCA/FEL/3D surfaces")
    print("SELF-TEST PASS: hostile token hosts, stale build/job/archive, wrong physics/cadence, premature stage release, and PCA/FEL configuration all fail closed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
