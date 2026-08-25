#!/usr/bin/env python3
"""Upload the analysis chain (scripts/config/builds) to the server release root.

Creates the analysis working layout inside the release directory so that
package_root == release root and realizations run_directory repNN resolve to
the existing production files (no 31 GB XTC copy). Never touches production
files. Prints a SHA-256 verification of every uploaded file.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.new_md_server import connect  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent  # MD project root

REMOTE_BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"

UPLOAD_SCRIPTS = [
    "make_analysis_trajectories.py",
    "prepare_primary_pbc_trajectories.py",
    "resume_primary_pbc_trajectories.py",
    "validate_md_outputs.py",
    "validate_qc_stationarity_report.py",
    "md_contract.py",
    "primary_postprocessing_common.py",
    "analyze_primary_structure_mdanalysis.py",
    "analyze_membrane_qc_mdanalysis.py",
    "gmx_energy_qc.py",
    "validate_primary_postprocessing.py",
    "run_primary_structure_memory_safe.py",
    "_remote_launch_primary_qc_memory_safe_recovery.sh",
    "seal_secondary_endpoint_all_three_gate.py",
    "prepare_secondary_endpoint_energy_inputs.py",
    "freeze_endpoint_energy_membrane_geometry.py",
    "build_endpoint_energy_cpu_migration_package.py",
    "plan_secondary_endpoint_resources.py",
    "collect_endpoint_cpu_inventory.py",
    "run_gmx_mmpbsa_canary.py",
    "run_secondary_endpoint_energy_cpu.py",
    "summarize_secondary_endpoint_energy.py",
    "validate_secondary_endpoint_energy_plan.py",
    "validate_secondary_endpoint_execution_defaults.py",
    "capture_gmx_mmpbsa_toolchain.py",
    "install_endpoint_preprocess_env.sh",
    "install_gmx_mmpbsa_1_6_5_cpu.sh",
    "install_gorder_1_5_0.sh",
    "install_fatslim_0_2_2.sh",
    "_remote_launch_gorder_install.sh",
    "_remote_launch_fatslim_install.sh",
    "_remote_resume_cpu_clone_work.sh",
    "_remote_pause_postproduction_for_cpu_switch.sh",
    "_remote_postproduction_status.sh",
    "build_primary_mapping_records.py",
    "build_membrane_mapping.py",
    "build_energy_terms_record.py",
    "build_primary_manifest.py",
    "seal_primary_postprocessing_manifest.py",
    "build_analysis_ndx.py",
    "build_study_manifest.py",
    "validate_preflight.py",
    "_remote_launch_rep02_pbc.sh",
    "_remote_launch_rep03_pbc.sh",
    "_sync_secondary_withdrawal_amendment_to_server.py",
]

UPLOAD_ROOT_FILES = [
    "MMGBSA_PBSA_SECONDARY_ANALYSIS_SENSITIVITY_WITHDRAWAL_20260822.md",
    "SECONDARY_ENDPOINT_ENERGY_EXECUTION_SUPPLEMENT_20260820.md",
]

UPLOAD_CONFIG = [
    "study_manifest.json",
    "analysis_plan.json",
    "primary_atom_mapping_contacts.json",
    "membrane_qc_mapping.json",
    "gromacs_energy_terms.json",
    "primary_postprocessing_manifest.json",
    "secondary_endpoint_energy_plan_v1.json",
    "secondary_endpoint_energy_execution_defaults_v1.json",
    "production_protocol_hmr4fs_303K_v1.json",
]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_remote_file(sftp, path: str) -> str:
    digest = hashlib.sha256()
    with sftp.open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    c = connect()
    sftp = c.open_sftp()
    try:
        for rel in ("scripts", "config", "builds"):
            sftp.mkdir(f"{REMOTE_BASE}/{rel}")
            print(f"mkdir {REMOTE_BASE}/{rel}")
    except OSError:
        pass
    results = []
    for name in UPLOAD_SCRIPTS:
        local = ROOT / "scripts" / name
        if not local.is_file():
            print(f"SKIP missing script: {name}")
            continue
        remote = f"{REMOTE_BASE}/scripts/{name}"
        sftp.put(str(local), remote)
        results.append((remote, sha256_file(local)))
        print(f"UP {name}")
    for name in UPLOAD_CONFIG:
        local = ROOT / "config" / name
        if not local.is_file():
            print(f"SKIP missing config: {name}")
            continue
        remote = f"{REMOTE_BASE}/config/{name}"
        sftp.put(str(local), remote)
        results.append((remote, sha256_file(local)))
        print(f"UP config/{name}")
    for name in UPLOAD_ROOT_FILES:
        local = ROOT / name
        if not local.is_file():
            print(f"SKIP missing root file: {name}")
            continue
        remote = f"{REMOTE_BASE}/{name}"
        sftp.put(str(local), remote)
        results.append((remote, sha256_file(local)))
        print(f"UP {name}")
    ndx = ROOT / "builds" / "analysis.ndx"
    sftp.put(str(ndx), f"{REMOTE_BASE}/builds/analysis.ndx")
    results.append((f"{REMOTE_BASE}/builds/analysis.ndx", sha256_file(ndx)))
    print("UP builds/analysis.ndx")
    mismatches = []
    for remote, expected in results:
        observed = sha256_remote_file(sftp, remote)
        if observed != expected:
            mismatches.append((remote, expected, observed))
    if mismatches:
        for remote, expected, observed in mismatches:
            print(f"HASH MISMATCH {remote}: expected={expected} observed={observed}")
        sftp.close()
        c.close()
        return 1
    sftp.close()
    c.close()
    print(f"\nUploaded and SHA-256 verified {len(results)} files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
