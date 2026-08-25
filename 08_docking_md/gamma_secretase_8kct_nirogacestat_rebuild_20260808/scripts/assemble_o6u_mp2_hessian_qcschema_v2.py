#!/usr/bin/env python3
"""Assemble and analyze the frozen O6U MP2 finite-difference Hessian."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from run_o6u_mp2_optimization_canary import sha256


PSI4_OPTIONS = {
    "basis": "6-31G(d)",
    "reference": "rhf",
    "scf_type": "df",
    "mp2_type": "df",
    "freeze_core": True,
    "guess": "sad",
    "e_convergence": 1.0e-8,
    "d_convergence": 1.0e-8,
    "maxiter": 200,
    "points": 3,
    "disp_size": 0.005,
}


def canonical(obj: object) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--manifest-sha256", required=True)
    ap.add_argument("--state", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    manifest_path = args.manifest.resolve()
    state_path = args.state.resolve()
    out = args.output_dir.resolve()
    if out.exists():
        raise SystemExit(f"Refusing to reuse output directory: {out}")
    if sha256(manifest_path) != args.manifest_sha256:
        raise SystemExit("Manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "pass_parallel_qcschema_plan_generated_no_qm_executed":
        raise SystemExit("Manifest is not the frozen no-QM plan")
    if state.get("status") != "pass_all_qcschema_gradients_validated":
        raise SystemExit("Parallel state has not passed all worker validations")
    if state.get("pass_count") != manifest.get("task_count") or state.get("failed_count") != 0:
        raise SystemExit("Parallel state count mismatch")
    if state.get("manifest", {}).get("sha256") != args.manifest_sha256:
        raise SystemExit("Parallel state does not bind the frozen manifest")
    tasks_meta = manifest["tasks"]
    if len(tasks_meta) != 445:
        raise SystemExit("Unexpected task count")

    from qcelemental.models import AtomicResult
    import psi4
    from psi4.driver import driver

    psi4.core.be_quiet()
    psi4.set_num_threads(1)
    psi4.set_memory("4 GiB")
    psi4.set_options(PSI4_OPTIONS)

    ref_input = json.loads(Path(tasks_meta[0]["input"]).read_text(encoding="utf-8"))
    ref_output = json.loads(Path(tasks_meta[0]["output"]).read_text(encoding="utf-8"))
    if tasks_meta[0]["label"] != "reference" or ref_output.get("success") is not True:
        raise SystemExit("Reference task is not valid")
    molecule = psi4.core.Molecule.from_schema(ref_input["molecule"])
    ref_gradient = np.asarray(ref_output["return_result"], dtype=float).reshape((-1, 3))
    if ref_gradient.shape != (76, 3) or not np.isfinite(ref_gradient).all():
        raise SystemExit("Reference gradient is invalid")
    plan = psi4.hessian(
        "mp2",
        molecule=molecule,
        dertype="gradient",
        return_plan=True,
        ref_gradient=psi4.core.Matrix.from_array(ref_gradient),
    )
    if list(plan.task_list) != [task["label"] for task in tasks_meta]:
        raise SystemExit("Reconstructed finite-difference labels/order differ")

    result_hashes = []
    for meta, (label, task) in zip(tasks_meta, plan.task_list.items(), strict=True):
        input_path = Path(meta["input"])
        output_path = Path(meta["output"])
        if sha256(input_path) != meta["input_sha256"]:
            raise SystemExit(f"Input hash mismatch: {label}")
        stored_input = json.loads(input_path.read_text(encoding="utf-8"))
        reconstructed_input = json.loads(task.plan().json())
        if canonical(stored_input) != canonical(reconstructed_input):
            raise SystemExit(f"Reconstructed QCSchema input differs: {label}")
        result_obj = AtomicResult.parse_obj(json.loads(output_path.read_text(encoding="utf-8")))
        if not result_obj.success or result_obj.driver.value != "gradient":
            raise SystemExit(f"Invalid result: {label}")
        grad = np.asarray(result_obj.return_result, dtype=float).reshape((-1, 3))
        if grad.shape != (76, 3) or not np.isfinite(grad).all():
            raise SystemExit(f"Invalid gradient: {label}")
        task.result = result_obj
        task.computed = True
        result_hashes.append({"label": label, "path": str(output_path), "sha256": sha256(output_path), "bytes": output_path.stat().st_size})

    out.mkdir(parents=True, exist_ok=False)
    analysis_log = out / "frequency_analysis.psi4.out"
    psi4.core.set_output_file(str(analysis_log), False)
    hessian_core, wfn = plan.get_psi_results(return_wfn=True)
    hessian = np.asarray(hessian_core, dtype=float)
    if hessian.shape != (228, 228) or not np.isfinite(hessian).all():
        raise SystemExit(f"Invalid assembled Hessian shape/content: {hessian.shape}")
    symmetry_error = float(np.max(np.abs(hessian - hessian.T)))
    psi4.core.set_variable("CURRENT ENERGY", float(ref_output["properties"]["return_energy"]))
    vibinfo = driver.vibanal_wfn(wfn, hess=hessian, molecule=molecule, project_trans=True, project_rot=True)
    psi4.core.flush_outfile()

    hessian_npy = out / "o6u_mp2_631gd_hessian.npy"
    hessian_txt = out / "o6u_mp2_631gd_hessian_hartree_per_bohr2.txt"
    np.save(hessian_npy, hessian, allow_pickle=False)
    np.savetxt(hessian_txt, hessian, fmt="%.16e")
    vib_json = out / "o6u_mp2_631gd_vibrational_analysis.json"
    vib_payload = {key: json.loads(value.json()) for key, value in vibinfo.items()}
    vib_json.write_text(json.dumps(vib_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

    omega = np.asarray(vibinfo["omega"].data, dtype=complex).reshape(-1)
    raw_frequencies = [{"real_cm-1": float(z.real), "imaginary_cm-1": float(abs(z.imag))} for z in omega]
    imaginary = [item for item in raw_frequencies if item["imaginary_cm-1"] > 0.0]
    if not imaginary:
        status = "pass_mp2_minimum_no_imaginary_modes"
        downstream_release = True
    else:
        status = "pending_imaginary_mode_scientific_audit"
        downstream_release = False

    report = {
        "schema_version": "2.0",
        "report_type": "o6u_mp2_hessian_minimum_character",
        "status": status,
        "production_approved": downstream_release,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": {"path": str(manifest_path), "sha256": args.manifest_sha256},
        "parallel_state": {"path": str(state_path), "sha256": sha256(state_path)},
        "task_count": len(tasks_meta),
        "result_hashes": result_hashes,
        "method": manifest["method"],
        "hessian_shape": list(hessian.shape),
        "hessian_symmetry_max_abs_hartree_per_bohr2": symmetry_error,
        "frequency_count": len(raw_frequencies),
        "frequencies": raw_frequencies,
        "imaginary_mode_count": len(imaginary),
        "imaginary_modes": imaginary,
        "decision_rule": "Pass automatically only with zero raw imaginary modes. Any imaginary mode remains pending explicit numerical and mode-shape audit; no arbitrary magnitude cutoff is applied here.",
        "artifacts": {
            "hessian_npy": {"path": str(hessian_npy), "sha256": sha256(hessian_npy)},
            "hessian_text": {"path": str(hessian_txt), "sha256": sha256(hessian_txt)},
            "vibrational_analysis": {"path": str(vib_json), "sha256": sha256(vib_json)},
            "analysis_log": {"path": str(analysis_log), "sha256": sha256(analysis_log)},
        },
        "compound_context": "This validates local minimum character only for the identity-, chirality-, charge-, and multiplicity-locked nirogacestat parameterization geometry.",
        "disease_context": "The result provides no SSc efficacy evidence and cannot replace disease-specific HES1-Notch validation.",
        "release_boundary": "Parameter fitting may proceed only when this report passes or a separately documented imaginary-mode audit authorizes a scientifically justified recovery.",
    }
    report_path = out / "O6U_MP2_HESSIAN_MINIMUM_CHARACTER.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": status, "report": str(report_path), "sha256": sha256(report_path), "imaginary_mode_count": len(imaginary)}, sort_keys=True))
    return 0 if downstream_release else 2


if __name__ == "__main__":
    raise SystemExit(main())
