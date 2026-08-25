#!/usr/bin/env python3
"""Create fixed-window, deterministic PBC-safe analysis trajectories.

The only allowed window is 200-500 ns from all three fixed 500 ns runs. No
operation removes, smooths, interpolates, or conditionally filters frames
inside that window. Any failed realization makes the analysis inconclusive.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from md_contract import artifact_path, canonical_json_sha256
from validate_qc_stationarity_report import (
    REQUIRED_CRITERIA,
    seal_report,
    validate_contract as validate_qc_report_contract,
    write_synthetic_criterion_evidence,
)


EXPECTED_STEPS = [
    "whole",
    "cluster_complex_if_required",
    "processed_first_frame_reference",
    "nojump",
    "center_and_rebox",
    "fit_analysis_selection",
    "fixed_window_extract",
]
SYSTEM_ID = "8kct_nirogacestat_native"
REALIZATION_IDS = ("rep01", "rep02", "rep03")


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def inside(root: Path, value: str, *, must_exist: bool = True) -> Path:
    path = (root / value).resolve()
    if path != root and root not in path.parents:
        raise ValueError(f"Path escapes package root: {value}")
    if must_exist and not path.exists():
        raise FileNotFoundError(path)
    return path


def parse_ndx(path: Path) -> dict[str, int]:
    groups: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            name = line[1:-1].strip()
            if not name or name in groups:
                raise ValueError(f"Empty or duplicate index group in {path}: {name!r}")
            groups[name] = len(groups)
    if not groups:
        raise ValueError(f"No index groups found: {path}")
    return groups


def command(
    gmx: str,
    tpr: Path,
    source: Path,
    output: Path,
    ndx: Path,
    mode: str,
    window_ns: tuple[float, float],
    processed_reference: Path | None = None,
) -> list[str]:
    if mode == "fixed_window_extract":
        return [
            gmx,
            "trjconv",
            "-f",
            str(source),
            "-o",
            str(output),
            "-b",
            f"{window_ns[0] * 1000.0:.3f}",
            "-e",
            f"{window_ns[1] * 1000.0:.3f}",
        ]
    structure = processed_reference if mode == "nojump" else tpr
    base = [gmx, "trjconv", "-s", str(structure), "-f", str(source), "-o", str(output), "-n", str(ndx)]
    if mode == "whole":
        return base + ["-pbc", "whole"]
    if mode == "cluster_complex_if_required":
        return base + ["-pbc", "cluster"]
    if mode == "processed_first_frame_reference":
        return base + ["-dump", "0"]
    if mode == "nojump":
        if processed_reference is None:
            raise ValueError("nojump requires the processed first-frame reference")
        return base + ["-pbc", "nojump"]
    if mode == "center_and_rebox":
        return base + ["-pbc", "mol", "-center", "-ur", "compact"]
    if mode == "fit_analysis_selection":
        return base + ["-fit", "rot+trans"]
    raise ValueError(mode)


def execute(argv: list[str], selections: list[int], cwd: Path, record_path: Path) -> None:
    selection_text = "".join(f"{value}\n" for value in selections)
    result = subprocess.run(
        argv,
        cwd=cwd,
        input=selection_text,
        text=True,
        capture_output=True,
        errors="replace",
        check=False,
    )
    record = {
        "argv": argv,
        "cwd": str(cwd),
        "selections": selections,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "finished_at_utc": utc_now(),
    }
    record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"Trajectory command failed; see {record_path}")


def parse_xvg_distance(path: Path) -> tuple[list[float], list[float]]:
    times: list[float] = []
    distances: list[float] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="strict").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"{path}:{line_number}: expected time and distance")
        time_ps = float(fields[0]); distance_nm = float(fields[1])
        if not (math.isfinite(time_ps) and math.isfinite(distance_nm)) or distance_nm < 0:
            raise ValueError(f"{path}:{line_number}: non-finite/negative distance row")
        if times and time_ps <= times[-1]:
            raise ValueError(f"{path}:{line_number}: times are not strictly increasing")
        times.append(time_ps); distances.append(distance_nm)
    if not times:
        raise ValueError(f"No distance rows found in {path}")
    return times, distances


def compare_minimum_image_distances(raw_xvg: Path, processed_xvg: Path, tolerance_nm: float) -> dict[str, Any]:
    raw_times, raw_distances = parse_xvg_distance(raw_xvg)
    processed_times, processed_distances = parse_xvg_distance(processed_xvg)
    errors: list[str] = []
    if raw_times != processed_times:
        errors.append("raw and processed minimum-distance series do not contain exactly matching frame times")
    differences: list[float] = []
    if not errors:
        differences = [abs(left - right) for left, right in zip(raw_distances, processed_distances)]
        if any(value > tolerance_nm for value in differences):
            errors.append(
                f"minimum-image protein-O6U heavy-atom distance changed by more than {tolerance_nm:.6f} nm"
            )
    return {
        "schema_version": "1.0",
        "metric": "matching_frame_minimum_image_protein_O6U_heavy_atom_distance",
        "frame_count": len(raw_times),
        "times_match_exactly": raw_times == processed_times,
        "time_ps_sha256": canonical_json_sha256(raw_times),
        "tolerance_nm": tolerance_nm,
        "maximum_absolute_difference_nm": max(differences) if differences else None,
        "status": "pass" if not errors else "fail",
        "errors": errors,
    }


def synthetic_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="single_system_trajectory_plan_") as temporary:
        root = Path(temporary)
        config = root / "config"
        reports = root / "reports"
        evidence_dir = root / "evidence"
        config.mkdir(); reports.mkdir(); evidence_dir.mkdir()
        index_path = root / "builds" / "analysis.ndx"
        index_path.parent.mkdir()
        index_path.write_text(
            "[ System ]\n1\n[ Protein_O6U ]\n1\n[ PSEN1_Core ]\n1\n"
            "[ Protein_Heavy ]\n1\n[ O6U_Heavy ]\n1\n",
            encoding="utf-8",
        )
        realizations: list[dict[str, Any]] = []
        for rid in REALIZATION_IDS:
            run = root / "runs" / SYSTEM_ID / rid
            work = run / "work"
            work.mkdir(parents=True)
            (work / "production.tpr").write_bytes(b"synthetic tpr")
            (work / "production.xtc").write_bytes(b"synthetic xtc")
            realizations.append({"id": rid, "run_directory": str(run.relative_to(root))})
        build_report = reports / "build.json"; build_report.write_text("{}\n", encoding="utf-8")
        manifest = {
            "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
            "simulation": {
                "analysis_window_ns": [200.0, 500.0],
            },
            "analysis": {"pbc_distance_invariance_tolerance_nm": 0.01},
            "systems": [{
                "id": SYSTEM_ID,
                "construction": {
                    "analysis_index": {"path": str(index_path.relative_to(root)).replace("\\", "/"), "sha256": sha256(index_path)},
                    "charmm_gui_archive": {"sha256": "a" * 64},
                    "build_validation_report": {"path": "reports/build.json", "sha256": sha256(build_report)},
                },
                "realizations": realizations,
            }],
        }
        manifest_path = config / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        plan = {
            "analysis_window_ns": [200.0, 500.0],
            "trajectory_processing": {
                "ordered_steps": EXPECTED_STEPS,
                "cluster_complex_required": True,
                "fixed_window_ns": [200.0, 500.0],
                "no_frame_removal_within_fixed_window": True,
                "no_smoothing": True,
                "no_interpolation": True,
                "pbc_distance_invariance_tolerance_nm": 0.01,
                "snapshot_times_ns": [200.0, 350.0, 500.0],
                "groups": {
                    "system": "System", "complex": "Protein_O6U", "fit": "PSEN1_Core",
                    "analysis_output": "Protein_O6U", "pbc_invariance_protein_heavy": "Protein_Heavy",
                    "pbc_invariance_ligand_heavy": "O6U_Heavy",
                },
            },
            "eligibility_gate": {
                "required_realization_ids": list(REALIZATION_IDS),
                "all_realizations_passed_qc_and_stationarity": True,
                "failure_policy": "inconclusive_if_any_realization_fails_qc_or_stationarity",
                "qc_and_stationarity_report_sha256": "BOUND_AFTER_SEALING",
            },
        }
        plan_path = config / "analysis_plan.json"
        plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        raw_runs = [{"realization_id": rid, "status": "pass", "artifacts": {}} for rid in REALIZATION_IDS]
        raw_report = {
            "schema_version": "2.0", "status": "pass", "phase": "production", "strict": True,
            "manifest_sha256": sha256(manifest_path),
            "construction_archive_sha256": "a" * 64, "runs": raw_runs,
        }
        raw_report_path = reports / "raw.json"
        raw_report_path.write_text(json.dumps(raw_report, indent=2) + "\n", encoding="utf-8")
        timestamp = "2026-08-08T00:00:00+00:00"
        results = []
        for rid, raw_run in zip(REALIZATION_IDS, raw_runs):
            evidence_rows = []
            for criterion in REQUIRED_CRITERIA:
                evidence_path = write_synthetic_criterion_evidence(
                    root, evidence_dir, rid, criterion, canonical_json_sha256(raw_run)
                )
                evidence_rows.append({
                    "criterion_id": criterion, "status": "pass",
                    "artifact": {"path": evidence_path.relative_to(root).as_posix(), "sha256": sha256(evidence_path)},
                    "evaluator_note": "Synthetic criterion evidence.",
                })
            results.append({
                "realization_id": rid, "raw_output_run_payload_sha256": canonical_json_sha256(raw_run),
                "qc_status": "pass", "stationarity_status": "pass",
                "evidence": evidence_rows,
            })
        qc_draft = {
            "schema_version": "1.0", "report_type": "production_qc_and_stationarity_attestation",
            "study_id": manifest["study_id"], "system_id": SYSTEM_ID, "construction_id": "build01",
            "status": "pass", "strict": True, "evaluated_at_utc": timestamp,
            "evaluator": {"name": "Test", "role": "test", "software_and_versions": ["synthetic"],
                          "declaration": "All prespecified evidence was reviewed without changing the window or selecting a realization."},
            "bindings": {
                "manifest": {"path": "config/manifest.json", "sha256": sha256(manifest_path)},
                "analysis_plan": {"path": "config/analysis_plan.json", "contract_sha256": "filled by plan"},
                "raw_output_validation_report": {"path": "reports/raw.json", "sha256": sha256(raw_report_path)},
                "build_validation_report": {"path": "reports/build.json", "sha256": sha256(build_report)},
                "construction_archive_sha256": "a" * 64, "validator_sha256": "UNSEALED",
            },
            "adjudication_contract": {
                "required_realization_ids": list(REALIZATION_IDS), "all_must_pass": True,
                "analysis_window_ns": [200.0, 500.0], "trajectory_exclusion_allowed": False,
                "analysis_cutoff_change_allowed": False, "outcome_dependent_extension_allowed": False,
                "ligand_behavior_used_for_replica_selection": False,
                "required_criterion_ids": list(REQUIRED_CRITERIA),
                "failure_policy": "inconclusive_if_any_realization_fails_qc_or_stationarity",
            },
            "results": results, "overall_decision": "pass",
            "approval": {"approval_status": "approved", "approver_name": "Approver", "approver_role": "test",
                         "approved_at_utc": timestamp, "signature_scheme": "sha256_canonical_json_checksum_attestation_v1",
                         "signed_payload_sha256": "UNSEALED"},
        }
        from md_contract import analysis_plan_contract_sha256
        qc_draft["bindings"]["analysis_plan"]["contract_sha256"] = analysis_plan_contract_sha256(plan)
        draft_path = reports / "qc_draft.json"; draft_path.write_text(json.dumps(qc_draft), encoding="utf-8")
        report_path = reports / "qc_sealed.json"
        seal_report(root, manifest_path, plan_path, raw_report_path, draft_path, report_path)
        plan["eligibility_gate"]["qc_and_stationarity_report_sha256"] = sha256(report_path)
        plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        command_line = [
            sys.executable, str(Path(__file__).resolve()), "--manifest", str(manifest_path),
            "--analysis-plan", str(plan_path), "--production-qc-report", str(report_path),
        ]
        passed = subprocess.run(command_line, capture_output=True, text=True, check=False)
        if passed.returncode != 0:
            raise RuntimeError(f"Synthetic trajectory dry run failed:\n{passed.stdout}\n{passed.stderr}")
        planned = json.loads(passed.stdout)
        if len(planned.get("runs", [])) != 3 or any(len(run.get("steps", [])) != 7 for run in planned["runs"]):
            raise RuntimeError("Synthetic trajectory dry run did not plan three complete pipelines")
        steps = planned["runs"][0]["steps"]
        fixed = steps[-1]["argv"]
        if "-b" not in fixed or fixed[fixed.index("-b") + 1] != "200000.000" or "-e" not in fixed or fixed[fixed.index("-e") + 1] != "500000.000":
            raise RuntimeError("Primary fixed-window command is incorrect")
        if steps[-1].get("selection_group_numbers") != [0] or "-s" in fixed or "-pbc" in fixed:
            raise RuntimeError("Fixed-window extraction must select the fitted trajectory group without another structural or PBC transform")
        fit_position = next(index for index, step in enumerate(steps) if step["mode"] == "fit_analysis_selection")
        if any("-pbc" in step["argv"] for step in steps[fit_position + 1:]):
            raise RuntimeError("A PBC operation was planned after coordinate fitting")
        if (root / "analysis").exists():
            raise RuntimeError("Plan-only mode mutated the analysis directory")
        pbc_plan = planned["runs"][0].get("pbc_distance_invariance", {})
        if len(pbc_plan.get("commands", [])) != 2 or pbc_plan.get("tolerance_nm") != 0.01:
            raise RuntimeError("Raw/processed PBC-invariance calculation was not planned")
        if pbc_plan.get("window_ns") != [0.0, 500.0]:
            raise RuntimeError("PBC invariance must cover the complete 0-500 ns trajectory")
        for pbc_command in pbc_plan["commands"]:
            if pbc_command[pbc_command.index("-b") + 1] != "0.000" or pbc_command[pbc_command.index("-e") + 1] != "500000.000":
                raise RuntimeError("PBC-invariance command does not cover exactly 0-500 ns")
        raw_xvg = reports / "raw.xvg"; processed_xvg = reports / "processed.xvg"
        raw_xvg.write_text("0 0.30\n100 0.40\n", encoding="utf-8")
        processed_xvg.write_text("0 0.309\n100 0.391\n", encoding="utf-8")
        if compare_minimum_image_distances(raw_xvg, processed_xvg, 0.01)["status"] != "pass":
            raise RuntimeError("A <=0.01 nm matching-frame PBC difference failed")
        processed_xvg.write_text("0 0.311\n100 0.40\n", encoding="utf-8")
        if compare_minimum_image_distances(raw_xvg, processed_xvg, 0.01)["status"] != "fail":
            raise RuntimeError("A >0.01 nm PBC difference did not fail")
        tampered = load_json(report_path)
        tampered["results"][0]["stationarity_status"] = "fail"
        report_path.write_text(json.dumps(tampered, indent=2) + "\n", encoding="utf-8")
        plan["eligibility_gate"]["qc_and_stationarity_report_sha256"] = sha256(report_path)
        plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        failed = subprocess.run(command_line, capture_output=True, text=True, check=False)
        if failed.returncode == 0 or "altered or is unsealed" not in (failed.stdout + failed.stderr):
            raise RuntimeError("Manually altered QC/stationarity report did not fail closed")
    print("SELF-TEST PASS: frozen PBC order/window, sealed three-run adjudication, exact raw/processed matching-frame 0.01-nm test, and no dry-run mutation.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--analysis-plan", type=Path)
    parser.add_argument("--production-qc-report", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--execute", action="store_true", help="Execute; default prints a plan")
    parser.add_argument("--self-test", action="store_true", help="Run a synthetic fail-closed trajectory-planning test")
    args = parser.parse_args()

    if args.self_test:
        if any(value is not None for value in (args.manifest, args.analysis_plan, args.production_qc_report)) or args.execute:
            raise SystemExit("--self-test cannot be combined with input or execution options")
        synthetic_self_test()
        return 0
    if args.manifest is None or args.analysis_plan is None or args.production_qc_report is None:
        parser.error("--manifest, --analysis-plan, and --production-qc-report are required unless --self-test is used")

    manifest_path = args.manifest.resolve()
    package_root = manifest_path.parent.parent.resolve()
    manifest = load_json(manifest_path)
    plan = load_json(args.analysis_plan.resolve())
    qc_report_path = args.production_qc_report.resolve()
    qc_report = load_json(qc_report_path)
    raw_binding = qc_report.get("bindings", {}).get("raw_output_validation_report", {})
    raw_output_report_path = inside(package_root, str(raw_binding.get("path", "")))
    qc_errors = validate_qc_report_contract(
        package_root,
        manifest_path,
        args.analysis_plan.resolve(),
        raw_output_report_path,
        qc_report_path,
    )
    if qc_errors:
        raise SystemExit(
            "QC/stationarity report validation failed; analysis is locked:\n" +
            "\n".join(f"- {error}" for error in qc_errors)
        )
    if qc_report.get("status") != "pass" or qc_report.get("overall_decision") != "pass":
        raise SystemExit("QC and stationarity were not passed by all three realizations; analysis is inconclusive")
    raw_output_report = load_json(raw_output_report_path)
    eligibility = plan.get("eligibility_gate", {})
    if eligibility.get("required_realization_ids") != list(REALIZATION_IDS):
        raise SystemExit("Analysis-plan eligibility gate must require rep01-rep03")
    if eligibility.get("all_realizations_passed_qc_and_stationarity") is not True:
        raise SystemExit("Analysis plan is not eligible because all three QC/stationarity results are not passed")
    if eligibility.get("failure_policy") != "inconclusive_if_any_realization_fails_qc_or_stationarity":
        raise SystemExit("Analysis-plan failure policy is not frozen correctly")
    if eligibility.get("qc_and_stationarity_report_sha256") != sha256(qc_report_path):
        raise SystemExit("Analysis plan is bound to another QC/stationarity report")
    processing = plan.get("trajectory_processing", {})
    if processing.get("ordered_steps") != EXPECTED_STEPS:
        raise SystemExit(f"trajectory_processing.ordered_steps must be {EXPECTED_STEPS}")
    for key in ("no_frame_removal_within_fixed_window", "no_smoothing", "no_interpolation"):
        if processing.get(key) is not True:
            raise SystemExit(f"trajectory_processing.{key} must be true")
    if processing.get("cluster_complex_required") is not True:
        raise SystemExit("cluster_complex_required must be frozen true for this membrane complex")
    manifest_tolerance = manifest.get("analysis", {}).get("pbc_distance_invariance_tolerance_nm")
    if manifest_tolerance != 0.01 or processing.get("pbc_distance_invariance_tolerance_nm") != manifest_tolerance:
        raise SystemExit("PBC distance-invariance tolerance must be 0.01 nm in both manifest and analysis plan")
    expected_window = manifest.get("simulation", {}).get("analysis_window_ns")
    window = processing.get("fixed_window_ns")
    if expected_window != [200.0, 500.0]:
        raise SystemExit("Manifest analysis window must be exactly 200-500 ns")
    if window != expected_window or plan.get("analysis_window_ns") != expected_window:
        raise SystemExit("Analysis plan does not match the fixed 200-500 ns window")
    window_ns = (float(window[0]), float(window[1]))
    expected_snapshots = [window_ns[0], (window_ns[0] + window_ns[1]) / 2.0, window_ns[1]]
    if processing.get("snapshot_times_ns") != expected_snapshots:
        raise SystemExit(f"snapshot_times_ns must be {expected_snapshots} for the selected window")
    if args.execute and shutil.which(args.gmx) is None:
        raise SystemExit(f"GROMACS executable not found: {args.gmx}")

    plans: list[dict[str, Any]] = []
    systems = manifest.get("systems", [])
    if not isinstance(systems, list) or len(systems) != 1 or systems[0].get("id") != SYSTEM_ID:
        raise SystemExit(f"Manifest must contain exactly one system: {SYSTEM_ID}")
    system = systems[0]
    construction = system.get("construction", {})
    realizations = system.get("realizations", [])
    if not isinstance(realizations, list) or [item.get("id") for item in realizations] != list(REALIZATION_IDS):
        raise SystemExit(f"Manifest realizations must be exactly {list(REALIZATION_IDS)}")
    raw_runs = raw_output_report.get("runs", [])
    passed = {item.get("realization_id") for item in raw_runs if item.get("status") == "pass"}
    if passed != set(REALIZATION_IDS):
        raise SystemExit("Production QC/stationarity report must pass rep01-rep03 together")
    group_names = processing.get("groups", {})
    required_group_keys = (
        "system", "complex", "fit", "analysis_output",
        "pbc_invariance_protein_heavy", "pbc_invariance_ligand_heavy",
    )
    if not all(isinstance(group_names.get(key), str) and "TODO" not in group_names[key] for key in required_group_keys):
        raise SystemExit(f"Unfrozen trajectory groups for {SYSTEM_ID}")
    ndx = artifact_path(package_root, construction.get("analysis_index"), "build01 analysis index")
    groups = parse_ndx(ndx)
    missing = [group_names[key] for key in required_group_keys if group_names[key] not in groups]
    if missing:
        raise SystemExit(f"Missing frozen index groups: {missing}")
    for realization in realizations:
        rid = str(realization.get("id"))
        run_dir = inside(package_root, str(realization["run_directory"]))
        work = run_dir / "work"
        tpr = work / "production.tpr"
        raw = work / "production.xtc"
        if not tpr.is_file() or not raw.is_file():
            raise FileNotFoundError(f"Missing production TPR/XTC for {SYSTEM_ID}/{rid}")
        output_dir = package_root / "analysis" / "trajectories" / SYSTEM_ID / rid
        paths = {
            "whole": output_dir / "01_whole.xtc",
            "cluster_complex_if_required": output_dir / "02_clustered_complex.xtc",
            "processed_first_frame_reference": output_dir / "03_processed_first_frame.gro",
            "nojump": output_dir / "04_nojump.xtc",
            "center_and_rebox": output_dir / "05_centered_reboxed.xtc",
            "fit_analysis_selection": output_dir / "06_fitted_analysis.xtc",
            "fixed_window_extract": output_dir / f"07_fixed_{window_ns[0]:g}_{window_ns[1]:g}ns.xtc",
        }
        invariance_paths = {
            "raw_xvg": output_dir / "09_raw_minimum_image_protein_O6U_heavy.xvg",
            "processed_xvg": output_dir / "10_processed_minimum_image_protein_O6U_heavy.xvg",
            "report": output_dir / "11_pbc_distance_invariance.json",
        }
        steps = [
            ("whole", raw, paths["whole"], [groups[group_names["system"]]]),
            ("cluster_complex_if_required", paths["whole"], paths["cluster_complex_if_required"], [groups[group_names["complex"]], groups[group_names["system"]]]),
            ("processed_first_frame_reference", paths["cluster_complex_if_required"], paths["processed_first_frame_reference"], [groups[group_names["system"]]]),
            ("nojump", paths["cluster_complex_if_required"], paths["nojump"], [groups[group_names["system"]]]),
            ("center_and_rebox", paths["nojump"], paths["center_and_rebox"], [groups[group_names["complex"]], groups[group_names["system"]]]),
            ("fit_analysis_selection", paths["center_and_rebox"], paths["fit_analysis_selection"], [groups[group_names["fit"]], groups[group_names["analysis_output"]]]),
            # The fitted trajectory already contains only the frozen analysis
            # selection.  With no topology supplied, trjconv exposes that
            # trajectory as group 0; select it explicitly so execution cannot
            # block or fail at an unanswered interactive prompt.
            ("fixed_window_extract", paths["fit_analysis_selection"], paths["fixed_window_extract"], [0]),
        ]
        run_plan: list[dict[str, Any]] = []
        for index, (mode, source, output, selections) in enumerate(steps, start=1):
            argv = command(
                args.gmx,
                tpr,
                source,
                output,
                ndx,
                mode,
                window_ns,
                paths["processed_first_frame_reference"],
            )
            run_plan.append({"mode": mode, "argv": argv, "selection_group_numbers": selections})
            if args.execute:
                output_dir.mkdir(parents=True, exist_ok=True)
                if output.exists():
                    raise SystemExit(f"Refusing to overwrite derived trajectory: {output}")
                execute(argv, selections, output_dir, output_dir / f"{index:02d}_{mode}.command.json")
        pbc_selections = [
            groups[group_names["pbc_invariance_protein_heavy"]],
            groups[group_names["pbc_invariance_ligand_heavy"]],
        ]
        pbc_commands = [
            [
                args.gmx, "mindist", "-s", str(tpr), "-f", str(raw), "-n", str(ndx),
                "-od", str(invariance_paths["raw_xvg"]), "-b", "0.000",
                "-e", "500000.000",
            ],
            [
                args.gmx, "mindist", "-s", str(tpr), "-f", str(paths["center_and_rebox"]), "-n", str(ndx),
                "-od", str(invariance_paths["processed_xvg"]), "-b", "0.000",
                "-e", "500000.000",
            ],
        ]
        pbc_result: dict[str, Any] | None = None
        if args.execute:
            for output in invariance_paths.values():
                if output.exists():
                    raise SystemExit(f"Refusing to overwrite PBC-invariance artifact: {output}")
            execute(pbc_commands[0], pbc_selections, output_dir, output_dir / "09_raw_mindist.command.json")
            execute(pbc_commands[1], pbc_selections, output_dir, output_dir / "10_processed_mindist.command.json")
            pbc_result = compare_minimum_image_distances(
                invariance_paths["raw_xvg"], invariance_paths["processed_xvg"], float(manifest_tolerance)
            )
            pbc_result.update({
                "system_id": SYSTEM_ID,
                "construction_id": "build01",
                "realization_id": rid,
                "analysis_window_ns": list(window_ns),
                "pbc_invariance_window_ns": [0.0, 500.0],
                "manifest_sha256": sha256(manifest_path),
                "analysis_plan_sha256": sha256(args.analysis_plan.resolve()),
                "raw_tpr": {"path": str(tpr), "sha256": sha256(tpr)},
                "raw_trajectory": {"path": str(raw), "sha256": sha256(raw)},
                "processed_trajectory": {
                    "path": str(paths["center_and_rebox"]),
                    "sha256": sha256(paths["center_and_rebox"]),
                },
                "analysis_index": {"path": str(ndx), "sha256": sha256(ndx)},
                "distance_series": {
                    "raw": {"path": str(invariance_paths["raw_xvg"]), "sha256": sha256(invariance_paths["raw_xvg"])},
                    "processed": {"path": str(invariance_paths["processed_xvg"]), "sha256": sha256(invariance_paths["processed_xvg"])},
                },
                "created_at_utc": utc_now(),
            })
            invariance_paths["report"].write_text(json.dumps(pbc_result, indent=2) + "\n", encoding="utf-8")
            if pbc_result["status"] != "pass":
                raise SystemExit(
                    f"PBC minimum-image distance invariance failed for {SYSTEM_ID}/{rid}; retained report: "
                    f"{invariance_paths['report']}"
                )
        if args.execute:
            integrity = subprocess.run(
                [args.gmx, "check", "-f", str(paths["fixed_window_extract"])],
                cwd=output_dir,
                capture_output=True,
                text=True,
                errors="replace",
                check=False,
            )
            integrity_record = {
                "argv": [args.gmx, "check", "-f", str(paths["fixed_window_extract"])],
                "returncode": integrity.returncode,
                "stdout": integrity.stdout,
                "stderr": integrity.stderr,
            }
            (output_dir / "08_gmx_check.json").write_text(json.dumps(integrity_record, indent=2) + "\n", encoding="utf-8")
            if integrity.returncode != 0:
                raise SystemExit(f"gmx check failed for {SYSTEM_ID}/{rid}")
            provenance = {
                "system_id": SYSTEM_ID,
                "construction_id": "build01",
                "realization_id": rid,
                "analysis_window_ns": list(window_ns),
                "pbc_distance_invariance_tolerance_nm": manifest_tolerance,
                "manifest_sha256": sha256(manifest_path),
                "analysis_plan_sha256": sha256(args.analysis_plan.resolve()),
                "raw_tpr_sha256": sha256(tpr),
                "raw_trajectory_sha256": sha256(raw),
                "analysis_index_sha256": sha256(ndx),
                "outputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
                "pbc_distance_invariance_report": {
                    "path": str(invariance_paths["report"]),
                    "sha256": sha256(invariance_paths["report"]),
                },
                "created_at_utc": utc_now(),
            }
            (output_dir / "trajectory_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        plans.append({
            "system_id": SYSTEM_ID,
            "realization_id": rid,
            "analysis_window_ns": list(window_ns),
            "pbc_distance_invariance_tolerance_nm": manifest_tolerance,
            "steps": run_plan,
            "pbc_distance_invariance": {
                "metric": "matching_frame_minimum_image_protein_O6U_heavy_atom_distance",
                "window_ns": [0.0, 500.0],
                "selection_group_numbers": pbc_selections,
                "commands": pbc_commands,
                "tolerance_nm": manifest_tolerance,
                "report": str(invariance_paths["report"]),
            },
        })
    if not plans:
        raise SystemExit("No production realizations found")
    print(json.dumps({"mode": "execute" if args.execute else "plan_only", "runs": plans}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
