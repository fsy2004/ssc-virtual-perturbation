#!/usr/bin/env python3
"""Seal and validate the evidence-bound all-realization QC/stationarity report."""

from __future__ import annotations

import argparse
import copy
import json
import re
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from md_contract import (
    REALIZATION_IDS,
    SYSTEM_ID,
    analysis_plan_contract_sha256,
    artifact_path,
    canonical_json_sha256,
    load_json,
    report_payload_sha256,
    sha256,
)


TODO_RE = re.compile(r"\b(?:TODO|TBD|PENDING|PLACEHOLDER)\b", re.IGNORECASE)
REQUIRED_CRITERIA = (
    "structural_qc_and_stationarity",
    "membrane_qc_and_stationarity",
    "thermodynamic_qc_and_stationarity",
    "spike_adjudication_completeness",
)
REQUIRED_CHECKS = {
    "structural_qc_and_stationarity": (
        "input_times_and_frames_complete", "five_fixed_200_500ns_blocks_evaluated",
        "linear_change_gate", "first_last_shift_gate", "adjacent_block_shift_gate", "change_point_shift_gate",
    ),
    "membrane_qc_and_stationarity": (
        "input_times_and_frames_complete", "five_fixed_200_500ns_blocks_evaluated",
        "linear_change_gate", "first_last_shift_gate", "adjacent_block_shift_gate", "change_point_shift_gate",
    ),
    "thermodynamic_qc_and_stationarity": (
        "input_times_and_frames_complete", "five_fixed_200_500ns_blocks_evaluated",
        "linear_change_gate", "first_last_shift_gate", "adjacent_block_shift_gate", "change_point_shift_gate",
    ),
    "spike_adjudication_completeness": (
        "all_frames_retained", "all_detected_flags_source_hashed",
        "all_detected_flags_adjudicated", "no_unresolved_technical_flags",
    ),
}


def contains_todo(value: Any) -> bool:
    if isinstance(value, str):
        return bool(TODO_RE.search(value))
    if isinstance(value, dict):
        return any(contains_todo(item) for item in value.values())
    if isinstance(value, list):
        return any(contains_todo(item) for item in value)
    return False


def parse_utc(value: Any, label: str, errors: list[str]) -> None:
    if not isinstance(value, str) or not value.endswith(("Z", "+00:00")):
        errors.append(f"{label} must be an explicit UTC ISO-8601 timestamp")
        return
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        errors.append(f"{label} is not valid ISO-8601")


def resolve_evidence(package_root: Path, artifact: Any, label: str, errors: list[str]) -> Path | None:
    try:
        return artifact_path(package_root, artifact, label)
    except (ValueError, FileNotFoundError) as exc:
        errors.append(str(exc))
        return None


def validate_criterion_evidence(
    package_root: Path,
    artifact: Any,
    rid: str,
    criterion: str,
    expected_status: str,
    raw_payload_sha256: str,
    errors: list[str],
) -> dict[str, str] | None:
    path = resolve_evidence(package_root, artifact, f"{rid} {criterion} evidence", errors)
    if path is None:
        return None
    try:
        evidence = load_json(path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"cannot parse {rid} {criterion} evidence: {exc}")
        return None
    expected_identity = {
        "schema_version": "1.0",
        "report_type": "prespecified_qc_stationarity_criterion_evidence",
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "system_id": SYSTEM_ID,
        "construction_id": "build01",
        "realization_id": rid,
        "criterion_id": criterion,
        "status": expected_status,
        "analysis_window_ns": [200.0, 500.0],
        "raw_output_run_payload_sha256": raw_payload_sha256,
        "all_input_frames_retained": True,
        "smoothing_applied": False,
        "interpolation_applied": False,
        "window_changed_after_review": False,
    }
    for key, expected in expected_identity.items():
        if evidence.get(key) != expected:
            errors.append(f"{rid} {criterion} evidence mismatch at {key}")
    contract_path = resolve_evidence(package_root, evidence.get("criterion_contract"), f"{rid} {criterion} contract", errors)
    contract_rules: dict[str, str] = {}
    if contract_path is not None:
        try:
            contract = load_json(contract_path)
            expected_contract_identity = {
                "schema_version": "1.0",
                "report_type": "prespecified_qc_stationarity_criterion_contract",
                "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
                "system_id": SYSTEM_ID,
                "criterion_id": criterion,
                "analysis_window_ns": [200.0, 500.0],
                "frozen_before_production_review": True,
                "approval_status": "approved",
            }
            for key, expected in expected_contract_identity.items():
                if contract.get(key) != expected:
                    errors.append(f"{rid} {criterion} contract mismatch at {key}")
            if contains_todo(contract):
                errors.append(f"{rid} {criterion} contract contains unresolved placeholders")
            parse_utc(contract.get("approved_at_utc"), f"{rid} {criterion} contract approved_at_utc", errors)
            rules = contract.get("check_rules")
            if not isinstance(rules, list) or [item.get("check_id") for item in rules if isinstance(item, dict)] != list(REQUIRED_CHECKS[criterion]):
                errors.append(f"{rid} {criterion} contract check rules differ from the frozen criterion schema")
            else:
                for item in rules:
                    rule = item.get("rule")
                    if not isinstance(rule, str) or not rule.strip() or contains_todo(rule):
                        errors.append(f"{rid} {criterion} contract rule {item.get('check_id')} is unresolved")
                    else:
                        contract_rules[item["check_id"]] = rule
            contract_integrity = contract.get("integrity")
            if not isinstance(contract_integrity, dict) or contract_integrity.get("payload_sha256") != report_payload_sha256(
                contract, ("integrity", "payload_sha256")
            ):
                errors.append(f"{rid} {criterion} contract payload checksum is invalid")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"cannot parse {rid} {criterion} contract: {exc}")
    sources = evidence.get("source_artifacts")
    if not isinstance(sources, list) or not sources:
        errors.append(f"{rid} {criterion} evidence lacks source artifacts")
    else:
        for index, source in enumerate(sources):
            resolve_evidence(package_root, source, f"{rid} {criterion} source {index}", errors)
    checks = evidence.get("checks")
    expected_checks = REQUIRED_CHECKS[criterion]
    if not isinstance(checks, list) or [item.get("check_id") for item in checks if isinstance(item, dict)] != list(expected_checks):
        errors.append(f"{rid} {criterion} checks must be exactly {list(expected_checks)}")
    else:
        for check in checks:
            if check.get("status") not in {"pass", "fail"}:
                errors.append(f"{rid} {criterion} check {check.get('check_id')} has invalid status")
            if "observed" not in check or contains_todo(check.get("observed")):
                errors.append(f"{rid} {criterion} check {check.get('check_id')} lacks an observed result")
            if not isinstance(check.get("rule"), str) or not check["rule"].strip() or contains_todo(check["rule"]):
                errors.append(f"{rid} {criterion} check {check.get('check_id')} lacks its frozen rule")
            elif contract_rules.get(check.get("check_id")) != check.get("rule"):
                errors.append(f"{rid} {criterion} check {check.get('check_id')} rule differs from its approved contract")
        computed_status = "pass" if all(check.get("status") == "pass" for check in checks) else "fail"
        if computed_status != expected_status:
            errors.append(f"{rid} {criterion} outer status differs from parsed checks")
    outcomes = evidence.get("criterion_outcomes")
    if not isinstance(outcomes, dict) or outcomes.get("qc_status") not in {"pass", "fail"}:
        errors.append(f"{rid} {criterion} lacks parsed QC outcome")
        outcomes = None
    elif criterion == "spike_adjudication_completeness":
        if outcomes.get("stationarity_status") != "not_applicable":
            errors.append(f"{rid} spike adjudication stationarity outcome must be not_applicable")
    elif outcomes.get("stationarity_status") not in {"pass", "fail"}:
        errors.append(f"{rid} {criterion} lacks parsed stationarity outcome")
    if isinstance(outcomes, dict):
        outcome_status = (
            "pass" if outcomes.get("qc_status") == "pass" and (
                criterion == "spike_adjudication_completeness" or outcomes.get("stationarity_status") == "pass"
            ) else "fail"
        )
        if outcome_status != expected_status:
            errors.append(f"{rid} {criterion} status differs from parsed QC/stationarity outcomes")
    integrity = evidence.get("integrity")
    if not isinstance(integrity, dict) or integrity.get("payload_sha256") != report_payload_sha256(
        evidence, ("integrity", "payload_sha256")
    ):
        errors.append(f"{rid} {criterion} evidence payload checksum is invalid")
    if contract_path is None:
        return None
    return outcomes


def validate_contract(
    package_root: Path,
    manifest_path: Path,
    plan_path: Path,
    raw_report_path: Path,
    report_path: Path,
    *,
    allow_unsealed: bool = False,
) -> list[str]:
    errors: list[str] = []
    manifest = load_json(manifest_path)
    plan = load_json(plan_path)
    raw_report = load_json(raw_report_path)
    report = load_json(report_path)
    if report.get("schema_version") != "1.0" or report.get("report_type") != "production_qc_and_stationarity_attestation":
        errors.append("QC/stationarity report schema or report_type is invalid")
    if report.get("study_id") != manifest.get("study_id"):
        errors.append("QC/stationarity report study_id differs from manifest")
    if report.get("system_id") != SYSTEM_ID or report.get("construction_id") != "build01":
        errors.append("QC/stationarity report system/build identity is invalid")
    if report.get("strict") is not True:
        errors.append("QC/stationarity report must be strict")
    if contains_todo(report):
        errors.append("QC/stationarity report contains unresolved TODO/TBD values")
    parse_utc(report.get("evaluated_at_utc"), "evaluated_at_utc", errors)
    evaluator = report.get("evaluator", {})
    for key in ("name", "role", "declaration"):
        if not isinstance(evaluator.get(key), str) or not evaluator[key].strip():
            errors.append(f"evaluator.{key} is unresolved")
    if not isinstance(evaluator.get("software_and_versions"), list) or not evaluator["software_and_versions"]:
        errors.append("evaluator.software_and_versions is unresolved")

    bindings = report.get("bindings", {})
    expected_manifest = {"path": manifest_path.relative_to(package_root).as_posix(), "sha256": sha256(manifest_path)}
    if bindings.get("manifest") != expected_manifest:
        errors.append("report is stale or bound to another manifest")
    expected_plan = {
        "path": plan_path.relative_to(package_root).as_posix(),
        "contract_sha256": analysis_plan_contract_sha256(plan),
    }
    if bindings.get("analysis_plan") != expected_plan:
        errors.append("report is stale or bound to another analysis-plan contract")
    expected_raw = {"path": raw_report_path.relative_to(package_root).as_posix(), "sha256": sha256(raw_report_path)}
    if bindings.get("raw_output_validation_report") != expected_raw:
        errors.append("report is stale or bound to another raw-output validation report")
    validator_hash = sha256(Path(__file__).resolve())
    if bindings.get("validator_sha256") != validator_hash:
        errors.append("report was sealed by a different validator version")

    systems = manifest.get("systems", [])
    construction = systems[0].get("construction", {}) if isinstance(systems, list) and len(systems) == 1 else {}
    build_artifact = construction.get("build_validation_report", {})
    try:
        build_path = artifact_path(package_root, build_artifact, "manifest build-validation report")
        expected_build = {"path": build_path.relative_to(package_root).as_posix(), "sha256": sha256(build_path)}
        if bindings.get("build_validation_report") != expected_build:
            errors.append("QC report is stale or bound to another build-validation report")
    except (ValueError, FileNotFoundError) as exc:
        errors.append(str(exc))
    archive_sha = construction.get("charmm_gui_archive", {}).get("sha256")
    if bindings.get("construction_archive_sha256") != archive_sha:
        errors.append("QC report archive SHA-256 differs from manifest")

    if raw_report.get("schema_version") not in {"1.0", "2.0"} or raw_report.get("status") != "pass":
        errors.append("bound raw-output report is not a passing recognized report")
    if raw_report.get("phase") != "production" or raw_report.get("strict") is not True:
        errors.append("bound raw-output report is not the strict completed-production audit")
    if raw_report.get("manifest_sha256") != sha256(manifest_path):
        errors.append("raw-output report is stale relative to manifest")
    if raw_report.get("construction_archive_sha256") != archive_sha:
        errors.append("raw-output report archive SHA-256 differs from manifest")
    raw_runs = raw_report.get("runs", [])
    if not isinstance(raw_runs, list) or [item.get("realization_id") for item in raw_runs] != list(REALIZATION_IDS):
        errors.append("raw-output report must contain rep01-rep03 in order")
        raw_by_id: dict[str, Any] = {}
    else:
        raw_by_id = {item["realization_id"]: item for item in raw_runs}
        if any(item.get("status") != "pass" for item in raw_runs):
            errors.append("all raw-output realization audits must pass")

    contract = report.get("adjudication_contract", {})
    expected_contract = {
        "required_realization_ids": list(REALIZATION_IDS),
        "all_must_pass": True,
        "analysis_window_ns": [200.0, 500.0],
        "trajectory_exclusion_allowed": False,
        "analysis_cutoff_change_allowed": False,
        "outcome_dependent_extension_allowed": False,
        "ligand_behavior_used_for_replica_selection": False,
        "required_criterion_ids": list(REQUIRED_CRITERIA),
        "failure_policy": "inconclusive_if_any_realization_fails_qc_or_stationarity",
    }
    if contract != expected_contract:
        errors.append("QC/stationarity adjudication contract differs from the frozen contract")
    results = report.get("results", [])
    if not isinstance(results, list) or [item.get("realization_id") for item in results] != list(REALIZATION_IDS):
        errors.append("QC/stationarity results must contain rep01-rep03 in order")
    else:
        for item in results:
            rid = item["realization_id"]
            raw_item = raw_by_id.get(rid)
            if raw_item is None or item.get("raw_output_run_payload_sha256") != canonical_json_sha256(raw_item):
                errors.append(f"{rid} adjudication is stale relative to its raw-output audit")
            if item.get("qc_status") not in {"pass", "fail"} or item.get("stationarity_status") not in {"pass", "fail"}:
                errors.append(f"{rid} QC/stationarity status is invalid")
            evidence = item.get("evidence")
            if not isinstance(evidence, list) or not evidence:
                errors.append(f"{rid} has no criterion-level evidence")
                continue
            criterion_ids: list[str] = []
            parsed_outcomes: list[tuple[str, dict[str, str] | None]] = []
            for index, evidence_item in enumerate(evidence):
                if not isinstance(evidence_item, dict):
                    errors.append(f"{rid} evidence {index} is not an object")
                    continue
                criterion = evidence_item.get("criterion_id")
                if not isinstance(criterion, str) or not criterion.strip():
                    errors.append(f"{rid} evidence {index} criterion_id is unresolved")
                else:
                    criterion_ids.append(criterion)
                if evidence_item.get("status") not in {"pass", "fail"}:
                    errors.append(f"{rid} evidence {index} status is invalid")
                if not isinstance(evidence_item.get("evaluator_note"), str) or not evidence_item["evaluator_note"].strip():
                    errors.append(f"{rid} evidence {index} evaluator_note is unresolved")
                if isinstance(criterion, str) and criterion in REQUIRED_CHECKS and evidence_item.get("status") in {"pass", "fail"}:
                    parsed_outcomes.append((criterion, validate_criterion_evidence(
                        package_root, evidence_item.get("artifact"), rid, criterion,
                        evidence_item["status"], canonical_json_sha256(raw_item), errors,
                    )))
            if len(criterion_ids) != len(set(criterion_ids)):
                errors.append(f"{rid} repeats a criterion_id")
            if criterion_ids != list(REQUIRED_CRITERIA):
                errors.append(f"{rid} criterion set/order must be exactly {list(REQUIRED_CRITERIA)}")
            parsed_by_id = {criterion: outcome for criterion, outcome in parsed_outcomes}
            expected_qc = "pass" if all(
                isinstance(parsed_by_id.get(criterion), dict) and parsed_by_id[criterion].get("qc_status") == "pass"
                for criterion in REQUIRED_CRITERIA
            ) else "fail"
            expected_stationarity = "pass" if all(
                isinstance(parsed_by_id.get(criterion), dict) and parsed_by_id[criterion].get("stationarity_status") == "pass"
                for criterion in REQUIRED_CRITERIA[:3]
            ) else "fail"
            if item.get("qc_status") != expected_qc:
                errors.append(f"{rid} qc_status differs from parsed criterion evidence")
            if item.get("stationarity_status") != expected_stationarity:
                errors.append(f"{rid} stationarity_status differs from parsed criterion evidence")

    all_pass = isinstance(results, list) and len(results) == 3 and all(
        item.get("qc_status") == "pass" and item.get("stationarity_status") == "pass" and
        all(evidence.get("status") == "pass" for evidence in item.get("evidence", []))
        for item in results if isinstance(item, dict)
    )
    expected_decision = "pass" if all_pass else "inconclusive"
    if report.get("overall_decision") != expected_decision or report.get("status") != expected_decision:
        errors.append("overall decision/status is inconsistent with realization and evidence results")

    approval = report.get("approval", {})
    if approval.get("approval_status") != "approved":
        errors.append("QC/stationarity report is not approved")
    for key in ("approver_name", "approver_role"):
        if not isinstance(approval.get(key), str) or not approval[key].strip():
            errors.append(f"approval.{key} is unresolved")
    parse_utc(approval.get("approved_at_utc"), "approval.approved_at_utc", errors)
    if approval.get("signature_scheme") != "sha256_canonical_json_checksum_attestation_v1":
        errors.append("approval signature scheme is invalid")
    expected_signature = report_payload_sha256(report, ("approval", "signed_payload_sha256"))
    if not allow_unsealed and approval.get("signed_payload_sha256") != expected_signature:
        errors.append("QC/stationarity signed payload checksum is invalid; report was altered or is unsealed")
    return errors


def seal_report(
    package_root: Path,
    manifest_path: Path,
    plan_path: Path,
    raw_report_path: Path,
    draft_path: Path,
    output_path: Path,
) -> None:
    if output_path.exists():
        raise ValueError(f"Refusing to overwrite sealed report: {output_path}")
    report = load_json(draft_path)
    report.setdefault("bindings", {})["validator_sha256"] = sha256(Path(__file__).resolve())
    report.setdefault("approval", {})["signed_payload_sha256"] = "UNSEALED"
    temporary = output_path.with_suffix(output_path.suffix + ".preseal.tmp")
    if temporary.exists():
        raise ValueError(f"Refusing to overwrite stale pre-seal file: {temporary}")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    try:
        errors = validate_contract(package_root, manifest_path, plan_path, raw_report_path, temporary, allow_unsealed=True)
        if errors:
            raise ValueError("Cannot seal invalid report:\n" + "\n".join(f"- {error}" for error in errors))
        report["approval"]["signed_payload_sha256"] = report_payload_sha256(
            report, ("approval", "signed_payload_sha256")
        )
        with output_path.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
    finally:
        temporary.unlink(missing_ok=True)


def write_synthetic_criterion_evidence(
    package_root: Path,
    evidence_dir: Path,
    rid: str,
    criterion: str,
    raw_payload_sha256: str,
) -> Path:
    """Build schema-valid evidence for self-tests only; never used by production CLI."""
    contract = evidence_dir / f"{criterion}.contract.json"
    if not contract.exists():
        rules = [
            {"check_id": check_id, "rule": "Frozen synthetic self-test rule."}
            for check_id in REQUIRED_CHECKS[criterion]
        ]
        contract_payload: dict[str, Any] = {
            "schema_version": "1.0", "report_type": "prespecified_qc_stationarity_criterion_contract",
            "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808", "system_id": SYSTEM_ID,
            "criterion_id": criterion, "analysis_window_ns": [200.0, 500.0],
            "frozen_before_production_review": True, "approval_status": "approved",
            "approved_at_utc": "2026-08-08T00:00:00+00:00", "check_rules": rules,
            "integrity": {"payload_sha256": "UNSEALED"},
        }
        contract_payload["integrity"]["payload_sha256"] = report_payload_sha256(
            contract_payload, ("integrity", "payload_sha256")
        )
        contract.write_text(json.dumps(contract_payload, indent=2) + "\n", encoding="utf-8")
    source = evidence_dir / f"{rid}_{criterion}.source.json"
    source.write_text('{"synthetic_source":true}\n', encoding="utf-8")
    checks = [
        {"check_id": check_id, "status": "pass", "observed": {"synthetic_pass": True},
         "rule": "Frozen synthetic self-test rule."}
        for check_id in REQUIRED_CHECKS[criterion]
    ]
    evidence: dict[str, Any] = {
        "schema_version": "1.0", "report_type": "prespecified_qc_stationarity_criterion_evidence",
        "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
        "system_id": SYSTEM_ID, "construction_id": "build01", "realization_id": rid,
        "criterion_id": criterion, "status": "pass", "analysis_window_ns": [200.0, 500.0],
        "raw_output_run_payload_sha256": raw_payload_sha256,
        "all_input_frames_retained": True, "smoothing_applied": False,
        "interpolation_applied": False, "window_changed_after_review": False,
        "criterion_contract": {"path": contract.relative_to(package_root).as_posix(), "sha256": sha256(contract)},
        "source_artifacts": [{"path": source.relative_to(package_root).as_posix(), "sha256": sha256(source)}],
        "checks": checks,
        "criterion_outcomes": {
            "qc_status": "pass",
            "stationarity_status": "not_applicable" if criterion == "spike_adjudication_completeness" else "pass",
        },
        "integrity": {"payload_sha256": "UNSEALED"},
    }
    evidence["integrity"]["payload_sha256"] = report_payload_sha256(evidence, ("integrity", "payload_sha256"))
    path = evidence_dir / f"{rid}_{criterion}.json"
    path.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    return path


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="qc_attestation_") as temporary:
        root = Path(temporary)
        config = root / "config"
        reports = root / "reports"
        evidence_dir = root / "evidence"
        config.mkdir(); reports.mkdir(); evidence_dir.mkdir()
        build = reports / "build.json"; build.write_text("{}\n", encoding="utf-8")
        manifest = {
            "study_id": "gamma_secretase_native_nirogacestat_rebuild_20260808",
            "systems": [{"id": SYSTEM_ID, "construction": {
                "id": "build01", "charmm_gui_archive": {"sha256": "a" * 64},
                "build_validation_report": {"path": "reports/build.json", "sha256": sha256(build)},
            }}],
        }
        manifest_path = config / "manifest.json"; manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        plan = {"eligibility_gate": {"all_realizations_passed_qc_and_stationarity": True, "qc_and_stationarity_report_sha256": "x"}}
        plan_path = config / "plan.json"; plan_path.write_text(json.dumps(plan), encoding="utf-8")
        raw_runs = [{"realization_id": rid, "status": "pass", "artifacts": {}} for rid in REALIZATION_IDS]
        raw = {
            "schema_version": "2.0", "status": "pass", "phase": "production", "strict": True,
            "manifest_sha256": sha256(manifest_path), "construction_archive_sha256": "a" * 64, "runs": raw_runs,
        }
        raw_path = reports / "raw.json"; raw_path.write_text(json.dumps(raw), encoding="utf-8")
        results = []
        for rid, raw_run in zip(REALIZATION_IDS, raw_runs):
            evidence_rows = []
            for criterion in REQUIRED_CRITERIA:
                evidence = write_synthetic_criterion_evidence(
                    root, evidence_dir, rid, criterion, canonical_json_sha256(raw_run)
                )
                evidence_rows.append({
                    "criterion_id": criterion, "status": "pass",
                    "artifact": {"path": evidence.relative_to(root).as_posix(), "sha256": sha256(evidence)},
                    "evaluator_note": "Synthetic criterion evidence.",
                })
            results.append({
                "realization_id": rid, "raw_output_run_payload_sha256": canonical_json_sha256(raw_run),
                "qc_status": "pass", "stationarity_status": "pass",
                "evidence": evidence_rows,
            })
        timestamp = "2026-08-08T00:00:00+00:00"
        draft = {
            "schema_version": "1.0", "report_type": "production_qc_and_stationarity_attestation",
            "study_id": manifest["study_id"], "system_id": SYSTEM_ID, "construction_id": "build01",
            "status": "pass", "strict": True, "evaluated_at_utc": timestamp,
            "evaluator": {"name": "Test", "role": "test", "software_and_versions": ["synthetic"],
                          "declaration": "All prespecified technical QC and stationarity evidence was reviewed without selection."},
            "bindings": {
                "manifest": {"path": "config/manifest.json", "sha256": sha256(manifest_path)},
                "analysis_plan": {"path": "config/plan.json", "contract_sha256": analysis_plan_contract_sha256(plan)},
                "raw_output_validation_report": {"path": "reports/raw.json", "sha256": sha256(raw_path)},
                "build_validation_report": {"path": "reports/build.json", "sha256": sha256(build)},
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
        draft_path = reports / "draft.json"; draft_path.write_text(json.dumps(draft), encoding="utf-8")
        sealed = reports / "sealed.json"
        seal_report(root, manifest_path, plan_path, raw_path, draft_path, sealed)
        if validate_contract(root, manifest_path, plan_path, raw_path, sealed):
            raise RuntimeError("valid sealed report failed")
        empty = evidence_dir / "empty.json"; empty.write_text("{}\n", encoding="utf-8")
        empty_draft = copy.deepcopy(draft)
        empty_draft["results"][0]["evidence"][0]["artifact"] = {
            "path": empty.relative_to(root).as_posix(), "sha256": sha256(empty)
        }
        empty_path = reports / "empty_draft.json"; empty_path.write_text(json.dumps(empty_draft), encoding="utf-8")
        try:
            seal_report(root, manifest_path, plan_path, raw_path, empty_path, reports / "must_not_seal.json")
        except ValueError:
            pass
        else:
            raise RuntimeError("empty arbitrary criterion evidence was sealed as PASS")
        tampered = load_json(sealed); tampered["results"][0]["stationarity_status"] = "fail"
        sealed.write_text(json.dumps(tampered), encoding="utf-8")
        if not validate_contract(root, manifest_path, plan_path, raw_path, sealed):
            raise RuntimeError("manually altered report passed")
    print("SELF-TEST PASS: sealed report binds fixed parsed criterion evidence for all three runs; empty evidence and later alteration fail closed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--analysis-plan", type=Path)
    parser.add_argument("--raw-output-report", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--seal-output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return 0
    if any(value is None for value in (args.manifest, args.analysis_plan, args.raw_output_report, args.report)):
        parser.error("--manifest, --analysis-plan, --raw-output-report, and --report are required")
    manifest_path = args.manifest.resolve(); package_root = manifest_path.parent.parent.resolve()
    plan_path = args.analysis_plan.resolve(); raw_path = args.raw_output_report.resolve(); report_path = args.report.resolve()
    if args.seal_output:
        seal_report(package_root, manifest_path, plan_path, raw_path, report_path, args.seal_output.resolve())
        print(json.dumps({"status": "sealed", "report": str(args.seal_output.resolve()), "sha256": sha256(args.seal_output.resolve())}, indent=2))
        return 0
    errors = validate_contract(package_root, manifest_path, plan_path, raw_path, report_path)
    for error in errors:
        print(f"ERROR: {error}")
    if errors:
        print(f"QC/stationarity validation failed with {len(errors)} error(s).")
        return 1
    print("QC/stationarity report validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
