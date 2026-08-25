#!/usr/bin/env python3
"""Run a frozen O6U MP2 Hessian QCSchema plan with bounded parallelism."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from run_o6u_mp2_optimization_canary import sha256


def atomic_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def validate_result(task: dict) -> dict:
    inp = Path(task["input"])
    out = Path(task["output"])
    if not out.is_file() or out.stat().st_size == 0:
        raise ValueError("missing_or_empty_output")
    source = json.loads(inp.read_text(encoding="utf-8"))
    result = json.loads(out.read_text(encoding="utf-8"))
    if result.get("success") is not True or result.get("schema_name") != "qcschema_output":
        raise ValueError("not_successful_qcschema_output")
    if result.get("driver") != "gradient" or source.get("driver") != "gradient":
        raise ValueError("driver_mismatch")
    if result.get("model") != source.get("model"):
        raise ValueError("model_mismatch")
    if result.get("molecule", {}).get("symbols") != source.get("molecule", {}).get("symbols"):
        raise ValueError("atom_identity_mismatch")
    grad = result.get("return_result")
    natom = len(source["molecule"]["symbols"])
    if not isinstance(grad, list) or len(grad) != 3 * natom:
        raise ValueError("gradient_shape_mismatch")
    return {
        "output_sha256": sha256(out),
        "output_bytes": out.stat().st_size,
        "success": True,
    }


def run_one(task: dict, psi4: str, threads: int, memory_gib: int) -> dict:
    inp = Path(task["input"])
    out = Path(task["output"])
    scratch = Path(task["scratch"])
    log = out.with_suffix(out.suffix + ".launcher.log")
    if out.exists() or log.exists():
        try:
            valid = validate_result(task)
        except Exception as exc:
            return {"index": task["index"], "label": task["label"], "status": "refused_existing_invalid", "error": str(exc)}
        return {"index": task["index"], "label": task["label"], "status": "skipped_existing_valid", **valid}
    scratch.mkdir(parents=True, exist_ok=False)
    command = [psi4, "--qcschema", "-i", str(inp), "-o", str(out), "-n", str(threads), "--memory", f"{memory_gib}GiB", "-s", str(scratch)]
    started = time.monotonic()
    proc = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    elapsed = time.monotonic() - started
    log.write_text(proc.stdout or "", encoding="utf-8", newline="\n")
    base = {"index": task["index"], "label": task["label"], "exit_code": proc.returncode, "elapsed_seconds": elapsed, "launcher_log": str(log), "launcher_log_sha256": sha256(log)}
    if proc.returncode != 0:
        return {**base, "status": "failed_exit"}
    try:
        valid = validate_result(task)
    except Exception as exc:
        return {**base, "status": "failed_validation", "error": str(exc)}
    scratch_bytes = sum(p.stat().st_size for p in scratch.rglob("*") if p.is_file())
    shutil.rmtree(scratch)
    return {**base, "status": "pass", "scratch_cleanup": "removed_after_validated_result", "scratch_bytes_removed": scratch_bytes, **valid}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--manifest-sha256", required=True)
    ap.add_argument("--psi4", required=True, type=Path)
    ap.add_argument("--max-parallel", type=int, default=8)
    ap.add_argument("--threads-per-job", type=int)
    ap.add_argument("--state-file", type=Path)
    args = ap.parse_args()
    manifest_path = args.manifest.resolve()
    if sha256(manifest_path) != args.manifest_sha256:
        raise SystemExit("Manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "pass_parallel_qcschema_plan_generated_no_qm_executed":
        raise SystemExit("Manifest is not a released no-QM plan")
    manifest_threads = int(manifest["threads_per_job"])
    threads = args.threads_per_job if args.threads_per_job is not None else manifest_threads
    memory_gib = int(manifest["memory_gib_per_job"])
    if threads not in {2, manifest_threads}:
        raise SystemExit("Threads per job must be the frozen manifest value or audited 2-thread scaling value")
    if not 1 <= args.max_parallel <= 12 or args.max_parallel * threads > 24 or args.max_parallel * memory_gib > 96:
        raise SystemExit("Parallel resources exceed effective 24-CPU/96-GiB cgroup-aware bounds")
    psi4 = args.psi4.resolve()
    if not psi4.is_file():
        raise SystemExit("Psi4 executable missing")
    tasks = manifest["tasks"]
    if len(tasks) != 445 or sorted(t["index"] for t in tasks) != list(range(445)):
        raise SystemExit("Unexpected task count or indices")
    for task in tasks:
        inp = Path(task["input"])
        if sha256(inp) != task["input_sha256"]:
            raise SystemExit(f"Input hash mismatch: {inp}")

    state_path = args.state_file.resolve() if args.state_file else manifest_path.parent / "O6U_MP2_HESSIAN_PARALLEL_STATE.json"
    if state_path.exists():
        raise SystemExit(f"Refusing to overwrite existing state file: {state_path}")
    state = {
        "schema_version": "2.0",
        "report_type": "o6u_mp2_hessian_parallel_state",
        "status": "running",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": {"path": str(manifest_path), "sha256": args.manifest_sha256},
        "max_parallel": args.max_parallel,
        "manifest_threads_per_job": manifest_threads,
        "threads_per_job": threads,
        "memory_gib_per_job": memory_gib,
        "results": [],
    }
    atomic_json(state_path, state)
    failed = False
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_parallel) as pool:
        future_map = {pool.submit(run_one, task, str(psi4), threads, memory_gib): task for task in tasks}
        for future in concurrent.futures.as_completed(future_map):
            try:
                result = future.result()
            except Exception as exc:
                task = future_map[future]
                result = {"index": task["index"], "label": task["label"], "status": "controller_exception", "error": repr(exc)}
            state["results"].append(result)
            state["results"].sort(key=lambda x: x["index"])
            if result["status"] not in {"pass", "skipped_existing_valid"}:
                failed = True
            state["completed_count"] = len(state["results"])
            state["pass_count"] = sum(r["status"] in {"pass", "skipped_existing_valid"} for r in state["results"])
            state["failed_count"] = state["completed_count"] - state["pass_count"]
            state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
            atomic_json(state_path, state)
    state["status"] = "fail_closed_worker_errors" if failed else "pass_all_qcschema_gradients_validated"
    state["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    atomic_json(state_path, state)
    print(json.dumps({"status": state["status"], "state": str(state_path), "sha256": sha256(state_path), "pass_count": state["pass_count"]}, sort_keys=True))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
