#!/usr/bin/env python3
"""Seal a finished 500 ns production run and start the next frozen replica."""

from __future__ import annotations

import argparse
import json
import shlex
from textwrap import dedent

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
REMOTE_PYTHON = "/root/miniconda3/bin/python3"
RELEASE_SHA = "a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa"
REPLICAS = ("rep01", "rep02", "rep03")


REMOTE_VALIDATE = dedent(
    r'''
    import datetime
    import hashlib
    import json
    import math
    import pathlib
    import re
    import subprocess
    import sys

    base = pathlib.Path(sys.argv[1])
    rep = sys.argv[2]
    gmx = sys.argv[3]
    expected_release_sha = sys.argv[4]

    def sha(path):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    blocking = re.compile(
        r"LINCS WARNING|Too many LINCS warnings|constraint warning|"
        r"SETTLE.*(?:error|constraint)|(?:^|[^A-Za-z])NaN(?:[^A-Za-z]|$)|"
        r"Fatal error|Segmentation fault", re.I | re.M,
    )
    release_path = base / "PRODUCTION_TPR_RELEASE.json"
    release_sha = sha(release_path)
    recorded_release_sha = release_path.with_suffix(".json.sha256").read_text(
        encoding="ascii"
    ).split()[0]
    release = json.loads(release_path.read_text(encoding="utf-8"))
    if release_sha != expected_release_sha or recorded_release_sha != release_sha:
        raise ValueError("production release hash mismatch")
    if release.get("status") != "pass":
        raise ValueError("production release is not passing")

    work = base / rep / "work"
    paths = {ext: work / f"production.{ext}" for ext in
             ("tpr", "log", "edr", "xtc", "cpt", "gro")}
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"missing or empty final artifact: {label}")

    expected_tpr_sha = release["replicas"][rep]["production_tpr_sha256"]
    if sha(paths["tpr"]) != expected_tpr_sha:
        raise ValueError("production TPR hash mismatch")

    log_text = paths["log"].read_text(encoding="utf-8", errors="replace")
    hit = blocking.search(log_text)
    if hit:
        raise ValueError(f"blocking final log pattern: {hit.group(0)!r}")
    progress = re.findall(r"^\s*(\d+)\s+([0-9]+(?:\.[0-9]+)?)\s*$", log_text, re.M)
    if not progress:
        raise ValueError("no production Step/Time rows")
    final_step, final_time_ps = int(progress[-1][0]), float(progress[-1][1])
    if final_step != 125000000 or abs(final_time_ps - 500000.0) > 1e-6:
        raise ValueError(f"unexpected final Step/Time: {final_step}, {final_time_ps}")
    if "Finished mdrun" not in log_text:
        raise ValueError("final log lacks Finished mdrun marker")

    checks = {}
    for label, command in {
        "edr": [gmx, "check", "-e", str(paths["edr"])],
        "xtc": [gmx, "check", "-f", str(paths["xtc"])],
        "gro": [gmx, "check", "-f", str(paths["gro"])],
    }.items():
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=900, check=False)
        if result.returncode != 0 or blocking.search(result.stdout):
            raise ValueError(f"gmx check failed for {label}")
        checks[label] = "pass"

    dump = subprocess.run([gmx, "dump", "-cp", str(paths["cpt"])], text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          timeout=300, check=False)
    if dump.returncode != 0 or blocking.search(dump.stdout):
        raise ValueError("checkpoint dump failed")
    cpt_step_match = re.search(r"\bstep\s*=\s*(\d+)", dump.stdout)
    cpt_time_match = re.search(r"(?:^|\n)\s*t\s*=\s*([0-9.eE+-]+)", dump.stdout)
    if not cpt_step_match or int(cpt_step_match.group(1)) != 125000000:
        raise ValueError("checkpoint does not bind final step")
    if not cpt_time_match or abs(float(cpt_time_match.group(1)) - 500000.0) > 1e-6:
        raise ValueError("checkpoint does not bind final time")
    checks["cpt"] = "pass"

    energy_xvg = work / "production_completion_energy.xvg"
    energy = subprocess.run(
        [gmx, "energy", "-f", str(paths["edr"]), "-o", str(energy_xvg)],
        input="Potential\nTotal-Energy\nTemperature\nPressure\nVolume\nDensity\n0\n",
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=900, check=False,
    )
    if energy.returncode != 0 or not energy_xvg.is_file():
        raise ValueError("final energy extraction failed")
    energy_rows = []
    for raw in energy_xvg.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw or raw[0] in "#@":
            continue
        values = [float(value) for value in raw.split()]
        if not values or not all(math.isfinite(value) for value in values):
            raise ValueError("non-finite value in final energy series")
        energy_rows.append(values)
    if not energy_rows or abs(energy_rows[-1][0] - 500000.0) > 1e-6:
        raise ValueError("final energy series does not reach exactly 500 ns")

    recovery_audits = sorted((base / "audit").glob(f"{rep}_cuda700_*/**/*"))
    recovery_files = [p for p in recovery_audits if p.is_file()]
    recovery = {
        "audit_files": len(recovery_files),
        "audit_sha256": {str(p.relative_to(base)): sha(p) for p in recovery_files},
    }

    artifacts = {}
    for path in (*paths.values(), energy_xvg):
        artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    report = {
        "schema_version": "1.0",
        "report_type": "production_500ns_completion",
        "status": "pass",
        "replica": rep,
        "validated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "final_step": final_step,
        "final_time_ps": final_time_ps,
        "production_release_sha256": release_sha,
        "production_tpr_sha256": expected_tpr_sha,
        "checks": {
            "blocking_log_scan": "pass",
            "finished_mdrun_marker": "pass",
            "exact_final_step_time": "pass",
            "gmx_readability": checks,
            "finite_energy_series": "pass",
            "checkpoint_final_step_time": "pass",
            "cuda_recovery_audit_bound": bool(recovery_files),
        },
        "recovery_audit": recovery,
        "artifacts": artifacts,
    }
    report_path = base / rep / "PRODUCTION_COMPLETION_500NS.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(report_path)
    report_sha = sha(report_path)
    report_path.with_suffix(".json.sha256").write_text(
        f"{report_sha}  {report_path.name}\n", encoding="ascii"
    )
    print(json.dumps({"status": "pass", "replica": rep,
                      "final_time_ps": final_time_ps,
                      "report": str(report_path), "sha256": report_sha}))
    '''
).strip()


REMOTE_START = dedent(
    r'''
    import hashlib
    import json
    import os
    import pathlib
    import subprocess
    import sys

    base = pathlib.Path(sys.argv[1])
    rep = sys.argv[2]
    gmx = sys.argv[3]
    expected_release_sha = sys.argv[4]

    def sha(path):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    release_path = base / "PRODUCTION_TPR_RELEASE.json"
    if sha(release_path) != expected_release_sha:
        raise ValueError("production release hash mismatch before launch")
    release = json.loads(release_path.read_text(encoding="utf-8"))
    work = base / rep / "work"
    tpr_sha = release["replicas"][rep]["production_tpr_sha256"]
    if sha(work / "production.tpr") != tpr_sha:
        raise ValueError("next production TPR hash mismatch")
    for name in ("production.cpt", "production.log", "production.edr",
                 "production.xtc", "production.gro"):
        if (work / name).exists():
            raise ValueError(f"refusing non-pristine first start: {name} exists")
    if subprocess.run(["pgrep", "-x", "gmx"], stdout=subprocess.DEVNULL,
                      stderr=subprocess.DEVNULL).returncode == 0:
        raise ValueError("GPU worker already active")

    stdout = (base / f"{rep}_production.stdout").open("ab")
    stderr = (base / f"{rep}_production.stderr").open("ab")
    env = os.environ.copy()
    env.update({"CUDA_VISIBLE_DEVICES": "0", "GMX_BIN": gmx,
                "MDRUN_ARGS": "-ntmpi 1 -ntomp 16 -pin on"})
    process = subprocess.Popen(
        [str(base / "run_frozen_production.sh"), rep, tpr_sha], cwd=base,
        stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
        env=env, start_new_session=True,
    )
    (base / f"{rep}_production.pid").write_text(f"{process.pid}\n", encoding="ascii")
    print(json.dumps({"status": "started", "replica": rep,
                      "runner_pid": process.pid, "tpr_sha256": tpr_sha,
                      "append": False}))
    '''
).strip()


def run_remote(client, script: str, args: list[str], timeout: int) -> dict:
    command = "{} -c {} {}".format(
        shlex.quote(REMOTE_PYTHON), shlex.quote(script),
        " ".join(shlex.quote(arg) for arg in args),
    )
    code, stdout, stderr = new_md_server.run(client, command, timeout=timeout)
    if code:
        raise RuntimeError((stderr or stdout)[-2000:])
    result = json.loads(stdout.strip().splitlines()[-1])
    print(json.dumps(result, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--completed", required=True, choices=REPLICAS)
    parser.add_argument("--start-next", choices=REPLICAS)
    args = parser.parse_args()
    compile(REMOTE_VALIDATE, "remote_validate_completion", "exec")
    compile(REMOTE_START, "remote_start_next", "exec")
    client = new_md_server.connect()
    try:
        run_remote(client, REMOTE_VALIDATE,
                   [BASE, args.completed, GMX, RELEASE_SHA], timeout=3600)
        if args.start_next:
            run_remote(client, REMOTE_START,
                       [BASE, args.start_next, GMX, RELEASE_SHA], timeout=300)
    finally:
        client.close()


if __name__ == "__main__":
    main()
