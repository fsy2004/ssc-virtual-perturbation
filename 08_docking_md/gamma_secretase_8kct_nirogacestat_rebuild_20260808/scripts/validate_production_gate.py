#!/usr/bin/env python3
"""Validate and seal a non-destructive production milestone gate."""

from __future__ import annotations

import argparse
import shlex
from textwrap import dedent

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
REMOTE_PYTHON = "/root/miniconda3/bin/python3"
REPLICAS = ("rep01", "rep02", "rep03")
GATES_PS = (100, 1000, 5000, 10000)


REMOTE_GATE = dedent(
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
    gate_ps = int(sys.argv[3])
    gmx = sys.argv[4]

    def sha(path):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    release_path = base / "PRODUCTION_TPR_RELEASE.json"
    release_sha = sha(release_path)
    recorded_release_sha = release_path.with_suffix(".json.sha256").read_text(encoding="ascii").split()[0]
    release = json.loads(release_path.read_text(encoding="utf-8"))
    if recorded_release_sha != release_sha or release.get("status") != "pass":
        raise ValueError("invalid production TPR release report")
    work = base / rep / "work"
    tpr = work / "production.tpr"
    expected_tpr_sha = release["replicas"][rep]["production_tpr_sha256"]
    if sha(tpr) != expected_tpr_sha:
        raise ValueError("production TPR hash mismatch")

    log = work / "production.log"
    edr = work / "production.edr"
    xtc = work / "production.xtc"
    for path in (log, edr, xtc):
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"missing or empty production artifact: {path}")
    log_text = log.read_text(encoding="utf-8", errors="replace")
    blocking = re.compile(
        r"LINCS WARNING|Too many LINCS warnings|constraint warning|"
        r"SETTLE.*(?:error|constraint)|(?:^|[^A-Za-z])NaN(?:[^A-Za-z]|$)|"
        r"Fatal error|Segmentation fault",
        re.I | re.M,
    )
    hit = blocking.search(log_text)
    if hit:
        raise ValueError(f"blocking production log pattern: {hit.group(0)!r}")
    progress = re.findall(r"^\s*(\d+)\s+([0-9]+(?:\.[0-9]+)?)\s*$", log_text, re.M)
    if not progress:
        raise ValueError("no production Step/Time progress rows")
    latest_step, latest_time_ps = int(progress[-1][0]), float(progress[-1][1])
    if latest_time_ps < gate_ps:
        raise ValueError(f"gate {gate_ps} ps not reached; latest {latest_time_ps} ps")

    tool_checks = {}
    for label, command in {
        "edr": [gmx, "check", "-e", str(edr)],
        "xtc": [gmx, "check", "-f", str(xtc)],
    }.items():
        checked = subprocess.run(
            command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=300, check=False,
        )
        if checked.returncode != 0 or blocking.search(checked.stdout):
            raise ValueError(f"gmx check failed for production {label}")
        tool_checks[label] = "pass"

    energy_xvg = work / f"production_gate_{gate_ps}ps_energy.xvg"
    selection = "Potential\nTotal-Energy\nTemperature\nPressure\nVolume\n0\n"
    energy = subprocess.run(
        [gmx, "energy", "-f", str(edr), "-o", str(energy_xvg)],
        input=selection, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=300, check=False,
    )
    if energy.returncode != 0 or not energy_xvg.is_file():
        raise ValueError("gmx energy extraction failed")
    rows = []
    for raw in energy_xvg.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw or raw[0] in "#@":
            continue
        values = [float(value) for value in raw.split()]
        if not values or not all(math.isfinite(value) for value in values):
            raise ValueError("non-finite value in production energy series")
        rows.append(values)
    if not rows or rows[-1][0] < gate_ps:
        raise ValueError("finite energy series does not reach gate")

    artifacts = {}
    for path in (tpr, log, edr, xtc, energy_xvg):
        artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    cpt = work / "production.cpt"
    if cpt.is_file() and cpt.stat().st_size > 0:
        artifacts[cpt.name] = {"bytes": cpt.stat().st_size, "sha256": sha(cpt)}

    report = {
        "schema_version": "1.0",
        "report_type": "production_milestone_gate",
        "status": "pass",
        "replica": rep,
        "gate_ps": gate_ps,
        "validated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "latest_step": latest_step,
        "latest_time_ps": latest_time_ps,
        "production_release_sha256": release_sha,
        "production_tpr_sha256": expected_tpr_sha,
        "checks": {
            "blocking_log_scan": "pass",
            "gmx_check": tool_checks,
            "finite_energy_series": "pass",
            "trajectory_and_energy_reach_gate": "pass",
            "checkpoint_present": cpt.is_file() and cpt.stat().st_size > 0,
        },
        "artifacts": artifacts,
    }
    report_path = base / rep / f"PRODUCTION_GATE_{gate_ps}ps.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    report_sha = sha(report_path)
    report_path.with_suffix(".json.sha256").write_text(
        f"{report_sha}  {report_path.name}\n", encoding="ascii"
    )
    print(json.dumps({"status": "pass", "replica": rep, "gate_ps": gate_ps,
                      "latest_time_ps": latest_time_ps, "report": str(report_path),
                      "sha256": report_sha}))
    '''
).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica", required=True, choices=REPLICAS)
    parser.add_argument("--gate-ps", required=True, type=int, choices=GATES_PS)
    args = parser.parse_args()
    compile(REMOTE_GATE, "remote_production_gate", "exec")
    client = new_md_server.connect()
    try:
        command = "{} -c {} {} {} {} {}".format(
            shlex.quote(REMOTE_PYTHON), shlex.quote(REMOTE_GATE), shlex.quote(BASE),
            shlex.quote(args.replica), args.gate_ps, shlex.quote(GMX),
        )
        code, stdout, stderr = new_md_server.run(client, command, timeout=900)
        if code:
            raise RuntimeError(f"production gate failed with exit code {code}: {(stderr or stdout)[-1600:]}")
        print(stdout, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
