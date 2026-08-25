#!/usr/bin/env python3
"""Release three frozen production TPRs and optionally start one realization."""

from __future__ import annotations

import argparse
import json
import shlex
from textwrap import dedent

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
REMOTE_PYTHON = "/root/miniconda3/bin/python3"
REPLICAS = ("rep01", "rep02", "rep03")


REMOTE_RELEASE = dedent(
    r'''
    import datetime
    import hashlib
    import json
    import os
    import pathlib
    import re
    import stat
    import subprocess
    import sys

    base = pathlib.Path(sys.argv[1])
    gmx = sys.argv[2]
    replicas = ("rep01", "rep02", "rep03")
    expected_root = {
        "archive": "5a421f28afee664b5a8919db5f415f1205f35200950117bb3a67fceaba544a98",
        "manifest": "f442e1411d6f355254d5783903d96f43998a0d5758b469088d91ac8add18aee5",
        "canary": "b036732167b51064b532c454f73137921a4c0cb1ae5a0a4be83d4261d9b0d5ee",
    }

    def sha(path):
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    def mdp_fields(path):
        fields = {}
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.split(";", 1)[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                fields[key.strip().lower().replace("_", "-")] = value.strip()
        return fields

    root_hashes = {
        "archive": sha(base.with_suffix(".tgz")),
        "manifest": sha(base / "RELEASE_MANIFEST.json"),
        "canary": sha(base / "CANARY_VALIDATION.json"),
    }
    if root_hashes != expected_root:
        raise ValueError(f"sealed root hash mismatch: {root_hashes}")

    equilibration_reports = {}
    for rep in replicas:
        report_path = base / rep / "EQUILIBRATION_VALIDATION.json"
        sha_path = report_path.with_suffix(".json.sha256")
        report_sha = sha(report_path)
        recorded_sha = sha_path.read_text(encoding="ascii").split()[0]
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if recorded_sha != report_sha or report.get("status") != "pass" or report.get("replica") != rep:
            raise ValueError(f"invalid equilibration report for {rep}")
        if report.get("sealed_root_hashes") != expected_root or len(report.get("stages", [])) != 6:
            raise ValueError(f"unbound equilibration report for {rep}")
        stage6 = report["stages"][5]["artifacts"]
        for extension in ("cpt", "gro", "edr", "log", "tpr"):
            path = base / rep / "work" / f"step6.6_equilibration.{extension}"
            if sha(path) != stage6[extension]["sha256"]:
                raise ValueError(f"{rep} step6.6 {extension} changed after validation")
        equilibration_reports[rep] = report_sha

    required = {
        "integrator": "md",
        "dt": "0.004",
        "nsteps": "125000000",
        "continuation": "yes",
        "gen-vel": "no",
        "tcoupl": "v-rescale",
        "pcoupl": "c-rescale",
        "pcoupltype": "semiisotropic",
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
    }
    tprs = {}
    for rep in replicas:
        mdp = base / rep / "production_500ns.mdp"
        fields = mdp_fields(mdp)
        for key, expected in required.items():
            observed = fields.get(key, "").lower()
            if observed != expected:
                raise ValueError(f"{rep} production MDP {key}={observed!r}, expected {expected!r}")
        ref_t = fields.get("ref-t", "").split()
        if ref_t != ["303.15", "303.15", "303.15"]:
            raise ValueError(f"{rep} production ref-t is not frozen 303.15 K")
        work = base / rep / "work"
        tpr = work / "production.tpr"
        mdout = work / "production_mdout.mdp"
        log = work / "production_grompp.log"
        command = [
            gmx, "grompp", "-f", str(mdp), "-o", str(tpr),
            "-po", str(mdout), "-c", str(work / "step6.6_equilibration.gro"),
            "-t", str(work / "step6.6_equilibration.cpt"),
            "-p", str(base / "common" / "topol.top"),
            "-n", str(base / "common" / "index.ndx"), "-maxwarn", "0",
        ]
        completed = subprocess.run(
            command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=600, check=False,
        )
        log.write_text(completed.stdout, encoding="utf-8")
        if completed.returncode != 0 or not tpr.is_file() or tpr.stat().st_size <= 0:
            raise ValueError(f"grompp failed for {rep}; see {log}")
        if re.search(r"Fatal error|There was 1 warning|There were [1-9][0-9]* warnings", completed.stdout, re.I):
            raise ValueError(f"grompp warning/fatal gate failed for {rep}")
        tprs[rep] = {
            "production_mdp_sha256": sha(mdp),
            "step6_6_gro_sha256": sha(work / "step6.6_equilibration.gro"),
            "step6_6_cpt_sha256": sha(work / "step6.6_equilibration.cpt"),
            "production_tpr_sha256": sha(tpr),
            "production_tpr_bytes": tpr.stat().st_size,
            "production_mdout_sha256": sha(mdout),
            "grompp_log_sha256": sha(log),
            "grompp_maxwarn": 0,
        }

    runner = base / "run_frozen_production.sh"
    runner_text = "\n".join(
        (
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            'rep="${1:?usage: run_frozen_production.sh rep01|rep02|rep03 expected_tpr_sha256}"',
            'expected="${2:?expected production TPR SHA256 required}"',
            'case "$rep" in rep01|rep02|rep03) ;; *) echo "invalid replica: $rep" >&2; exit 2;; esac',
            'root="$(cd "$(dirname "$0")" && pwd)"',
            'work="$root/$rep/work"',
            'gmx_bin="${GMX_BIN:-gmx}"',
            'read -r -a mdrun_args <<< "${MDRUN_ARGS:-}"',
            'actual="$(sha256sum "$work/production.tpr" | awk \'{print $1}\')"',
            '[[ "$actual" == "$expected" ]] || { echo "production TPR hash mismatch" >&2; exit 3; }',
            "restart_args=()",
            'if [[ -s "$work/production.cpt" ]]; then restart_args=(-cpi "$work/production.cpt" -append); fi',
            '"$gmx_bin" mdrun "${mdrun_args[@]}" -s "$work/production.tpr" -deffnm "$work/production" "${restart_args[@]}"',
            "if grep -Ein 'LINCS WARNING|Too many LINCS warnings|constraint warning|SETTLE.*(error|constraint)|(^|[^A-Za-z])NaN([^A-Za-z]|$)|Fatal error|Segmentation fault' \"$work/production.log\" >/dev/null; then",
            '  echo "blocking dynamics warning in $work/production.log" >&2; exit 4',
            "fi",
            "",
        )
    )
    temporary_runner = runner.with_suffix(".sh.tmp")
    temporary_runner.write_text(runner_text, encoding="utf-8", newline="\n")
    temporary_runner.chmod(0o755)
    temporary_runner.replace(runner)
    runner.chmod(runner.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    report = {
        "schema_version": "1.0",
        "report_type": "three_replica_production_tpr_release",
        "status": "pass",
        "released_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "release": str(base),
        "sealed_root_hashes": root_hashes,
        "equilibration_report_sha256": equilibration_reports,
        "protocol": {
            "temperature_K": 303.15,
            "dt_ps": 0.004,
            "steps": 125000000,
            "duration_ns": 500.0,
            "continuation": True,
            "generate_velocities": False,
            "grompp_maxwarn": 0,
        },
        "runner_sha256": sha(runner),
        "replicas": tprs,
    }
    report_path = base / "PRODUCTION_TPR_RELEASE.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    report_sha = sha(report_path)
    report_path.with_suffix(".json.sha256").write_text(
        f"{report_sha}  {report_path.name}\n", encoding="ascii"
    )
    print(json.dumps({"status": "pass", "report": str(report_path), "sha256": report_sha, "replicas": tprs}))
    '''
).strip()


def checked(client, command: str, timeout: int = 300) -> str:
    code, stdout, stderr = new_md_server.run(client, command, timeout=timeout)
    if code:
        message = (stderr or stdout)[-1600:]
        raise RuntimeError(f"remote command failed with exit code {code}: {message}")
    return stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", choices=REPLICAS)
    args = parser.parse_args()
    client = new_md_server.connect()
    try:
        release_command = "{} -c {} {} {}".format(
            shlex.quote(REMOTE_PYTHON), shlex.quote(REMOTE_RELEASE),
            shlex.quote(BASE), shlex.quote(GMX),
        )
        output = checked(client, release_command, timeout=1800)
        print(output, end="")
        release = json.loads(output.strip().splitlines()[-1])
        if not args.start:
            return
        rep = args.start
        tpr_sha = release["replicas"][rep]["production_tpr_sha256"]
        launch = f"""
set -eu
cd {BASE}
if pgrep -x gmx >/dev/null 2>&1; then echo 'GPU worker already active' >&2; exit 4; fi
if [ -s {rep}_production.pid ]; then
  prior=$(cat {rep}_production.pid)
  if [ -d /proc/$prior ]; then echo '{rep} production already running' >&2; exit 5; fi
fi
nohup env CUDA_VISIBLE_DEVICES=0 GMX_BIN={GMX} MDRUN_ARGS='-ntmpi 1 -ntomp 16 -pin on' \
  ./run_frozen_production.sh {rep} {tpr_sha} > {rep}_production.stdout \
  2> {rep}_production.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > {rep}_production.pid
echo REPLICA={rep}
echo PRODUCTION_RUNNER_PID=$pid
echo TPR_SHA256={tpr_sha}
"""
        print(checked(client, launch, timeout=120), end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
