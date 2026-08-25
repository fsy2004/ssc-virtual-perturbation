#!/usr/bin/env python3
"""Validate one completed equilibration branch and optionally launch the next."""

from __future__ import annotations

import argparse
import shlex
from textwrap import dedent

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
REMOTE_PYTHON = "/root/miniconda3/bin/python3"
REPLICAS = ("rep01", "rep02", "rep03")


REMOTE_VALIDATOR = dedent(
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
                fields[key.strip().lower()] = value.strip()
        return fields

    def finite_gro(path, expected_atoms=330583):
        with path.open("r", encoding="utf-8") as handle:
            handle.readline()
            atom_count = int(handle.readline().strip())
            if atom_count != expected_atoms:
                raise ValueError(f"{path}: atom count {atom_count} != {expected_atoms}")
            for atom_index in range(atom_count):
                line = handle.readline()
                if not line:
                    raise ValueError(f"{path}: truncated at atom {atom_index + 1}")
                values = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
                if not all(math.isfinite(value) for value in values):
                    raise ValueError(f"{path}: non-finite coordinate at atom {atom_index + 1}")
            box = [float(value) for value in handle.readline().split()]
            if len(box) not in (3, 9) or not all(math.isfinite(value) for value in box):
                raise ValueError(f"{path}: invalid box")
        return atom_count

    archive = base.with_suffix(".tgz")
    root_hashes = {
        "archive": sha(archive),
        "manifest": sha(base / "RELEASE_MANIFEST.json"),
        "canary": sha(base / "CANARY_VALIDATION.json"),
    }
    if root_hashes != expected_root:
        raise ValueError(f"sealed root hash mismatch: {root_hashes}")

    blocking = re.compile(
        r"LINCS WARNING|Too many LINCS warnings|constraint warning|"
        r"SETTLE.*(?:error|constraint)|(?:^|[^A-Za-z])NaN(?:[^A-Za-z]|$)|"
        r"Fatal error|Segmentation fault",
        re.I | re.M,
    )
    immutable = [
        base / "common" / "topol.top",
        base / "common" / "index.ndx",
        base / "common" / "minimized.gro",
        base / "common" / "step5_input.gro",
        base / "production_protocol.json",
    ]
    immutable.extend(base / rep / f"step6.{stage}_equilibration.mdp" for stage in range(1, 7))
    immutable_hashes = {str(path.relative_to(base)): sha(path) for path in immutable}

    stages = []
    work = base / rep / "work"
    for stage in range(1, 7):
        prefix = f"step6.{stage}_equilibration"
        mdp = base / rep / f"{prefix}.mdp"
        fields = mdp_fields(mdp)
        expected_steps = int(fields["nsteps"])
        dt_ps = float(fields["dt"])
        if expected_steps not in (125000, 250000):
            raise ValueError(f"{mdp}: unexpected nsteps={expected_steps}")
        log = work / f"{prefix}.log"
        text = log.read_text(encoding="utf-8", errors="replace")
        if "Finished mdrun on rank 0" not in text:
            raise ValueError(f"{log}: missing normal completion marker")
        if not re.search(rf"Writing checkpoint, step\s+{expected_steps}\b", text):
            raise ValueError(f"{log}: missing expected final checkpoint step {expected_steps}")
        hit = blocking.search(text)
        if hit:
            raise ValueError(f"{log}: blocking pattern {hit.group(0)!r}")
        artifacts = {}
        for extension in ("tpr", "cpt", "gro", "edr", "log"):
            path = work / f"{prefix}.{extension}"
            if not path.is_file() or path.stat().st_size <= 0:
                raise ValueError(f"missing or empty artifact: {path}")
            artifacts[extension] = {"bytes": path.stat().st_size, "sha256": sha(path)}
        atom_count = finite_gro(work / f"{prefix}.gro")
        checked = subprocess.run(
            [gmx, "check", "-e", str(work / f"{prefix}.edr")],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=120,
            check=False,
        )
        if checked.returncode != 0 or blocking.search(checked.stdout):
            raise ValueError(f"gmx check failed for {prefix}.edr")
        stages.append(
            {
                "stage": stage,
                "expected_steps": expected_steps,
                "dt_ps": dt_ps,
                "expected_time_ps": expected_steps * dt_ps,
                "atom_count": atom_count,
                "gmx_check_edr": "pass",
                "artifacts": artifacts,
            }
        )

    report = {
        "schema_version": "1.0",
        "report_type": "hash_bound_equilibration_validation",
        "release": str(base),
        "replica": rep,
        "validated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "pass",
        "sealed_root_hashes": root_hashes,
        "immutable_input_hashes": immutable_hashes,
        "checks": {
            "all_six_stages_complete": True,
            "expected_steps_and_times": True,
            "finite_coordinates_and_box": True,
            "energy_files_gmx_check": True,
            "blocking_log_scan": True,
            "required_artifacts_nonempty": True,
        },
        "stages": stages,
    }
    report_path = base / rep / "EQUILIBRATION_VALIDATION.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    report_sha = sha(report_path)
    sha_path = report_path.with_suffix(".json.sha256")
    sha_path.write_text(f"{report_sha}  {report_path.name}\n", encoding="ascii")
    print(json.dumps({"status": "pass", "replica": rep, "report": str(report_path), "sha256": report_sha}))
    '''
).strip()


def checked(client, command: str, timeout: int = 300) -> str:
    code, stdout, stderr = new_md_server.run(client, command, timeout=timeout)
    if code:
        message = (stderr or stdout)[-1200:]
        raise RuntimeError(f"remote command failed with exit code {code}: {message}")
    return stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica", required=True, choices=REPLICAS)
    parser.add_argument("--start-next", choices=REPLICAS)
    args = parser.parse_args()
    if args.start_next and REPLICAS.index(args.start_next) != REPLICAS.index(args.replica) + 1:
        raise ValueError("--start-next must be the immediately following replica")

    client = new_md_server.connect()
    try:
        validation_command = "{} -c {} {} {} {}".format(
            shlex.quote(REMOTE_PYTHON),
            shlex.quote(REMOTE_VALIDATOR),
            shlex.quote(BASE),
            shlex.quote(args.replica),
            shlex.quote(GMX),
        )
        print(checked(client, validation_command, timeout=900), end="")
        if not args.start_next:
            return
        next_rep = args.start_next
        launch = f"""
set -eu
cd {BASE}
if pgrep -x gmx >/dev/null 2>&1; then echo 'GPU worker already active' >&2; exit 4; fi
if [ -s {next_rep}_equilibrate.pid ]; then
  prior=$(cat {next_rep}_equilibrate.pid)
  if [ -d /proc/$prior ]; then echo '{next_rep} already running' >&2; exit 5; fi
fi
nohup env CUDA_VISIBLE_DEVICES=0 GMX_BIN={GMX} MDRUN_ARGS='-ntmpi 1 -ntomp 16 -pin on' \
  ./run_replica.sh {next_rep} equilibrate > {next_rep}_equilibrate.stdout \
  2> {next_rep}_equilibrate.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > {next_rep}_equilibrate.pid
echo REPLICA={next_rep}
echo RUNNER_PID=$pid
"""
        print(checked(client, launch, timeout=120), end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
