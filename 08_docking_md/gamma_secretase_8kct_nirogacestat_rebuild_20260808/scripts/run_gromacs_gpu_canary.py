#!/usr/bin/env python3
"""Run a tiny, non-scientific GROMACS CUDA offload canary.

This verifies the server/runtime before any project equilibration or production.
It deliberately uses a disposable water box and must never be cited as study data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_text(path: Path, value: str) -> None:
    path.write_text(value, encoding="utf-8", newline="\n")


def run(argv: list[str], cwd: Path, stem: str, env: dict[str, str] | None = None) -> dict:
    completed = subprocess.run(argv, cwd=cwd, text=True, capture_output=True, env=env)
    stdout_path = cwd / f"{stem}.stdout.log"
    stderr_path = cwd / f"{stem}.stderr.log"
    write_text(stdout_path, completed.stdout)
    write_text(stderr_path, completed.stderr)
    record = {
        "argv": argv,
        "cwd": str(cwd.resolve()),
        "returncode": completed.returncode,
        "stdout": stdout_path.name,
        "stderr": stderr_path.name,
        "stdout_sha256": sha256(stdout_path),
        "stderr_sha256": sha256(stderr_path),
    }
    if completed.returncode != 0:
        raise RuntimeError(json.dumps(record, indent=2))
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmx", required=True, help="Absolute GROMACS executable")
    parser.add_argument("--spc216", required=True, help="GROMACS spc216.gro path")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    gmx = Path(args.gmx).resolve()
    spc216 = Path(args.spc216).resolve()
    out = Path(args.output_dir).resolve()
    if not gmx.is_file() or not os.access(gmx, os.X_OK):
        raise SystemExit(f"GROMACS executable unavailable: {gmx}")
    if not spc216.is_file():
        raise SystemExit(f"Water coordinate source unavailable: {spc216}")
    if args.threads < 1:
        raise SystemExit("--threads must be positive")

    if out.exists():
        if any(out.iterdir()):
            raise SystemExit(f"Refusing to overwrite non-empty output directory: {out}")
    else:
        out.mkdir(parents=True)

    topology = """; Disposable environment canary only -- not study data.
#include "amber99sb-ildn.ff/forcefield.itp"
#include "amber99sb-ildn.ff/tip3p.itp"

[ system ]
Disposable CUDA water canary

[ molecules ]
"""
    em_mdp = """; Disposable environment canary only -- not study data.
integrator               = steep
nsteps                   = 1000
emtol                    = 500.0
emstep                   = 0.01
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
rvdw                     = 1.0
pbc                      = xyz
constraints              = h-bonds
constraint-algorithm     = lincs
"""
    nvt_mdp = """; Disposable environment canary only -- not study data.
integrator               = md
dt                       = 0.002
nsteps                   = 2500
continuation             = no
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
rvdw                     = 1.0
DispCorr                 = EnerPres
pbc                      = xyz
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = no
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = 20260809
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstenergy                = 100
nstlog                   = 100
nstxout-compressed       = 500
compressed-x-precision   = 1000
"""
    write_text(out / "topol.top", topology)
    write_text(out / "em.mdp", em_mdp)
    write_text(out / "nvt.mdp", nvt_mdp)

    runtime_env = dict(os.environ)
    runtime_env["OMP_NUM_THREADS"] = str(args.threads)
    commands: list[dict] = []
    commands.append(run([
        str(gmx), "solvate", "-cs", str(spc216), "-box", "3", "3", "3",
        "-o", "water.gro", "-p", "topol.top",
    ], out, "01_solvate", runtime_env))
    commands.append(run([
        str(gmx), "grompp", "-f", "em.mdp", "-c", "water.gro", "-p", "topol.top",
        "-o", "em.tpr", "-po", "em_mdout.mdp", "-maxwarn", "0",
    ], out, "02_em_grompp", runtime_env))
    commands.append(run([
        str(gmx), "mdrun", "-deffnm", "em", "-ntmpi", "1", "-ntomp", str(args.threads),
        # GROMACS does not support PME GPU offload for the steep integrator.
        # The subsequent dynamical NVT stage is the CUDA/offload test.
        "-nb", "cpu", "-pme", "cpu", "-bonded", "cpu", "-pin", "on",
    ], out, "03_em_mdrun", runtime_env))
    commands.append(run([
        str(gmx), "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-p", "topol.top",
        "-o", "nvt.tpr", "-po", "nvt_mdout.mdp", "-maxwarn", "0",
    ], out, "04_nvt_grompp", runtime_env))
    commands.append(run([
        str(gmx), "mdrun", "-deffnm", "nvt", "-ntmpi", "1", "-ntomp", str(args.threads),
        # Rigid water has no implemented bonded term to offload. Project-system
        # bonded offload is validated later by the project-specific canary.
        "-nb", "gpu", "-pme", "gpu", "-bonded", "cpu", "-update", "gpu",
        "-gpu_id", "0", "-pin", "on",
    ], out, "05_nvt_mdrun", runtime_env))

    combined = "\n".join(
        p.read_text(encoding="utf-8", errors="replace")
        for p in [
            out / "03_em_mdrun.stdout.log", out / "03_em_mdrun.stderr.log", out / "em.log",
            out / "05_nvt_mdrun.stdout.log", out / "05_nvt_mdrun.stderr.log", out / "nvt.log",
        ] if p.exists()
    )
    lower = combined.lower()
    forbidden_patterns = {
        "fatal error": r"fatal error",
        "lincs warning": r"lincs warning",
        "segmentation fault": r"segmentation fault",
        "nan": r"\bnan\b",
    }
    forbidden = [label for label, pattern in forbidden_patterns.items() if re.search(pattern, lower)]
    gpu_evidence = {
        "mentions_cuda": "cuda" in lower,
        "mentions_gpu": "gpu" in lower,
        "mentions_rtx_5090": "5090" in lower,
        "nonbonded_gpu": ("nonbonded" in lower and "gpu" in lower),
        "pme_gpu": ("pme" in lower and "gpu" in lower),
        "update_gpu": ("update" in lower and "gpu" in lower),
        "finished": ("finished mdrun" in lower or "performance:" in lower),
    }
    required = ["mentions_cuda", "mentions_gpu", "nonbonded_gpu", "pme_gpu", "finished"]
    status = "pass" if not forbidden and all(gpu_evidence[key] for key in required) else "fail"

    artifacts = {}
    for path in sorted(p for p in out.iterdir() if p.is_file() and p.name != "GROMACS_GPU_CANARY_REPORT.json"):
        artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    report = {
        "schema_version": 1,
        "status": status,
        "scope": "non-scientific disposable GROMACS CUDA environment canary",
        "explicit_exclusions": ["study system", "parameter validation", "equilibration", "production", "scientific evidence"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "python": sys.version,
        "gmx": str(gmx),
        "spc216": str(spc216),
        "threads": args.threads,
        "runtime_environment_overrides": {"OMP_NUM_THREADS": str(args.threads)},
        "commands": commands,
        "gpu_evidence": gpu_evidence,
        "forbidden_log_tokens_found": forbidden,
        "artifacts": artifacts,
    }
    report_path = out / "GROMACS_GPU_CANARY_REPORT.json"
    write_text(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": status, "report": str(report_path), "report_sha256": sha256(report_path)}, indent=2))
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
