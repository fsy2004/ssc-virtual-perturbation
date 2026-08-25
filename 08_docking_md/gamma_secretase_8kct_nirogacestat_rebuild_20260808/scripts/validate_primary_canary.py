#!/usr/bin/env python3
"""Independently validate the completed three-replica primary canary."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path


HELPER = Path(r"C:\Users\fsy\AppData\Local\Temp\run_primary_smoke_8687168316.py")
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v1"
GMX = "/root/GROMACS-2025.2/bin/gmx"
REPLICAS = ("rep01", "rep02", "rep03")
BLOCKING = re.compile(
    r"LINCS WARNING|Too many LINCS warnings|constraint warning|"
    r"SETTLE.*(?:error|constraint)|(?:^|[^A-Za-z])NaN(?:[^A-Za-z]|$)|"
    r"Fatal error|Segmentation fault",
    re.IGNORECASE,
)


def load_helper():
    spec = importlib.util.spec_from_file_location("primary_ssh", HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load primary SSH helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checked(module, client, command: str, timeout: int = 300):
    code, stdout, stderr = module.run(client, command, timeout=timeout)
    if code:
        raise RuntimeError(f"Remote command failed with exit code {code}: {stderr[-500:]}")
    return stdout


def read_text(sftp, path: str) -> str:
    with sftp.open(path, "r") as handle:
        return handle.read().decode("utf-8", errors="replace")


def main():
    module = load_helper()
    client = module.connect()
    try:
        required = " ".join(
            f"canary/{rep}/canary.{ext}"
            for rep in REPLICAS
            for ext in ("tpr", "gro", "edr", "log", "cpt", "xtc")
        )
        checked(module, client, f"cd {BASE} && for f in {required}; do test -s \"$f\"; done")

        for rep in REPLICAS:
            checked(
                module,
                client,
                f"""
cd {BASE}/canary/{rep}
printf 'Potential\nKinetic-En.\nTotal-Energy\nTemperature\nPressure\nBox-X\nBox-Y\nBox-Z\n0\n' |
  {GMX} energy -f canary.edr -o energy_summary.xvg > energy_summary.stdout 2> energy_summary.stderr
{GMX} check -f canary.gro > check_gro.stdout 2> check_gro.stderr
""",
            )

        sha_text = checked(
            module,
            client,
            f"cd {BASE} && sha256sum "
            + " ".join(
                f"canary/{rep}/canary.{ext}"
                for rep in REPLICAS
                for ext in ("tpr", "gro", "edr", "log", "cpt", "xtc")
            ),
        )
        hashes = {line.split()[1]: line.split()[0] for line in sha_text.splitlines()}

        report = {
            "schema": "o6u-three-replica-canary-validation-v1",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "release_root": BASE,
            "protocol": {
                "duration_ps": 5.0,
                "dt_ps": 0.001,
                "steps": 5000,
                "temperature_K": 303.15,
                "stage": "step6.1 restrained equilibration",
            },
            "false_positive_note": (
                "The initial text scanner matched the normal configuration field "
                "epsilon-rf = inf. This is not a nonfinite trajectory value; numerical "
                "finiteness is assessed from energy records and completed coordinates."
            ),
            "replicas": {},
            "pass": True,
        }

        with client.open_sftp() as sftp:
            for rep in REPLICAS:
                root = f"{BASE}/canary/{rep}"
                log = read_text(sftp, f"{root}/canary.log")
                stderr = read_text(sftp, f"{root}/mdrun.stderr")
                xvg = read_text(sftp, f"{root}/energy_summary.xvg")
                check = read_text(sftp, f"{root}/check_gro.stderr")
                matches = [m.group(0) for m in BLOCKING.finditer(log + "\n" + stderr)]
                data_lines = [
                    line for line in xvg.splitlines()
                    if line.strip() and not line.lstrip().startswith(("#", "@"))
                ]
                if not data_lines:
                    raise RuntimeError(f"No energy data rows for {rep}")
                values = [float(token) for token in data_lines[-1].split()]
                finite = all(math.isfinite(value) for value in values)
                finished = "Finished mdrun on rank 0" in log and "step 5000" in log
                # gmx check returned zero before this parse; a GRO without stored
                # velocities is acceptable because checkpoint continuity is used.
                gro_parse_pass = bool(check.strip())
                replica_pass = finite and finished and not matches and gro_parse_pass
                report["replicas"][rep] = {
                    "pass": replica_pass,
                    "finished_step_5000": finished,
                    "blocking_log_matches": matches,
                    "final_energy_row": values,
                    "final_energy_row_all_finite": finite,
                    "gmx_check_completed": gro_parse_pass,
                    "hashes": {
                        ext: hashes[f"canary/{rep}/canary.{ext}"]
                        for ext in ("tpr", "gro", "edr", "log", "cpt", "xtc")
                    },
                }
                report["pass"] = report["pass"] and replica_pass

            gro_hashes = [report["replicas"][rep]["hashes"]["gro"] for rep in REPLICAS]
            report["distinct_final_gro_hashes"] = len(set(gro_hashes)) == len(gro_hashes)
            report["pass"] = report["pass"] and report["distinct_final_gro_hashes"]
            payload = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
            report["report_payload_sha256_before_embedding"] = hashlib.sha256(payload).hexdigest()
            payload = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
            part = f"{BASE}/canary/CANARY_VALIDATION.json.part"
            final = f"{BASE}/canary/CANARY_VALIDATION.json"
            with sftp.open(part, "wb") as handle:
                handle.write(payload)
            sftp.rename(part, final)

        if not report["pass"]:
            raise RuntimeError("Canary validation failed")
        checked(
            module,
            client,
            f"cd {BASE}/canary && sha256sum CANARY_VALIDATION.json > CANARY_VALIDATION.sha256 && touch COMPLETE",
        )
        report_sha = checked(
            module,
            client,
            f"sha256sum {BASE}/canary/CANARY_VALIDATION.json",
        ).split()[0]
        print(json.dumps({"pass": True, "report_sha256": report_sha, "replicas": report["replicas"]}, indent=2))
    finally:
        client.close()


if __name__ == "__main__":
    main()
