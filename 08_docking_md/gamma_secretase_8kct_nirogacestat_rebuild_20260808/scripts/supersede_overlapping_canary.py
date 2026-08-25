#!/usr/bin/env python3
"""Gracefully supersede the audit-only canary with overlapping CPU pinning.

Credentials are read only by the existing local SSH helper and are never
printed.  This script validates exact PIDs/comm names, sends SIGINT only to
the expected GROMACS processes, waits for clean exit, and archives the old
canary artifacts without deleting them.
"""

from __future__ import annotations

import importlib.util
import time
from pathlib import Path


HELPER = Path(r"C:\Users\fsy\AppData\Local\Temp\run_primary_smoke_8687168316.py")
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v1"
GMX_PIDS = (496662, 496663, 496664)
BASH_PIDS = (496649, 496652, 496659, 496660, 496661)
AUDIT = f"{BASE}/audit_canary_superseded_overlapping_pinning_20260816"


def load_helper():
    spec = importlib.util.spec_from_file_location("primary_ssh", HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load primary SSH helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(module, client, command: str, timeout: int = 120):
    code, stdout, stderr = module.run(client, command, timeout=timeout)
    if code:
        raise RuntimeError(f"Remote command failed with exit code {code}: {stderr[-500:]}")
    return stdout


def describe(module, client):
    pids = " ".join(map(str, GMX_PIDS + BASH_PIDS))
    command = f"""
for p in {pids}; do
  if [ -d /proc/$p ]; then
    state=$(awk '/^State:/{{print $2}}' /proc/$p/status)
    comm=$(cat /proc/$p/comm)
    ppid=$(awk '/^PPid:/{{print $2}}' /proc/$p/status)
    echo PID=$p PPID=$ppid STATE=$state COMM=$comm
  fi
done
"""
    return run(module, client, command)


def main():
    module = load_helper()
    client = module.connect()
    try:
        print(describe(module, client), end="")

        for pid in GMX_PIDS:
            command = f"""
if [ -d /proc/{pid} ]; then
  test \"$(cat /proc/{pid}/comm)\" = gmx
  kill -INT {pid}
fi
"""
            run(module, client, command)

        deadline = time.time() + 45
        while time.time() < deadline:
            alive = run(
                module,
                client,
                "for p in " + " ".join(map(str, GMX_PIDS)) + "; do [ -d /proc/$p ] && echo $p; done; true",
            ).split()
            if not alive:
                break
            time.sleep(3)
        else:
            raise RuntimeError(f"GROMACS processes did not exit after SIGINT: {alive}")

        for pid in BASH_PIDS:
            command = f"""
if [ -d /proc/{pid} ]; then
  test \"$(cat /proc/{pid}/comm)\" = bash
  kill -TERM {pid}
fi
"""
            run(module, client, command)

        command = f"""
set -eu
test ! -e {AUDIT}
mkdir -p {AUDIT}
for name in canary canary.stdout canary.stderr canary.pid; do
  if [ -e {BASE}/$name ]; then mv {BASE}/$name {AUDIT}/; fi
done
if [ -e /root/canary.pid ]; then mv /root/canary.pid {AUDIT}/root_canary.pid; fi
printf '%s\n' 'superseded_reason=all_three_replicas_were_pinned_to_cpus_0-7' > {AUDIT}/SUPERSEDED.txt
printf '%s\n' 'scientific_status=audit_only_technical_supersession' >> {AUDIT}/SUPERSEDED.txt
find {AUDIT} -maxdepth 2 -type f -printf '%P\n' | sort
"""
        print(run(module, client, command), end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
