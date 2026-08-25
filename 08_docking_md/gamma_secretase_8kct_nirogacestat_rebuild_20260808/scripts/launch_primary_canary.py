#!/usr/bin/env python3
"""Upload and launch the three-replica canary on the retained primary node."""

from __future__ import annotations

import importlib.util
import time
from pathlib import Path


HELPER = Path(r"C:\Users\fsy\AppData\Local\Temp\run_primary_smoke_8687168316.py")
RUNNER = Path(__file__).with_name("run_three_replica_canary.sh")
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v1"
GMX = "/root/GROMACS-2025.2/bin/gmx"


def load_helper():
    spec = importlib.util.spec_from_file_location("primary_ssh", HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load primary SSH helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checked(module, client, command: str, timeout: int = 120):
    code, stdout, stderr = module.run(client, command, timeout=timeout)
    if code:
        raise RuntimeError(f"Remote command failed with exit code {code}: {stderr[-500:]}")
    return stdout


def main():
    module = load_helper()
    client = module.connect()
    try:
        checked(module, client, f"test -d {BASE} && mkdir -p {BASE}/scripts")
        with client.open_sftp() as sftp:
            sftp.put(str(RUNNER), f"{BASE}/scripts/run_three_replica_canary.sh")
        checked(
            module,
            client,
            f"""
set -eu
cd {BASE}
test ! -e canary
chmod 0755 scripts/run_three_replica_canary.sh
nohup bash scripts/run_three_replica_canary.sh {BASE} {GMX} > canary.stdout 2> canary.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > canary.pid
echo RUNNER_PID=$pid
""",
        )
        time.sleep(5)
        summary = checked(
            module,
            client,
            f"""
set -eu
cd {BASE}
runner=$(cat canary.pid)
if [ -d /proc/$runner ]; then
  echo RUNNER_PID=$runner STATE=$(awk '/^State:/{{print $2}}' /proc/$runner/status) COMM=$(cat /proc/$runner/comm)
fi
for p in $(pgrep -x gmx || true); do
  ppid=$(awk '/^PPid:/{{print $2}}' /proc/$p/status)
  ancestor=$ppid
  is_child=0
  while [ "$ancestor" -gt 1 ] && [ -d /proc/$ancestor ]; do
    if [ "$ancestor" = "$runner" ]; then is_child=1; break; fi
    ancestor=$(awk '/^PPid:/{{print $2}}' /proc/$ancestor/status)
  done
  if [ "$is_child" -eq 1 ]; then
    cpus=$(for t in /proc/$p/task/*; do awk '/Cpus_allowed_list:/{{print $2}}' $t/status; done | sort -u | tr '\n' ',')
    echo GMX_PID=$p PPID=$ppid STATE=$(awk '/^State:/{{print $2}}' /proc/$p/status) CPUSETS=$cpus
  fi
done
""",
        )
        print(summary, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
