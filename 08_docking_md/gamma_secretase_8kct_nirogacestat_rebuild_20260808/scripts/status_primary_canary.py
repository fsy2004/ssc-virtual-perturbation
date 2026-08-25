#!/usr/bin/env python3
"""Report non-sensitive health and completion details for the primary canary."""

from __future__ import annotations

import importlib.util
from pathlib import Path


HELPER = Path(r"C:\Users\fsy\AppData\Local\Temp\run_primary_smoke_8687168316.py")
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v1"


def load_helper():
    spec = importlib.util.spec_from_file_location("primary_ssh", HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load primary SSH helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    module = load_helper()
    client = module.connect()
    try:
        command = f"""
set -eu
cd {BASE}
runner=$(cat canary.pid)
echo COMPLETE=$([ -e canary/COMPLETE ] && echo yes || echo no)
if [ -d /proc/$runner ]; then
  echo RUNNER_PID=$runner STATE=$(awk '/^State:/{{print $2}}' /proc/$runner/status) COMM=$(cat /proc/$runner/comm)
fi
for p in $(pgrep -x gmx || true); do
  ancestor=$(awk '/^PPid:/{{print $2}}' /proc/$p/status)
  is_child=0
  while [ "$ancestor" -gt 1 ] && [ -d /proc/$ancestor ]; do
    if [ "$ancestor" = "$runner" ]; then is_child=1; break; fi
    ancestor=$(awk '/^PPid:/{{print $2}}' /proc/$ancestor/status)
  done
  if [ "$is_child" -eq 1 ]; then
    cpus=$(for t in /proc/$p/task/*; do awk '/Cpus_allowed_list:/{{print $2}}' $t/status; done | sort -u | tr '\n' ',')
    echo GMX_PID=$p STATE=$(awk '/^State:/{{print $2}}' /proc/$p/status) CPUSETS=$cpus
  fi
done
for rep in rep01 rep02 rep03; do
  echo REPLICA=$rep
  grep -E 'Pinning threads|starting mdrun|Finished mdrun|Performance:|Time:' canary/$rep/mdrun.stderr 2>/dev/null | tail -8 || true
  if [ -s canary/$rep/canary.log ]; then
    grep -E '^Step|Writing checkpoint|Finished mdrun' canary/$rep/canary.log | tail -4 || true
    stat -c 'LOG_SIZE=%s LOG_MTIME=%y' canary/$rep/canary.log
  fi
  if [ -s canary/$rep/blocking_scan.txt ]; then
    echo BLOCKING_SCAN
    sed -n '1,20p' canary/$rep/blocking_scan.txt
  fi
done
echo RUNNER_STDERR
tail -20 canary.stderr 2>/dev/null || true
echo OUTPUT_FILES
find canary -maxdepth 2 -type f -printf '%P %s\n' | sort
df -h / /root/autodl-tmp | tail -n +2
if [ -e /sys/fs/cgroup/memory.events ]; then echo MEMORY_EVENTS; cat /sys/fs/cgroup/memory.events; fi
"""
        code, stdout, stderr = module.run(client, command, timeout=120)
        if code:
            raise RuntimeError(f"Status command failed with exit code {code}: {stderr[-500:]}")
        print(stdout, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
