#!/usr/bin/env python3
"""Sanitized integrity, resource, and GPU status for the retained primary."""

from __future__ import annotations

import importlib.util
from pathlib import Path


HELPER = Path(r"C:\Users\fsy\AppData\Local\Temp\run_primary_smoke_8687168316.py")
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v2"
ARCHIVE = "/root/autodl-tmp/o6u_md_release_3x500ns_v2.tgz"


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
echo INTEGRITY
sha256sum {ARCHIVE} {BASE}/RELEASE_MANIFEST.json {BASE}/CANARY_VALIDATION.json
echo EXACT_PIDFILES
for file in /root/autodl-tmp/o6u_md_release_3x500ns_v1/canary.pid; do
  [ -s "$file" ] || continue
  pid=$(cat "$file")
  if [ -d /proc/$pid ]; then
    echo FILE=$file PID=$pid STATE=$(awk '/^State:/{{print $2}}' /proc/$pid/status) COMM=$(cat /proc/$pid/comm)
  else
    echo FILE=$file PID=$pid STATE=exited
  fi
done
echo FILESYSTEMS
df -h / /root/autodl-tmp | tail -n +2
echo MEMORY_EVENTS
if [ -e /sys/fs/cgroup/memory.events ]; then cat /sys/fs/cgroup/memory.events; else echo unavailable; fi
echo GPU
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=index,name,memory.total,memory.free,utilization.gpu --format=csv,noheader
else
  echo unavailable
fi
"""
        code, stdout, stderr = module.run(client, command, timeout=120)
        if code:
            raise RuntimeError(f"Status command failed with exit code {code}: {stderr[-500:]}")
        print(stdout, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
