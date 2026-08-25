#!/usr/bin/env python3
"""Read-only sanitized readiness probe for the current MD server."""

from __future__ import annotations

import new_md_server


def main() -> None:
    client = new_md_server.connect()
    try:
        command = r"""
set -eu
echo HOSTNAME=$(hostname)
echo KERNEL=$(uname -srm)
echo CPU_COUNT=$(getconf _NPROCESSORS_ONLN)
echo MEMORY
awk '/MemTotal:|MemAvailable:/{print $1,$2,$3}' /proc/meminfo
echo FILESYSTEMS
df -h / /root/autodl-tmp 2>/dev/null | tail -n +2 || df -h /
echo GPU
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=index,name,memory.total,memory.free,utilization.gpu,compute_mode --format=csv,noheader
else
  echo unavailable
fi
echo GROMACS
for candidate in /root/GROMACS-2025.2/bin/gmx /usr/local/gromacs/bin/gmx /usr/bin/gmx; do
  if [ -x "$candidate" ]; then "$candidate" --version | sed -n '1,8p'; break; fi
done
echo ACTIVE_MDRUN
for p in $(pgrep -x gmx 2>/dev/null || true); do
  echo PID=$p STATE=$(awk '/^State:/{print $2}' /proc/$p/status) COMM=$(cat /proc/$p/comm)
done
echo MEMORY_EVENTS
if [ -e /sys/fs/cgroup/memory.events ]; then cat /sys/fs/cgroup/memory.events; else echo unavailable; fi
"""
        code, stdout, stderr = new_md_server.run(client, command, timeout=120)
        if code:
            raise RuntimeError(f"Probe failed with exit code {code}: {stderr[-500:]}")
        print(f"ENDPOINT={new_md_server.endpoint_label()}")
        print(stdout, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
