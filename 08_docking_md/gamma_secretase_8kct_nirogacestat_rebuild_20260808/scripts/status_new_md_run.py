#!/usr/bin/env python3
"""Sanitized integrity and progress monitor for the current GPU MD release."""

from __future__ import annotations

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
ARCHIVE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4.tgz"


def main() -> None:
    client = new_md_server.connect()
    try:
        command = f"""
set -eu
cd {BASE}
echo ENDPOINT={new_md_server.endpoint_label()}
echo INTEGRITY
sha256sum {ARCHIVE} RELEASE_MANIFEST.json CANARY_VALIDATION.json
echo RUNNERS
for rep in rep01 rep02 rep03; do
  file=${{rep}}_equilibrate.pid
  [ -s "$file" ] || continue
  pid=$(cat "$file")
  if [ -d /proc/$pid ]; then
    echo REP=$rep PID=$pid STATE=$(awk '/^State:/{{print $2}}' /proc/$pid/status) COMM=$(cat /proc/$pid/comm)
  else
    echo REP=$rep PID=$pid STATE=exited
  fi
done
echo PRODUCTION_RUNNERS
for rep in rep01 rep02 rep03; do
  file=${{rep}}_production.pid
  [ -s "$file" ] || continue
  pid=$(cat "$file")
  if [ -d /proc/$pid ]; then
    echo REP=$rep PID=$pid STATE=$(awk '/^State:/{{print $2}}' /proc/$pid/status) COMM=$(cat /proc/$pid/comm)
  else
    echo REP=$rep PID=$pid STATE=exited
  fi
done
echo MDRUNS
for p in $(pgrep -x gmx 2>/dev/null || true); do
  echo PID=$p STATE=$(awk '/^State:/{{print $2}}' /proc/$p/status) COMM=$(cat /proc/$p/comm)
done
echo STAGES
for rep in rep01 rep02 rep03; do
  work=$rep/work
  [ -d "$work" ] || continue
  echo REP=$rep
  for stage in 1 2 3 4 5 6; do
    prefix=step6.${{stage}}_equilibration
    if [ -s "$work/${{prefix}}.log" ]; then
      finished=no
      grep -q 'Finished mdrun on rank 0' "$work/${{prefix}}.log" && finished=yes
      step=$(grep -E '^(Writing checkpoint, step|Step +Time)' "$work/${{prefix}}.log" | tail -1 || true)
      bytes=$(stat -c %s "$work/${{prefix}}.log")
      mtime=$(stat -c %y "$work/${{prefix}}.log")
      echo STAGE=$stage FINISHED=$finished LOG_BYTES=$bytes LOG_MTIME="$mtime" PROGRESS="$step"
    fi
  done
done
echo PRODUCTION
if [ -s PRODUCTION_TPR_RELEASE.json ]; then
  sha256sum PRODUCTION_TPR_RELEASE.json
fi
for rep in rep01 rep02 rep03; do
  work=$rep/work
  log=$work/production.log
  [ -s "$log" ] || continue
  finished=no
  grep -q 'Finished mdrun on rank 0' "$log" && finished=yes
  progress=$(awk '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+([.][0-9]+)?$/ {{line=$1 " " $2}} END {{print line}}' "$log")
  time_ps=$(printf '%s\n' "$progress" | awk '{{print $2}}')
  bytes=$(stat -c %s "$log")
  mtime=$(stat -c %y "$log")
  cpt_bytes=0; xtc_bytes=0; edr_bytes=0
  [ -e "$work/production.cpt" ] && cpt_bytes=$(stat -c %s "$work/production.cpt")
  [ -e "$work/production.xtc" ] && xtc_bytes=$(stat -c %s "$work/production.xtc")
  [ -e "$work/production.edr" ] && edr_bytes=$(stat -c %s "$work/production.edr")
  gates="100ps:pending,1ns:pending,5ns:pending,10ns:pending"
  if [ -n "$time_ps" ]; then
    gates=$(awk -v t="$time_ps" 'BEGIN {{printf "100ps:%s,1ns:%s,5ns:%s,10ns:%s", (t>=100?"reached":"pending"), (t>=1000?"reached":"pending"), (t>=5000?"reached":"pending"), (t>=10000?"reached":"pending")}}')
  fi
  echo REP=$rep FINISHED=$finished LOG_BYTES=$bytes LOG_MTIME="$mtime" PROGRESS_STEP_TIME_PS="$progress" GATES=$gates CPT_BYTES=$cpt_bytes XTC_BYTES=$xtc_bytes EDR_BYTES=$edr_bytes
done
echo BLOCKING_SCAN
grep -ERin 'LINCS WARNING|Too many LINCS warnings|constraint warning|SETTLE.*(error|constraint)|(^|[^A-Za-z])NaN([^A-Za-z]|$)|Fatal error|Segmentation fault' rep*/work/*.log rep*_equilibrate.stderr 2>/dev/null || true
echo GPU
nvidia-smi --query-gpu=index,name,memory.total,memory.free,utilization.gpu,temperature.gpu,power.draw --format=csv,noheader
nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader 2>/dev/null || true
echo FILESYSTEMS
df -h / /root/autodl-tmp | tail -n +2
echo MEMORY_EVENTS
if [ -e /sys/fs/cgroup/memory.events ]; then cat /sys/fs/cgroup/memory.events; else echo unavailable; fi
"""
        code, stdout, stderr = new_md_server.run(client, command, timeout=120)
        if code:
            raise RuntimeError(f"Status failed with exit code {code}: {stderr[-800:]}")
        print(stdout, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
