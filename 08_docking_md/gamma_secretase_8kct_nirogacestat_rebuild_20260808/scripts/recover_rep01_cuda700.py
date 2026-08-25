from __future__ import annotations

import json
import shlex
import time

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
TPR_SHA = "fd11c7287d5670c81ccb44fcb5b4215344726989f66f4f55db33643ba618678f"
AUDIT = f"{BASE}/audit/rep01_cuda700_20260818T024129p0800"


def run(client, command: str, timeout: int = 300) -> tuple[int, str, str]:
    return new_md_server.run(client, command, timeout=timeout)


def main() -> None:
    client = new_md_server.connect()
    try:
        prepare = f"""
set -eu
base={shlex.quote(BASE)}
audit={shlex.quote(AUDIT)}
mkdir -p "$audit/canary"
cp -p "$base/rep01_production.stderr" "$audit/runner.stderr"
cp -p "$base/rep01_production.stdout" "$audit/runner.stdout"
cp -p "$base/rep01/work/production.cpt" "$audit/production.cpt"
tail -n 200 "$base/rep01/work/production.log" > "$audit/production.log.tail"
stat -c '%n %s %y' \
  "$base/rep01/work/production.cpt" \
  "$base/rep01/work/production_prev.cpt" \
  "$base/rep01/work/production.log" \
  "$base/rep01/work/production.edr" \
  "$base/rep01/work/production.xtc" > "$audit/ARTIFACTS.txt"
sha256sum "$audit/runner.stderr" "$audit/runner.stdout" \
  "$audit/production.cpt" "$audit/production.log.tail" \
  > "$audit/EVIDENCE.sha256"
actual=$(sha256sum "$base/rep01/work/production.tpr" | awk '{{print $1}}')
test "$actual" = {shlex.quote(TPR_SHA)}
{shlex.quote(GMX)} dump -cp "$base/rep01/work/production.cpt" \
  > "$audit/CHECKPOINT_DUMP.txt" 2> "$audit/CHECKPOINT_DUMP.stderr"
cp -p "$base/rep01/work/step6.1_equilibration.tpr" "$audit/canary/cuda_health.tpr"
"""
        code, _, stderr = run(client, prepare, timeout=180)
        if code:
            raise RuntimeError(f"audit preparation failed: {stderr[-800:]}")

        canary = f"""
set +e
cd {shlex.quote(AUDIT + '/canary')}
timeout 180s env CUDA_VISIBLE_DEVICES=0 {shlex.quote(GMX)} mdrun \
  -ntmpi 1 -ntomp 16 -pin on -s cuda_health.tpr \
  -deffnm cuda_health -nsteps 10 > mdrun.stdout 2> mdrun.stderr
rc=$?
sha256sum cuda_health.tpr mdrun.stdout mdrun.stderr \
  cuda_health.log cuda_health.edr cuda_health.trr \
  > {shlex.quote(AUDIT + '/CANARY.sha256')} 2>/dev/null
printf '%s' "$rc"
exit 0
"""
        code, stdout, stderr = run(client, canary, timeout=240)
        if code:
            raise RuntimeError(f"canary wrapper failed: {stderr[-800:]}")
        canary_exit = int(stdout.strip())
        if canary_exit != 0:
            print(json.dumps({"status": "canary_failed", "exit": canary_exit, "audit": AUDIT}))
            raise SystemExit(2)

        launch = f"""
set -eu
cd {shlex.quote(BASE)}
if pgrep -x gmx >/dev/null 2>&1; then
  echo active_gmx
  exit 4
fi
nohup env CUDA_VISIBLE_DEVICES=0 GMX_BIN={shlex.quote(GMX)} \
  MDRUN_ARGS='-ntmpi 1 -ntomp 16 -pin on' \
  ./run_frozen_production.sh rep01 {shlex.quote(TPR_SHA)} \
  > rep01_production.stdout 2> rep01_production.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > rep01_production.pid
printf '%s' "$pid"
"""
        code, stdout, stderr = run(client, launch, timeout=120)
        if code:
            raise RuntimeError(f"restart launch failed: {stderr[-800:] or stdout[-800:]}")
        runner_pid = int(stdout.strip())
        time.sleep(15)
        verify = f"""
set -eu
cd {shlex.quote(BASE)}
runner={runner_pid}
test -d /proc/$runner
gmx_pid=$(pgrep -x gmx | head -n 1)
test -n "$gmx_pid"
state=$(awk '/^State:/ {{print $2}}' /proc/$gmx_pid/status)
log_bytes=$(stat -c %s rep01/work/production.log)
printf '%s %s %s' "$gmx_pid" "$state" "$log_bytes"
"""
        code, stdout, stderr = run(client, verify, timeout=120)
        if code:
            raise RuntimeError(f"restart verification failed: {stderr[-800:] or stdout[-800:]}")
        gmx_pid, state, log_bytes = stdout.strip().split()
        print(
            json.dumps(
                {
                    "status": "restarted",
                    "audit": AUDIT,
                    "checkpoint_step": 75914750,
                    "checkpoint_time_ps": 303659.0,
                    "runner_pid": runner_pid,
                    "gmx_pid": int(gmx_pid),
                    "gmx_state": state,
                    "log_bytes": int(log_bytes),
                },
                sort_keys=True,
            )
        )
    finally:
        client.close()


if __name__ == "__main__":
    main()
