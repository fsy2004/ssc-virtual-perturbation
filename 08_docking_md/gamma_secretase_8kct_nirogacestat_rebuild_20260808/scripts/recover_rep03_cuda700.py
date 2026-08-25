from __future__ import annotations

import json
import shlex
import time

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"
TPR_SHA = "abe512a4971c6cc26a61c3d9fbea39df8b40074265e1d30b11a488bd3ffac9ad"
AUDIT = f"{BASE}/audit/rep03_cuda700_20260821T153927p0800"
EXPECTED_CHECKPOINT_STEP = 62981800
EXPECTED_CHECKPOINT_TIME_PS = 251927.2


def run(client, command: str, timeout: int = 300) -> tuple[int, str, str]:
    return new_md_server.run(client, command, timeout=timeout)


def main() -> None:
    client = new_md_server.connect()
    try:
        prepare = f"""
set -eu
base={shlex.quote(BASE)}
audit={shlex.quote(AUDIT)}
test ! -e "$audit"
mkdir -p "$audit/canary"
cp -p "$base/rep03_production.stderr" "$audit/runner.stderr"
cp -p "$base/rep03_production.stdout" "$audit/runner.stdout"
cp -p "$base/rep03/work/production.cpt" "$audit/production.cpt"
cp -p "$base/rep03/work/production_prev.cpt" "$audit/production_prev.cpt"
tail -n 240 "$base/rep03/work/production.log" > "$audit/production.log.tail"
stat -c '%n %s %y' \
  "$base/rep03/work/production.cpt" \
  "$base/rep03/work/production_prev.cpt" \
  "$base/rep03/work/production.log" \
  "$base/rep03/work/production.edr" \
  "$base/rep03/work/production.xtc" > "$audit/ARTIFACTS.txt"
sha256sum "$audit/runner.stderr" "$audit/runner.stdout" \
  "$audit/production.cpt" "$audit/production_prev.cpt" \
  "$audit/production.log.tail" > "$audit/EVIDENCE.sha256"
actual=$(sha256sum "$base/rep03/work/production.tpr" | awk '{{print $1}}')
test "$actual" = {shlex.quote(TPR_SHA)}
{shlex.quote(GMX)} dump -cp "$base/rep03/work/production.cpt" \
  > "$audit/CHECKPOINT_DUMP.txt" 2> "$audit/CHECKPOINT_DUMP.stderr"
grep -q '^step = {EXPECTED_CHECKPOINT_STEP}$' "$audit/CHECKPOINT_DUMP.txt"
grep -q '^t = {EXPECTED_CHECKPOINT_TIME_PS:.6f}$' "$audit/CHECKPOINT_DUMP.txt"
dmesg --ctime 2>/dev/null | grep -Ei 'NVRM|Xid|oom|killed process' \
  > "$audit/KERNEL_GPU_EVENTS.txt" || true
cp -p "$base/rep03/work/step6.1_equilibration.tpr" "$audit/canary/cuda_health.tpr"
"""
        code, _, stderr = run(client, prepare, timeout=180)
        if code:
            raise RuntimeError(f"audit preparation failed: {stderr[-1200:]}")

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
            raise RuntimeError(f"canary wrapper failed: {stderr[-1200:]}")
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
  ./run_frozen_production.sh rep03 {shlex.quote(TPR_SHA)} \
  > rep03_production.stdout 2> rep03_production.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > rep03_production.pid
printf '%s' "$pid"
"""
        code, stdout, stderr = run(client, launch, timeout=120)
        if code:
            raise RuntimeError(f"restart launch failed: {stderr[-1200:] or stdout[-1200:]}")
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
log_bytes=$(stat -c %s rep03/work/production.log)
grep -q 'continuing from step {EXPECTED_CHECKPOINT_STEP}, {EXPECTED_CHECKPOINT_TIME_PS:.1f} ps' rep03_production.stderr
grep -q -- '-cpi .*production.cpt -append' rep03_production.stderr
printf '%s %s %s' "$gmx_pid" "$state" "$log_bytes"
"""
        code, stdout, stderr = run(client, verify, timeout=120)
        if code:
            raise RuntimeError(f"restart verification failed: {stderr[-1200:] or stdout[-1200:]}")
        gmx_pid, state, log_bytes = stdout.strip().split()
        print(
            json.dumps(
                {
                    "status": "restarted",
                    "audit": AUDIT,
                    "checkpoint_step": EXPECTED_CHECKPOINT_STEP,
                    "checkpoint_time_ps": EXPECTED_CHECKPOINT_TIME_PS,
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
