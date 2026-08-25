#!/usr/bin/env python3
"""Resident rep03 production monitor (one process until completion).

Loop every INTERVAL seconds: run the hash-bound status check, print one
summary line to stdout, append to the log, and RETRY internally on transient
SSH failures (up to RETRY_MAX attempts, RETRY_WAIT s apart). The process exits
ONLY on:
  - rep03 FINISHED=yes (COMPLETE signal), or
  - a real alarm: frozen hash mismatch, cgroup OOM/oom_kill > 0, disk below
    threshold, gmx gone before completion, or repeated status failures.

Everything healthy is logged and printed; no action is taken otherwise. This
removes the heartbeat-chain fragility: one background job covers the whole
remaining production period.

--once performs a single check with retries (manual validation).
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
STATUS = PROJECT / "scripts" / "status_new_md_run.py"
LOG_DIR = PROJECT / "reports" / "rep03_monitor"
LOG = LOG_DIR / "rep03_monitor.log"
LATEST = LOG_DIR / "latest_status.txt"

INTERVAL = 1800          # seconds between checks (30 min)
STALL_LIMIT = 3          # consecutive no-progress checks before alerting
DISK_FREE_MIN_GB = 40.0  # alert threshold for /root/autodl-tmp free space
MAX_LOOPS = 120          # ~60 h at 30-min intervals
RETRY_MAX = 3            # transient SSH retries per check
RETRY_WAIT = 60          # seconds between retries

FROZEN_HASHES = {
    "5a421f28afee664b5a8919db5f415f1205f35200950117bb3a67fceaba544a98": "release_archive",
    "f442e1411d6f355254d5783903d96f43998a0d5758b469088d91ac8add18aee5": "release_manifest",
    "b036732167b51064b532c454f73137921a4c0cb1ae5a0a4be83d4261d9b0d5ee": "canary_validation",
    "a6e41f920f5af4860b7452c4cbdb2afeed8243bf65fb23b4fd6730e3ebbca4aa": "production_tpr_release",
}


def run_status() -> tuple[bool, str]:
    try:
        r = subprocess.run(
            [sys.executable, str(STATUS)],
            capture_output=True,
            text=True,
            timeout=900,
            cwd=str(PROJECT),
        )
        return r.returncode == 0, (r.stdout or "") + (r.stderr or "")
    except Exception as exc:  # noqa: BLE001
        return False, f"LOCAL_ERROR: {exc!r}"


def parse(text: str) -> dict:
    d: dict = {}
    m = re.search(r"ENDPOINT=(\S+)", text)
    d["endpoint"] = m.group(1) if m else None
    m = re.search(r"REP=rep03 FINISHED=(\w+) LOG_BYTES=(\d+) LOG_MTIME=\"?([^\"]+?)\"?\s+PROGRESS_STEP_TIME_PS=(\d+) ([\d.]+)", text)
    if m:
        d["rep03_finished"] = m.group(1) == "yes"
        d["rep03_step"] = int(m.group(4))
        d["rep03_ps"] = float(m.group(5))
        d["rep03_log_bytes"] = int(m.group(2))
    m = re.search(r"MDRUNS\s*\nPID=(\d+) STATE=(\w+) COMM=gmx", text)
    d["gmx_pid"] = m.group(1) if m else None
    d["gmx_state"] = m.group(2) if m else None
    m = re.search(r"(\d+), (NVIDIA[^,]+), (\d+) MiB, (\d+) MiB, (\d+) %, (\d+), ([\d.]+) W", text)
    if m:
        d["gpu_index"], d["gpu_name"], d["gpu_total"], d["gpu_used"], d["gpu_util"], d["gpu_temp"], d["gpu_power"] = (
            m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6)), float(m.group(7)))
    m = re.search(r"/dev/md0\s+(\d+)G\s+(\d+)G\s+(\d+)G\s+(\d+)%", text)
    if m:
        d["disk_total_g"], d["disk_used_g"], d["disk_free_g"], d["disk_pct"] = (
            int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))
    m = re.search(r"oom_kill (\d+)", text)
    d["oom_kill"] = int(m.group(1)) if m else None
    hashes = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) == 2 and re.fullmatch(r"[0-9a-f]{64}", parts[0]):
            hashes[parts[0]] = parts[1]
    d["hashes"] = hashes
    return d


def log_line(msg: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")


def write_latest(summary: str, alarm: str = "") -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    LATEST.write_text(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {summary}\n{alarm}\n", encoding="utf-8")


def evaluate(d: dict) -> tuple[int, str]:
    """Return (exit_code, alarm_text); 0 = healthy or COMPLETE handled separately."""
    seen = {h: f for h, f in FROZEN_HASHES.items() if h in d.get("hashes", {})}
    if len(seen) != len(FROZEN_HASHES):
        return 3, f"ALERT: integrity hash set changed; expected {len(FROZEN_HASHES)} frozen hashes, saw {len(seen)}."
    if d.get("oom_kill", 0) > 0:
        return 4, f"ALERT: cgroup oom_kill={d.get('oom_kill')} > 0."
    if d.get("disk_free_g", 999) is not None and d.get("disk_free_g", 999) < DISK_FREE_MIN_GB:
        return 5, f"ALERT: free disk {d.get('disk_free_g')}G < {DISK_FREE_MIN_GB}G."
    if d.get("rep03_finished"):
        return 0, ""  # COMPLETE handled by callers
    if d.get("gmx_pid") is None or d.get("gmx_state") != "R":
        return 6, "ALERT: gmx process not running while rep03 incomplete."
    return 0, ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interval", type=int, default=INTERVAL, help="seconds between checks")
    parser.add_argument("--once", action="store_true", help="single check then exit (validation)")
    args = parser.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    if args.once:
        for attempt in range(1, RETRY_MAX + 1):
            ok, text = run_status()
            if ok and "ENDPOINT=" in text:
                d = parse(text)
                summary = (f"check endpoint={d.get('endpoint')} finished={d.get('rep03_finished')} "
                           f"ps={d.get('rep03_ps')} gmx={d.get('gmx_pid')}/{d.get('gmx_state')} "
                           f"gpu={d.get('gpu_util')}% temp={d.get('gpu_temp')}C "
                           f"disk_free={d.get('disk_free_g')}G oom_kill={d.get('oom_kill')}")
                log_line(summary)
                print(summary)
                code, alarm = evaluate(d)
                if alarm:
                    print(alarm)
                return code
            log_line(f"attempt{attempt} failed: {text[:200]}")
            if attempt < RETRY_MAX:
                time.sleep(RETRY_WAIT)
        print("ALERT: status check failed after retries")
        return 2

    log_line(f"RESIDENT MONITOR START interval={args.interval}s stall_limit={STALL_LIMIT} retry={RETRY_MAX}x{RETRY_WAIT}s")
    print(f"RESIDENT MONITOR START interval={args.interval}s", flush=True)
    failures = 0
    stall = 0
    last_ps: float | None = None
    for loop in range(1, MAX_LOOPS + 1):
        d: dict = {}
        for attempt in range(1, RETRY_MAX + 1):
            ok, text = run_status()
            if ok and "ENDPOINT=" in text:
                d = parse(text)
                break
            log_line(f"check#{loop} attempt{attempt} failed: {text[:200]}")
            if attempt < RETRY_MAX:
                time.sleep(RETRY_WAIT)
        else:
            failures += 1
            if failures >= 2:
                msg = "ALERT: status check failed twice consecutively after retries."
                log_line(msg)
                print(msg, flush=True)
                return 2
            log_line(f"check#{loop}: transient failure (count={failures}); continuing")
            time.sleep(min(600, args.interval))
            continue
        failures = 0
        d_ps = d.get("rep03_ps")
        summary = (f"check#{loop} endpoint={d.get('endpoint')} finished={d.get('rep03_finished')} "
                   f"ps={d_ps} gmx={d.get('gmx_pid')}/{d.get('gmx_state')} "
                   f"gpu={d.get('gpu_name')} util={d.get('gpu_util')}% temp={d.get('gpu_temp')}C "
                   f"disk_free={d.get('disk_free_g')}G oom_kill={d.get('oom_kill')}")
        log_line(summary)
        write_latest(summary)
        print(summary, flush=True)

        code, alarm = evaluate(d)
        if code != 0:
            log_line(alarm)
            write_latest(summary, alarm)
            print(alarm, flush=True)
            return code
        if d.get("rep03_finished"):
            print(f"COMPLETE: rep03 finished at {d_ps} ps.", flush=True)
            return 0

        if last_ps is not None and d_ps == last_ps:
            stall += 1
            if stall >= STALL_LIMIT:
                msg = f"ALERT: rep03 progress stalled ({stall} cycles at {d_ps} ps)."
                log_line(msg)
                print(msg, flush=True)
                return 7
        else:
            stall = 0
        last_ps = d_ps

        time.sleep(args.interval)
    print("ALERT: monitor hit MAX_LOOPS without completion.")
    return 8


if __name__ == "__main__":
    sys.exit(main())
