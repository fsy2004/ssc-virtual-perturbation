#!/usr/bin/env python3
"""Deploy the sealed v3 release and start rep01 equilibration on the current GPU."""

from __future__ import annotations

import hashlib
import os
import time
from pathlib import Path

import new_md_server


LOCAL_ARCHIVE = Path(os.environ["LOCALAPPDATA"]) / "Temp" / "o6u_md_release_3x500ns_v4.tgz"
ARCHIVE_SHA = "5a421f28afee664b5a8919db5f415f1205f35200950117bb3a67fceaba544a98"
MANIFEST_SHA = "f442e1411d6f355254d5783903d96f43998a0d5758b469088d91ac8add18aee5"
REMOTE_ARCHIVE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4.tgz"
BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked(client, command: str, timeout: int = 300) -> str:
    code, stdout, stderr = new_md_server.run(client, command, timeout=timeout)
    if code:
        raise RuntimeError(f"Remote command failed with exit code {code}: {stderr[-800:]}")
    return stdout


def main() -> None:
    if sha256(LOCAL_ARCHIVE) != ARCHIVE_SHA:
        raise ValueError("Local v3 archive hash mismatch")
    client = new_md_server.connect()
    try:
        present = checked(client, f"if [ -s {REMOTE_ARCHIVE} ]; then echo yes; else echo no; fi").strip()
        if present != "yes":
            upload = REMOTE_ARCHIVE + ".uploading"
            with client.open_sftp() as sftp:
                sftp.put(str(LOCAL_ARCHIVE), upload)
                sftp.rename(upload, REMOTE_ARCHIVE)

        archive_hash = checked(client, f"sha256sum {REMOTE_ARCHIVE}").split()[0]
        if archive_hash != ARCHIVE_SHA:
            raise ValueError("Remote v3 archive hash mismatch")

        checked(
            client,
            f"""
set -eu
cd /root/autodl-tmp
tar -tzf {REMOTE_ARCHIVE} >/dev/null
if [ ! -d {BASE} ]; then tar -xzf {REMOTE_ARCHIVE}; fi
test "$(sha256sum {BASE}/RELEASE_MANIFEST.json | awk '{{print $1}}')" = {MANIFEST_SHA}
/root/miniconda3/bin/python - <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path('{BASE}')
data=json.loads((root/'RELEASE_MANIFEST.json').read_text())
bad=[]
for item in data['artifacts']:
    path=root/item['path']
    digest=hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    if digest != item['sha256'] or (path.is_file() and path.stat().st_size != item['bytes']):
        bad.append(item['path'])
print('ARTIFACTS',len(data['artifacts']))
print('MISMATCHES',len(bad))
if bad:
    print('BAD',*bad)
    sys.exit(3)
PY
chmod 0755 {BASE}/run_replica.sh
""",
            timeout=300,
        )

        launched = checked(
            client,
            f"""
set -eu
cd {BASE}
if [ -s rep01_equilibrate.pid ]; then
  prior=$(cat rep01_equilibrate.pid)
  if [ -d /proc/$prior ]; then echo already_running; exit 4; fi
fi
nohup env CUDA_VISIBLE_DEVICES=0 GMX_BIN={GMX} MDRUN_ARGS='-ntmpi 1 -ntomp 16 -pin on' \
  ./run_replica.sh rep01 equilibrate > rep01_equilibrate.stdout 2> rep01_equilibrate.stderr < /dev/null &
pid=$!
printf '%s\n' "$pid" > rep01_equilibrate.pid
echo RUNNER_PID=$pid
""",
        )
        print(f"ENDPOINT={new_md_server.endpoint_label()}")
        print(f"ARCHIVE_SHA256={archive_hash}")
        print(launched, end="")
        time.sleep(12)
        status = checked(
            client,
            f"""
set -eu
cd {BASE}
runner=$(cat rep01_equilibrate.pid)
if [ -d /proc/$runner ]; then
  echo RUNNER_PID=$runner STATE=$(awk '/^State:/{{print $2}}' /proc/$runner/status) COMM=$(cat /proc/$runner/comm)
else
  echo RUNNER_PID=$runner STATE=exited
fi
for p in $(pgrep -x gmx 2>/dev/null || true); do
  ancestor=$(awk '/^PPid:/{{print $2}}' /proc/$p/status)
  is_child=0
  while [ "$ancestor" -gt 1 ] && [ -d /proc/$ancestor ]; do
    if [ "$ancestor" = "$runner" ]; then is_child=1; break; fi
    ancestor=$(awk '/^PPid:/{{print $2}}' /proc/$ancestor/status)
  done
  if [ "$is_child" -eq 1 ]; then
    echo CHILD_PID=$p STATE=$(awk '/^State:/{{print $2}}' /proc/$p/status) COMM=$(cat /proc/$p/comm)
  fi
done
tail -20 rep01_equilibrate.stderr 2>/dev/null || true
echo GPU_PROCESSES
nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader 2>/dev/null || true
""",
        )
        print(status, end="")
    finally:
        client.close()


if __name__ == "__main__":
    main()
