from __future__ import annotations

import sys

import new_md_server


BASE = "/root/autodl-tmp/o6u_md_release_3x500ns_v4"
GMX = "/root/GROMACS-2025.2/bin/gmx"


def main() -> None:
    client = new_md_server.connect()
    try:
        command = f"""
set -u
base={BASE}
echo STDERR
tail -n 120 "$base/rep03_production.stderr"
echo STDOUT
tail -n 80 "$base/rep03_production.stdout"
echo LOGTAIL
tail -n 120 "$base/rep03/work/production.log"
echo FILES
stat -c '%n %s %y' \
  "$base/rep03/work/production.cpt" \
  "$base/rep03/work/production_prev.cpt" \
  "$base/rep03/work/production.log" \
  "$base/rep03/work/production.edr" \
  "$base/rep03/work/production.xtc"
echo TPR
sha256sum "$base/rep03/work/production.tpr"
echo CPT
{GMX} dump -cp "$base/rep03/work/production.cpt" 2>/dev/null \
  | grep -E '^(step|t =|time =)' | head -n 10
echo PREV_CPT
{GMX} dump -cp "$base/rep03/work/production_prev.cpt" 2>/dev/null \
  | grep -E '^(step|t =|time =)' | head -n 10
echo KERNEL_GPU
dmesg --ctime 2>/dev/null | grep -Ei 'NVRM|Xid|oom|killed process' | tail -n 40 || true
"""
        code, stdout, stderr = new_md_server.run(client, command, timeout=300)
        print(stdout)
        if stderr:
            print(stderr, file=sys.stderr)
        raise SystemExit(code)
    finally:
        client.close()


if __name__ == "__main__":
    main()
