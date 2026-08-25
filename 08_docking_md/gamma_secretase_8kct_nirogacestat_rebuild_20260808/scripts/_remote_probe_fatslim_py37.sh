#!/usr/bin/env bash
set -euo pipefail

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
/root/miniconda3/bin/conda search --override-channels \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main \
  'python=3.7.12' --json > "$tmp/python.json"
/root/miniconda3/bin/python - "$tmp/python.json" <<'PY'
import json
import sys

payload = json.load(open(sys.argv[1], encoding="utf-8"))
print([(item.get("version"), item.get("build")) for item in payload.get("python", [])])
PY
