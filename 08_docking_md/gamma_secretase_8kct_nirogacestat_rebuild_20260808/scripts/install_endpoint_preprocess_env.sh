#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 ABSOLUTE_ENV_PREFIX" >&2
  exit 2
fi

PREFIX="$1"
if [[ "$PREFIX" != /* ]]; then
  echo "environment prefix must be absolute" >&2
  exit 2
fi
if [[ -e "$PREFIX" ]]; then
  echo "refusing to overwrite existing environment prefix: $PREFIX" >&2
  exit 3
fi

if [[ -n "${CONDA_EXE:-}" ]]; then
  CONDA_FRONTEND="$CONDA_EXE"
elif command -v micromamba >/dev/null 2>&1; then
  CONDA_FRONTEND="$(command -v micromamba)"
elif command -v conda >/dev/null 2>&1; then
  CONDA_FRONTEND="$(command -v conda)"
else
  echo "conda or micromamba is required" >&2
  exit 4
fi

"$CONDA_FRONTEND" create -y -p "$PREFIX" -c conda-forge python=3.11.8 pip
"$CONDA_FRONTEND" run -p "$PREFIX" python -m pip install --no-cache-dir \
  numpy==1.26.4 MDAnalysis==2.10.0
"$CONDA_FRONTEND" run -p "$PREFIX" python -c \
  "import MDAnalysis, numpy; assert MDAnalysis.__version__ == '2.10.0'; assert numpy.__version__ == '1.26.4'"

echo "isolated endpoint preprocessing environment created at $PREFIX"
