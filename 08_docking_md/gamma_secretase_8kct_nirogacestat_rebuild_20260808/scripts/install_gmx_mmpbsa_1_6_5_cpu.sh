#!/usr/bin/env bash
set -euo pipefail

GMX_MMPBSA_COMMIT="64e994c71aaff315f3c82dd0852919aecb1ab62e"
CONDA_FORGE_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"

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

"$CONDA_FRONTEND" create -y -p "$PREFIX" --override-channels -c "$CONDA_FORGE_MIRROR" \
  python=3.11.8 ambertools=23.3 gromacs=2023.4 pocl openmpi=4.1.6 c-compiler pip git

"$CONDA_FRONTEND" run -p "$PREFIX" python -m pip install --no-cache-dir \
  mpi4py==4.0.1 \
  numpy==1.26.4 pandas==1.5.3 matplotlib==3.7.3 seaborn==0.11.2 \
  scipy==1.14.1 ParmEd==4.3.0 tqdm==4.67.1

source_archive="$PREFIX/gmx_MMPBSA-${GMX_MMPBSA_COMMIT}.tar.gz"
source_dir="$PREFIX/gmx_MMPBSA-source"
curl -L --http1.1 --retry 10 --retry-all-errors --retry-delay 5 \
  --connect-timeout 20 -fsS \
  "https://codeload.github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tar.gz/${GMX_MMPBSA_COMMIT}" \
  -o "$source_archive"
mkdir "$source_dir"
tar -xzf "$source_archive" --strip-components=1 -C "$source_dir"
sha256sum "$source_archive" > "$source_archive.sha256"

"$CONDA_FRONTEND" run -p "$PREFIX" python -m pip install --no-cache-dir --no-deps \
  "$source_dir"

"$CONDA_FRONTEND" run -p "$PREFIX" python -c \
  "from pathlib import Path; Path(r'$PREFIX/GMX_MMPBSA_SOURCE_COMMIT.txt').write_text('$GMX_MMPBSA_COMMIT\\n', encoding='ascii')"

"$CONDA_FRONTEND" run -p "$PREFIX" python -c \
  "import importlib.metadata as m; assert m.version('gmx-MMPBSA') == '1.6.5'"

echo "isolated CPU endpoint-energy environment created at $PREFIX"
