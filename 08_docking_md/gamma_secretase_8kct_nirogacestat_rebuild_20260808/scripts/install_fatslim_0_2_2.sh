#!/usr/bin/env bash
set -euo pipefail

FATSLIM_VERSION=0.2.2
FATSLIM_COMMIT=ad79df027b62f10edf8e7d65298b13088d46f151
CONDA_FORGE_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"

if [[ $# -ne 1 ]]; then
  echo "usage: $0 ABSOLUTE_BASE_PREFIX" >&2
  exit 2
fi
base=$1
[[ "$base" == /* ]] || { echo 'base prefix must be absolute' >&2; exit 2; }
[[ ! -e "$base" ]] || { echo "refusing existing base prefix: $base" >&2; exit 3; }

mkdir -p "$base"
export CONDA_PKGS_DIRS="$base/conda-pkgs"
env_prefix="$base/env"
source_dir="$base/source"

/root/miniconda3/bin/conda create -y -p "$env_prefix" --override-channels -c "$CONDA_FORGE_MIRROR" \
  python=3.7.12 pip

/root/miniconda3/bin/conda run -p "$env_prefix" python -m pip install --no-cache-dir \
  pip==23.3.2 numpy==1.21.6 cython==0.29.36 pytest==7.4.4

source_archive="$base/fatslim-v${FATSLIM_VERSION}-source.tar.gz"
curl -L --http1.1 --retry 10 --retry-all-errors --retry-delay 5 \
  --connect-timeout 20 -fsS \
  "https://codeload.github.com/FATSLiM/fatslim/tar.gz/${FATSLIM_COMMIT}" \
  -o "$source_archive"
mkdir "$source_dir"
tar -xzf "$source_archive" --strip-components=1 -C "$source_dir"
sha256sum "$source_archive" > "$source_archive.sha256"

/root/miniconda3/bin/conda run -p "$env_prefix" python -m pip install \
  --no-build-isolation --no-deps "$source_dir"
/root/miniconda3/bin/conda run -p "$env_prefix" fatslim version > "$base/fatslim-version.txt"
/root/miniconda3/bin/conda run -p "$env_prefix" fatslim self-test > "$base/fatslim-self-test.txt" 2>&1
grep -Eq '0\.2\.2|version 0\.2\.2' "$base/fatslim-version.txt"

printf '%s\n' "$FATSLIM_COMMIT" > "$base/FATSLIM_SOURCE_COMMIT.txt"
sha256sum "$env_prefix/bin/fatslim" > "$base/fatslim-entrypoint.sha256"
/root/miniconda3/bin/conda list -p "$env_prefix" --explicit > "$base/conda-explicit.txt"
/root/miniconda3/bin/conda list -p "$env_prefix" --json > "$base/conda-list.json"
echo "isolated FATSLiM $FATSLIM_VERSION toolchain created at $base"
