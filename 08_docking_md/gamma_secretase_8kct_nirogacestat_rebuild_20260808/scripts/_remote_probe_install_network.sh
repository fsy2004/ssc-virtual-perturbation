#!/usr/bin/env bash
set -u

probe() {
  label=$1
  url=$2
  code=$(curl -L --http1.1 --connect-timeout 10 --max-time 25 -sS -o /dev/null -w '%{http_code}' "$url" 2>&1)
  rc=$?
  printf '%s\trc=%s\thttp=%s\n' "$label" "$rc" "$code"
}

probe conda_forge https://conda.anaconda.org/conda-forge/linux-64/repodata.json.zst
probe tuna_conda_forge https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json.zst
probe bfsu_conda_forge https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json.zst
probe github_repo https://github.com/VachaLab/gorder.git
probe github_codeload https://codeload.github.com/VachaLab/gorder/tar.gz/1beece37dc58a819be0a20b3ec691ef6cade365d
probe rust_static https://static.rust-lang.org/dist/channel-rust-1.82.0.toml.sha256
