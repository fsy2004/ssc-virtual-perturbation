#!/usr/bin/env bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
prefix=/root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822_v5
audit="$root/audit/toolchain_install/20260822_gmx_mmpbsa_retry_v5/20260822T112115Z"
conda=/root/miniconda3/bin/conda
expected_commit=64e994c71aaff315f3c82dd0852919aecb1ab62e

if [[ -r /proc/2978/cmdline ]] && tr '\0' ' ' < /proc/2978/cmdline | grep -q 'gmx_mmpbsa'; then
  echo 'gmx_MMPBSA installer is still active' >&2
  exit 3
fi
[[ -d "$prefix" ]]
[[ "$(tr -d '\r\n' < "$prefix/GMX_MMPBSA_SOURCE_COMMIT.txt")" == "$expected_commit" ]]
(cd "$prefix" && sha256sum -c "$(basename "$prefix/gmx_MMPBSA-${expected_commit}.tar.gz.sha256")")

"$conda" run -p "$prefix" python -c \
  "import importlib.metadata as m, mpi4py; assert m.version('gmx-MMPBSA') == '1.6.5'; assert m.version('mpi4py') == '4.0.1'"
"$conda" run -p "$prefix" gmx_MMPBSA --version > "$audit/gmx_MMPBSA-version.txt"
"$conda" run -p "$prefix" mpirun --version > "$audit/mpirun-version.txt"
"$conda" run -p "$prefix" mpicc --showme:command > "$audit/mpicc-command.txt"
"$conda" run -p "$prefix" gmx --version > "$audit/gromacs-version.txt"
"$conda" list -p "$prefix" --json > "$audit/conda-list.json"
"$conda" list -p "$prefix" --explicit > "$audit/conda-explicit.txt"
"$conda" run -p "$prefix" python -m pip freeze --all > "$audit/pip-freeze.txt"

: > "$audit/ENTRYPOINTS.sha256"
for executable in gmx_MMPBSA mpicc mpirun gmx MMPBSA.py; do
  if [[ -f "$prefix/bin/$executable" ]]; then
    sha256sum "$prefix/bin/$executable" >> "$audit/ENTRYPOINTS.sha256"
  fi
done
sha256sum \
  "$prefix/GMX_MMPBSA_SOURCE_COMMIT.txt" \
  "$prefix/gmx_MMPBSA-${expected_commit}.tar.gz" \
  "$audit/gmx_MMPBSA-version.txt" \
  "$audit/mpirun-version.txt" \
  "$audit/mpicc-command.txt" \
  "$audit/gromacs-version.txt" \
  "$audit/conda-list.json" \
  "$audit/conda-explicit.txt" \
  "$audit/pip-freeze.txt" \
  "$audit/ENTRYPOINTS.sha256" > "$audit/ACCEPTANCE_ARTIFACTS.sha256"
printf 'status\tpass\nprefix\t%s\nsource_commit\t%s\ngmx_mmpbsa_version\t1.6.5\nmpi4py_version\t4.0.1\n' \
  "$prefix" "$expected_commit" > "$audit/TOOLCHAIN_ACCEPTANCE.tsv"
sha256sum "$audit/TOOLCHAIN_ACCEPTANCE.tsv" > "$audit/TOOLCHAIN_ACCEPTANCE.tsv.sha256"
cat "$audit/TOOLCHAIN_ACCEPTANCE.tsv"
