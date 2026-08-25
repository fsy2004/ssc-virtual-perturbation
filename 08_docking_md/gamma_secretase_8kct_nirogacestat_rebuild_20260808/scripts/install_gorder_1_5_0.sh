#!/usr/bin/env bash
set -euo pipefail

GORDER_VERSION=1.5.0
GORDER_COMMIT=1beece37dc58a819be0a20b3ec691ef6cade365d
RUST_VERSION=1.87.0

if [[ $# -ne 1 ]]; then
  echo "usage: $0 ABSOLUTE_BASE_PREFIX" >&2
  exit 2
fi
base=$1
[[ "$base" == /* ]] || { echo 'base prefix must be absolute' >&2; exit 2; }
[[ ! -e "$base" ]] || { echo "refusing existing base prefix: $base" >&2; exit 3; }

mkdir -p "$base"
export CARGO_HOME="$base/cargo-home"
export RUSTUP_HOME="$base/rustup-home"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-96}"
rustup_init="$base/rustup-init"

curl --proto '=https' --tlsv1.2 --http1.1 --retry 10 --retry-all-errors --retry-delay 5 \
  --connect-timeout 20 -fsSL \
  https://static.rust-lang.org/rustup/dist/x86_64-unknown-linux-gnu/rustup-init \
  -o "$rustup_init"
chmod 0755 "$rustup_init"
sha256sum "$rustup_init" > "$base/rustup-init.sha256"
"$rustup_init" -y --no-modify-path --profile minimal --default-toolchain "$RUST_VERSION"

source_dir="$base/source"
source_archive="$base/gorder-v${GORDER_VERSION}-source.tar.gz"
curl -L --http1.1 --retry 10 --retry-all-errors --retry-delay 5 \
  --connect-timeout 20 -fsS \
  "https://codeload.github.com/VachaLab/gorder/tar.gz/${GORDER_COMMIT}" \
  -o "$source_archive"
mkdir "$source_dir"
tar -xzf "$source_archive" --strip-components=1 -C "$source_dir"
sha256sum "$source_archive" > "$source_archive.sha256"

"$CARGO_HOME/bin/cargo" install --locked --path "$source_dir" --root "$base/install"
version_output=$($base/install/bin/gorder --version)
[[ "$version_output" == "gorder $GORDER_VERSION" ]] || {
  echo "unexpected gorder version output: $version_output" >&2
  exit 5
}

"$CARGO_HOME/bin/rustc" --version --verbose > "$base/rustc-version.txt"
"$CARGO_HOME/bin/cargo" --version --verbose > "$base/cargo-version.txt"
printf '%s\n' "$version_output" > "$base/gorder-version.txt"
printf '%s\n' "$GORDER_COMMIT" > "$base/GORDER_SOURCE_COMMIT.txt"
sha256sum "$base/install/bin/gorder" > "$base/gorder-binary.sha256"
echo "isolated gorder $GORDER_VERSION toolchain created at $base"
