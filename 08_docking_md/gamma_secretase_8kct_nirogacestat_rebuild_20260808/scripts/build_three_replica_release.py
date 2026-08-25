#!/usr/bin/env python3
"""Build a fail-closed three-realization 8KCT--O6U GROMACS release."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any


REPLICA_IDS = ("rep01", "rep02", "rep03")
EQUILIBRATION_FILES = tuple(f"step6.{index}_equilibration.mdp" for index in range(1, 7))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_mdp(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().lower().replace("_", "-")] = value.strip()
    return values


def validate_realizations(protocol: dict[str, Any]) -> None:
    realizations = protocol.get("realizations")
    if not isinstance(realizations, list) or [item.get("id") for item in realizations] != list(REPLICA_IDS):
        raise ValueError("realizations must be exactly rep01, rep02, rep03")
    velocity = [int(item["velocity_seed"]) for item in realizations]
    thermostat = [int(item["thermostat_seed"]) for item in realizations]
    if len(set(velocity)) != 3:
        raise ValueError("velocity seeds must be distinct")
    if len(set(thermostat)) != 3:
        raise ValueError("thermostat seeds must be distinct")
    if any(value <= 0 or value >= 2**31 for value in velocity + thermostat):
        raise ValueError("all seeds must be positive signed 32-bit integers")


def validate_production_mdp(text: str) -> dict[str, str]:
    values = parse_mdp(text)
    expected = {
        "integrator": "md",
        "dt": "0.004",
        "nsteps": "125000000",
        "tcoupl": "v-rescale",
        "pcoupl": "C-rescale",
        "pcoupltype": "semiisotropic",
        "constraints": "h-bonds",
        "constraint-algorithm": "LINCS",
        "continuation": "yes",
        "gen-vel": "no",
    }
    for key, wanted in expected.items():
        actual = values.get(key)
        if actual is None or actual.lower() != wanted.lower():
            if key == "gen-vel":
                raise ValueError("gen-vel must be no in production")
            raise ValueError(f"{key} must be {wanted}, found {actual}")
    temperatures = [float(value) for value in values.get("ref-t", "").split()]
    if temperatures != [303.15, 303.15, 303.15]:
        raise ValueError("production ref-t must be exactly three 303.15 K values")
    if values.get("tc-grps", "").split() != ["SOLU", "MEMB", "SOLV"]:
        raise ValueError("production tc-grps must be SOLU MEMB SOLV")
    return values


def _replace_assignment(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=\s*)[^;\r\n]*(.*)$")
    rendered, count = pattern.subn(lambda match: f"{key:<24}= {value}{match.group(2)}", text, count=1)
    if count != 1:
        raise ValueError(f"expected exactly one {key} assignment")
    return rendered


def render_step61(source: str, *, velocity_seed: int, thermostat_seed: int) -> str:
    values = parse_mdp(source)
    required = {
        "integrator": "md", "dt": "0.001", "nsteps": "125000",
        "tcoupl": "v-rescale", "constraints": "h-bonds", "gen-vel": "yes",
        "gen-temp": "303.15",
    }
    for key, wanted in required.items():
        if values.get(key, "").lower() != wanted.lower():
            raise ValueError(f"step6.1 {key} must be {wanted}")
    rendered = _replace_assignment(source, "gen-seed", str(velocity_seed))
    if "ld-seed" in parse_mdp(rendered):
        rendered = _replace_assignment(rendered, "ld-seed", str(thermostat_seed))
    else:
        rendered = rendered.rstrip() + f"\nld-seed                 = {thermostat_seed}\n"
    return rendered


def ensure_new_destination(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"release destination already exists: {path}")


def _write_runner(destination: Path) -> None:
    runner = r'''#!/usr/bin/env bash
set -euo pipefail
rep="${1:?usage: run_replica.sh rep01|rep02|rep03 [equilibrate|produce]}"
mode="${2:-equilibrate}"
case "$rep" in rep01|rep02|rep03) ;; *) echo "invalid replica: $rep" >&2; exit 2;; esac
case "$mode" in equilibrate|produce) ;; *) echo "invalid mode: $mode" >&2; exit 2;; esac
gmx_bin="${GMX_BIN:-gmx}"
read -r -a mdrun_args <<< "${MDRUN_ARGS:-}"
root="$(cd "$(dirname "$0")" && pwd)"
common="$root/common"
work="$root/$rep/work"
mkdir -p "$work"

scan_log() {
  local log="$1"
  if grep -Ein 'LINCS WARNING|Too many LINCS warnings|constraint warning|SETTLE.*(error|constraint)|(^|[^A-Za-z])NaN([^A-Za-z]|$)|Fatal error|Segmentation fault' "$log" >/dev/null; then
    echo "blocking dynamics warning in $log" >&2
    exit 3
  fi
}

if [[ "$mode" == "equilibrate" ]]; then
  for stage in 1 2 3 4 5 6; do
    prefix="step6.${stage}_equilibration"
    if [[ "$stage" == 1 ]]; then
      "$gmx_bin" grompp -f "$root/$rep/${prefix}.mdp" -o "$work/${prefix}.tpr" \
        -c "$common/minimized.gro" -r "$common/step5_input.gro" \
        -p "$common/topol.top" -n "$common/index.ndx" -maxwarn 0
    else
      prev="step6.$((stage-1))_equilibration"
      test -s "$work/${prev}.cpt"
      "$gmx_bin" grompp -f "$root/$rep/${prefix}.mdp" -o "$work/${prefix}.tpr" \
        -c "$work/${prev}.gro" -t "$work/${prev}.cpt" -r "$common/step5_input.gro" \
        -p "$common/topol.top" -n "$common/index.ndx" -maxwarn 0
    fi
    restart_args=()
    if [[ -s "$work/${prefix}.cpt" ]]; then
      restart_args=(-cpi "$work/${prefix}.cpt" -append)
    fi
    "$gmx_bin" mdrun "${mdrun_args[@]}" -deffnm "$work/${prefix}" "${restart_args[@]}"
    scan_log "$work/${prefix}.log"
  done
  exit 0
fi

test -s "$work/step6.6_equilibration.cpt"
"$gmx_bin" grompp -f "$root/$rep/production_500ns.mdp" -o "$work/production.tpr" \
  -c "$work/step6.6_equilibration.gro" -t "$work/step6.6_equilibration.cpt" \
  -p "$common/topol.top" -n "$common/index.ndx" -maxwarn 0
restart_args=()
if [[ -s "$work/production.cpt" ]]; then
  restart_args=(-cpi "$work/production.cpt" -append)
fi
"$gmx_bin" mdrun "${mdrun_args[@]}" -deffnm "$work/production" "${restart_args[@]}"
scan_log "$work/production.log"
'''
    path = destination / "run_replica.sh"
    path.write_text(runner, encoding="utf-8", newline="\n")


def build_release(source: Path, minimized: Path, protocol_path: Path, production_path: Path, destination: Path) -> None:
    ensure_new_destination(destination)
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    validate_realizations(protocol)
    production_text = production_path.read_text(encoding="utf-8")
    validate_production_mdp(production_text)
    expected = protocol["system"]
    if sha256(source / "step5_input.gro") != expected["step5_input_gro_sha256"]:
        raise ValueError("step5_input.gro hash mismatch")
    if sha256(minimized) != expected["minimized_gro_sha256"]:
        raise ValueError("minimized GRO hash mismatch")
    for name in ("topol.top", "index.ndx", *EQUILIBRATION_FILES):
        if not (source / name).is_file():
            raise FileNotFoundError(source / name)

    destination.mkdir(parents=True)
    common = destination / "common"
    common.mkdir()
    for item in source.iterdir():
        target = common / item.name
        if item.is_dir():
            shutil.copytree(item, target)
        elif item.is_file():
            shutil.copy2(item, target)
    shutil.copy2(minimized, common / "minimized.gro")
    shutil.copy2(protocol_path, destination / "production_protocol.json")
    docs_source = Path(__file__).resolve().parent.parent / "docs"
    for name in ("GPU_LAUNCH_README.md", "RELEASE_LIMITATIONS.md"):
        source_doc = docs_source / name
        if not source_doc.is_file():
            raise FileNotFoundError(source_doc)
        shutil.copy2(source_doc, destination / name)

    source_step61 = (source / EQUILIBRATION_FILES[0]).read_text(encoding="utf-8")
    by_id = {item["id"]: item for item in protocol["realizations"]}
    for replica_id in REPLICA_IDS:
        replica = destination / replica_id
        replica.mkdir()
        item = by_id[replica_id]
        rendered = render_step61(
            source_step61,
            velocity_seed=int(item["velocity_seed"]),
            thermostat_seed=int(item["thermostat_seed"]),
        )
        (replica / EQUILIBRATION_FILES[0]).write_text(rendered, encoding="utf-8", newline="\n")
        for name in EQUILIBRATION_FILES[1:]:
            shutil.copy2(source / name, replica / name)
        shutil.copy2(production_path, replica / "production_500ns.mdp")
    _write_runner(destination)

    artifacts = []
    for path in sorted(destination.rglob("*")):
        if path.is_file():
            artifacts.append({
                "path": path.relative_to(destination).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    manifest = {
        "schema_version": "1.0",
        "protocol_id": protocol["protocol_id"],
        "archive_sha256": protocol["system"]["archive_sha256"],
        "realizations": protocol["realizations"],
        "artifacts": artifacts,
    }
    (destination / "RELEASE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-gromacs", type=Path, required=True)
    parser.add_argument("--minimized-gro", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--production-mdp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_release(
        args.source_gromacs.resolve(), args.minimized_gro.resolve(),
        args.protocol.resolve(), args.production_mdp.resolve(), args.output.resolve()
    )


if __name__ == "__main__":
    main()
