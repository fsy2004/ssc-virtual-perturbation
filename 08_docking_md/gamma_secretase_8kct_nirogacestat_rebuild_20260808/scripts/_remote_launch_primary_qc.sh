#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
manifest="$root/config/primary_postprocessing_manifest.approved.json"
expected_manifest_sha="0862ac9b1cad7f75a27e0dbda61d6eb4757f312c15ed7952bbe1da55e93818cb"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
runtime_dir="$root/audit/postproduction_runtime"
qc_parent="$root/analysis/primary_postprocessing"

mkdir -p "$runtime_dir" "$qc_parent"
cd "$root"

actual_manifest_sha="$(sha256sum "$manifest" | awk '{print $1}')"
if [[ "$actual_manifest_sha" != "$expected_manifest_sha" ]]; then
  echo "REFUSE: approved manifest SHA drift: $actual_manifest_sha" >&2
  exit 11
fi

if grep -R -n -E 'TODO|REPLACE_ME' "$manifest" "$manifest.sha256" >/tmp/primary_manifest_placeholder_hits.txt 2>/dev/null; then
  echo "REFUSE: approved manifest still contains placeholders" >&2
  cat /tmp/primary_manifest_placeholder_hits.txt >&2
  exit 12
fi

if [[ ! -x "$py" ]]; then
  echo "REFUSE: expected analysis Python is not executable: $py" >&2
  exit 13
fi

"$py" - <<'PY' >/tmp/primary_qc_existing_processes.txt
import os
from pathlib import Path

patterns = (
    "analyze_primary_structure_mdanalysis.py",
    "analyze_membrane_qc_mdanalysis.py",
    "gmx_energy_qc.py",
    "validate_primary_postprocessing.py",
)
skip = {os.getpid(), os.getppid()}
for proc in Path("/proc").iterdir():
    if not proc.name.isdigit():
        continue
    pid = int(proc.name)
    if pid in skip:
        continue
    try:
        cmd = proc.joinpath("cmdline").read_bytes().replace(b"\0", b" ").decode("utf-8", "replace").strip()
    except OSError:
        continue
    if cmd and any(p in cmd for p in patterns):
        print(f"{pid} {cmd}")
PY
if [[ -s /tmp/primary_qc_existing_processes.txt ]]; then
  echo "REFUSE: primary QC-related process already active" >&2
  cat /tmp/primary_qc_existing_processes.txt >&2
  exit 14
fi

run_id="$(date -u +%Y%m%dT%H%M%SZ)_primary_qc_v1"
out="$qc_parent/$run_id"
if [[ -e "$out" ]]; then
  echo "REFUSE: output root already exists: $out" >&2
  exit 15
fi

runner="$runtime_dir/run_${run_id}.sh"
runner_log="$runtime_dir/run_${run_id}.log"
runner_status="$runtime_dir/run_${run_id}.status.tsv"
runner_json="$runtime_dir/run_${run_id}.launch.json"

cat > "$runner" <<'RUNNER'
#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
manifest="$root/config/primary_postprocessing_manifest.approved.json"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
out="${O6U_PRIMARY_QC_OUT:?missing output root}"
status="${O6U_PRIMARY_QC_STATUS:?missing status path}"

mkdir -p "$out"/logs "$out"/runtime
cd "$root"

write_status() {
  local stage="$1"
  local state="$2"
  local code="$3"
  local ts
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '%s\t%s\t%s\t%s\n' "$ts" "$stage" "$state" "$code" >> "$status"
}

printf 'timestamp_utc\tstage\tstate\texit_code\n' > "$status"
write_status runner started 0

set +e
"$py" scripts/analyze_primary_structure_mdanalysis.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/structural.log" 2>&1 &
struct_pid=$!

"$py" scripts/gmx_energy_qc.py \
  --manifest "$manifest" \
  --output-root "$out" \
  --mode extract \
  > "$out/logs/energy.log" 2>&1 &
energy_pid=$!

wait "$struct_pid"
struct_code=$?
write_status structural finished "$struct_code"

wait "$energy_pid"
energy_code=$?
write_status energy finished "$energy_code"
set -e

if [[ "$struct_code" -ne 0 || "$energy_code" -ne 0 ]]; then
  write_status runner failed_pre_membrane 20
  exit 20
fi

set +e
"$py" scripts/analyze_membrane_qc_mdanalysis.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/membrane.log" 2>&1
membrane_code=$?
write_status membrane finished "$membrane_code"
set -e

if [[ "$membrane_code" -ne 0 ]]; then
  write_status runner failed_pre_validation 30
  exit 30
fi

set +e
"$py" scripts/validate_primary_postprocessing.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/validate_primary.log" 2>&1
validate_code=$?
write_status validate_primary finished "$validate_code"
set -e

if [[ "$validate_code" -ne 0 ]]; then
  write_status runner validation_failed "$validate_code"
  exit "$validate_code"
fi

write_status runner complete 0
RUNNER

chmod 700 "$runner"
mkdir -p "$out"

cat > "$runner_json" <<EOF
{
  "status": "launched",
  "run_id": "$run_id",
  "root": "$root",
  "manifest": "$manifest",
  "manifest_sha256": "$actual_manifest_sha",
  "output_root": "$out",
  "runner": "$runner",
  "runner_log": "$runner_log",
  "status_tsv": "$runner_status",
  "policy": "structural_and_energy_parallel_then_membrane_then_primary_validation",
  "raw_trajectory_policy": "read_only_no_overwrite",
  "launched_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
sha256sum "$runner_json" > "$runner_json.sha256"

O6U_PRIMARY_QC_OUT="$out" O6U_PRIMARY_QC_STATUS="$runner_status" \
  nohup bash "$runner" > "$runner_log" 2>&1 &
runner_pid="$!"

cat <<EOF
PRIMARY_QC_LAUNCHED
run_id=$run_id
pid=$runner_pid
output_root=$out
runner_log=$runner_log
status_tsv=$runner_status
manifest_sha256=$actual_manifest_sha
EOF
