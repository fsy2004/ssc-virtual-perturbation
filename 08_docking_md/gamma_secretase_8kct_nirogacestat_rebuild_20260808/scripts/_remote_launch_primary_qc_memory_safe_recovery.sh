#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
manifest="$root/config/primary_postprocessing_manifest.approved.json"
expected_manifest_sha="0862ac9b1cad7f75a27e0dbda61d6eb4757f312c15ed7952bbe1da55e93818cb"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
runtime_dir="$root/audit/postproduction_runtime"
qc_parent="$root/analysis/primary_postprocessing"
failed_run_id="20260822T144701Z_primary_qc_v1"
failed_status="$runtime_dir/run_${failed_run_id}.status.tsv"
recovery_status="$runtime_dir/recover_${failed_run_id}_after_legend_fix.status.tsv"

mkdir -p "$runtime_dir" "$qc_parent"
cd "$root"

actual_manifest_sha="$(sha256sum "$manifest" | awk '{print $1}')"
if [[ "$actual_manifest_sha" != "$expected_manifest_sha" ]]; then
  echo "REFUSE: approved manifest SHA drift: $actual_manifest_sha" >&2
  exit 51
fi

if ! grep -q $'\tstructural\tfinished\t137' "$failed_status"; then
  echo "REFUSE: expected failed structural exit 137 evidence is absent" >&2
  cat "$failed_status" >&2 || true
  exit 52
fi

if ! grep -q $'\trefuse\tstructural_not_passed\t41' "$recovery_status"; then
  echo "REFUSE: expected prior recovery refusal evidence is absent" >&2
  cat "$recovery_status" >&2 || true
  exit 53
fi

if grep -R -n -E 'TODO|REPLACE_ME' "$manifest" "$manifest.sha256" >/tmp/primary_manifest_placeholder_hits.txt 2>/dev/null; then
  echo "REFUSE: approved manifest still contains placeholders" >&2
  cat /tmp/primary_manifest_placeholder_hits.txt >&2
  exit 54
fi

if [[ ! -x "$py" ]]; then
  echo "REFUSE: expected analysis Python is not executable: $py" >&2
  exit 55
fi

"$py" - <<'PY' >/tmp/primary_qc_existing_processes.txt
import os
from pathlib import Path

patterns = (
    "run_primary_structure_memory_safe.py",
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
    if cmd and any(pattern in cmd for pattern in patterns):
        print(f"{pid} {cmd}")
PY
if [[ -s /tmp/primary_qc_existing_processes.txt ]]; then
  echo "REFUSE: primary QC-related process already active" >&2
  cat /tmp/primary_qc_existing_processes.txt >&2
  exit 56
fi

run_id="$(date -u +%Y%m%dT%H%M%SZ)_primary_qc_memory_safe_v1"
out="$qc_parent/$run_id"
if [[ -e "$out" ]]; then
  echo "REFUSE: output root already exists: $out" >&2
  exit 57
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

record_resources() {
  local label="$1"
  {
    printf '%s\n' "$label"
    date -Is
    cat /sys/fs/cgroup/memory.current 2>/dev/null || true
    cat /sys/fs/cgroup/memory.max 2>/dev/null || true
    cat /sys/fs/cgroup/memory.events 2>/dev/null || true
    free -h || true
    df -h /root/autodl-tmp || true
  } >> "$out/runtime/resource_snapshots.log"
}

printf 'timestamp_utc\tstage\tstate\texit_code\n' > "$status"
write_status runner started 0
record_resources before_structural

set +e
"$py" scripts/run_primary_structure_memory_safe.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/structural_memory_safe.log" 2>&1
struct_code=$?
write_status structural_memory_safe finished "$struct_code"
set -e
record_resources after_structural

if [[ "$struct_code" -ne 0 ]]; then
  write_status runner failed_pre_energy "$struct_code"
  exit "$struct_code"
fi

set +e
"$py" scripts/gmx_energy_qc.py \
  --manifest "$manifest" \
  --output-root "$out" \
  --mode extract \
  > "$out/logs/energy.log" 2>&1
energy_code=$?
write_status energy finished "$energy_code"
set -e
record_resources after_energy

if [[ "$energy_code" -ne 0 ]]; then
  write_status runner failed_pre_membrane "$energy_code"
  exit "$energy_code"
fi

set +e
"$py" scripts/analyze_membrane_qc_mdanalysis.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/membrane.log" 2>&1
membrane_code=$?
write_status membrane finished "$membrane_code"
set -e
record_resources after_membrane

if [[ "$membrane_code" -ne 0 ]]; then
  write_status runner failed_pre_validation "$membrane_code"
  exit "$membrane_code"
fi

set +e
"$py" scripts/validate_primary_postprocessing.py \
  --manifest "$manifest" \
  --output-root "$out" \
  > "$out/logs/validate_primary.log" 2>&1
validate_code=$?
write_status validate_primary finished "$validate_code"
set -e
record_resources after_validate

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
  "failed_prior_run_id": "$failed_run_id",
  "output_root": "$out",
  "runner": "$runner",
  "runner_log": "$runner_log",
  "status_tsv": "$runner_status",
  "policy": "memory_safe_structural_then_fixed_energy_then_membrane_then_primary_validation",
  "raw_trajectory_policy": "read_only_no_overwrite",
  "secondary_endpoint_energy_policy": "closed_until_primary_passes",
  "cgroup_memory_max_bytes": "$(cat /sys/fs/cgroup/memory.max 2>/dev/null || true)",
  "cgroup_memory_current_bytes": "$(cat /sys/fs/cgroup/memory.current 2>/dev/null || true)",
  "launched_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
sha256sum "$runner_json" > "$runner_json.sha256"

O6U_PRIMARY_QC_OUT="$out" O6U_PRIMARY_QC_STATUS="$runner_status" \
  nohup bash "$runner" > "$runner_log" 2>&1 &
runner_pid="$!"

cat <<EOF
PRIMARY_QC_MEMORY_SAFE_RECOVERY_LAUNCHED
run_id=$run_id
pid=$runner_pid
output_root=$out
runner_log=$runner_log
status_tsv=$runner_status
manifest_sha256=$actual_manifest_sha
policy=memory_safe_structural_then_fixed_energy_then_membrane_then_primary_validation
EOF
