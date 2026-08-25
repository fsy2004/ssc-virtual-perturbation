#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
run_id="20260822T144701Z_primary_qc_v1"
out="$root/analysis/primary_postprocessing/$run_id"
runtime="$root/audit/postproduction_runtime"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
expected_manifest_sha="0862ac9b1cad7f75a27e0dbda61d6eb4757f312c15ed7952bbe1da55e93818cb"
manifest="$root/config/primary_postprocessing_manifest.approved.json"
waiter="$runtime/recover_${run_id}_after_legend_fix.sh"
waiter_log="$runtime/recover_${run_id}_after_legend_fix.log"
waiter_status="$runtime/recover_${run_id}_after_legend_fix.status.tsv"

mkdir -p "$runtime"
cd "$root"

actual_manifest_sha="$(sha256sum "$manifest" | awk '{print $1}')"
if [[ "$actual_manifest_sha" != "$expected_manifest_sha" ]]; then
  echo "REFUSE: approved manifest SHA drift before recovery: $actual_manifest_sha" >&2
  exit 21
fi

"$py" - <<'PY' >/tmp/primary_qc_recovery_existing.txt
import os
from pathlib import Path

needle = "recover_20260822T144701Z_primary_qc_v1_after_legend_fix.sh"
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
    if needle in cmd:
        print(f"{pid} {cmd}")
PY
if [[ -s /tmp/primary_qc_recovery_existing.txt ]]; then
  echo "REFUSE: recovery waiter already active" >&2
  cat /tmp/primary_qc_recovery_existing.txt >&2
  exit 22
fi

cat > "$waiter" <<'WAITER'
#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
run_id="20260822T144701Z_primary_qc_v1"
out="$root/analysis/primary_postprocessing/$run_id"
runtime="$root/audit/postproduction_runtime"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
status="$runtime/run_${run_id}.status.tsv"
recovery_status="$runtime/recover_${run_id}_after_legend_fix.status.tsv"

write_status() {
  local stage="$1"
  local state="$2"
  local code="$3"
  local ts
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '%s\t%s\t%s\t%s\n' "$ts" "$stage" "$state" "$code" >> "$recovery_status"
}

printf 'timestamp_utc\tstage\tstate\texit_code\n' > "$recovery_status"
write_status waiter started 0

cd "$root"
"$py" scripts/gmx_energy_qc.py --self-test > "$out/logs/energy_parser_fix_selftest.log" 2>&1
write_status parser_fix_selftest pass 0

for _ in $(seq 1 720); do
  if ! pgrep -f "run_${run_id}.sh" >/dev/null 2>&1; then
    break
  fi
  sleep 30
done

if pgrep -f "run_${run_id}.sh" >/dev/null 2>&1; then
  write_status original_runner still_active_timeout 40
  exit 40
fi
write_status original_runner exited 0

if grep -q $'\trunner\tcomplete\t0' "$status"; then
  write_status no_recovery_needed original_complete 0
  exit 0
fi

if ! grep -q $'\tstructural\tfinished\t0' "$status"; then
  write_status refuse structural_not_passed 41
  exit 41
fi

if ! grep -q $'\tenergy\tfinished\t2' "$status"; then
  write_status refuse unexpected_energy_status 42
  exit 42
fi

if ! grep -q "Kinetic-En." "$out/logs/energy.log" || ! grep -q "Kinetic En." "$out/logs/energy.log"; then
  write_status refuse unexpected_energy_failure 43
  exit 43
fi

archive="$root/audit/primary_qc_recovery/$(date -u +%Y%m%dT%H%M%SZ)_kinetic_legend_parser_fix"
mkdir -p "$archive"
cp -a "$out/logs/energy.log" "$archive/energy_failed_before_parser_fix.log"
cp -a "$status" "$archive/original_runner_status.tsv"
if [[ -d "$out/energy_qc" ]]; then
  mv "$out/energy_qc" "$archive/energy_qc_failed_before_parser_fix"
fi
sha256sum "$archive"/* > "$archive/SHA256SUMS.txt" 2>/dev/null || true
write_status failure_evidence_archived "$archive" 0

set +e
"$py" scripts/gmx_energy_qc.py \
  --manifest "$root/config/primary_postprocessing_manifest.approved.json" \
  --output-root "$out" \
  --mode extract \
  > "$out/logs/energy_after_legend_fix.log" 2>&1
energy_code=$?
write_status energy_after_legend_fix finished "$energy_code"
set -e

if [[ "$energy_code" -ne 0 ]]; then
  write_status recovery_failed_pre_membrane "$energy_code"
  exit "$energy_code"
fi

set +e
"$py" scripts/analyze_membrane_qc_mdanalysis.py \
  --manifest "$root/config/primary_postprocessing_manifest.approved.json" \
  --output-root "$out" \
  > "$out/logs/membrane.log" 2>&1
membrane_code=$?
write_status membrane finished "$membrane_code"
set -e

if [[ "$membrane_code" -ne 0 ]]; then
  write_status recovery_failed_pre_validation "$membrane_code"
  exit "$membrane_code"
fi

set +e
"$py" scripts/validate_primary_postprocessing.py \
  --manifest "$root/config/primary_postprocessing_manifest.approved.json" \
  --output-root "$out" \
  > "$out/logs/validate_primary.log" 2>&1
validate_code=$?
write_status validate_primary finished "$validate_code"
set -e

if [[ "$validate_code" -ne 0 ]]; then
  write_status recovery_validation_failed "$validate_code"
  exit "$validate_code"
fi

write_status recovery complete 0
WAITER

chmod 700 "$waiter"
nohup bash "$waiter" > "$waiter_log" 2>&1 &
pid="$!"

cat <<EOF
PRIMARY_QC_RECOVERY_WAITER_LAUNCHED
run_id=$run_id
pid=$pid
waiter=$waiter
waiter_log=$waiter_log
waiter_status=$waiter_status
policy=wait_for_original_runner_then_resume_only_if_structural_passed_and_energy_failure_matches_legend_parser
EOF
