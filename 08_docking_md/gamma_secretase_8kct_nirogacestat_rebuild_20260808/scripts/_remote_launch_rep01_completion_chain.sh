#!/usr/bin/env bash
set -euo pipefail

root="/root/autodl-tmp/o6u_md_release_3x500ns_v4"
manifest="$root/config/primary_postprocessing_manifest.approved.json"
py="/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python"
qc_parent="$root/analysis/primary_postprocessing"
failed_run="20260822T162655Z_primary_qc_parallel_v1"
runtime_dir="$root/audit/postproduction_runtime"

run_id="$(date -u +%Y%m%dT%H%M%SZ)_primary_qc_rep01_completion"
out="$qc_parent/$run_id"
status="$runtime_dir/run_${run_id}.status.tsv"

mkdir -p "$out/logs" "$out/runtime" "$out/structural_analysis"
cd "$root"

write_status() { local stage="$1" state="$2" code="$3" ts; ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; printf '%s\t%s\t%s\t%s\n' "$ts" "$stage" "$state" "$code" >> "$status"; }
printf 'timestamp_utc\tstage\tstate\texit_code\n' > "$status"
write_status runner started 0

# ---- Guard: manifest sha ----
actual_manifest_sha="$(sha256sum "$manifest" | awk '{print $1}')"
if [[ "$actual_manifest_sha" != "0862ac9b1cad7f75a27e0dbda61d6eb4757f312c15ed7952bbe1da55e93818cb" ]]; then
  echo "REFUSE: manifest SHA drift $actual_manifest_sha" >&2; write_status runner refused_manifest 51; exit 51
fi

# ---- 1. Re-run rep01 structural child (memory-safe child mode) ----
set +e
"$py" scripts/run_primary_structure_memory_safe.py \
  --manifest "$manifest" --output-root "$out" --child-realization rep01 \
  > "$out/logs/structural_rep01_completion.log" 2>&1
rep01_code=$?
set -e
write_status structural_rep01_completion finished "$rep01_code"
if [[ "$rep01_code" -ne 0 ]]; then write_status runner failed_rep01 "$rep01_code"; exit "$rep01_code"; fi

# ---- 2. Copy rep02/rep03 summaries from the completed parallel run ----
for rep in rep02 rep03; do
  src="$qc_parent/$failed_run/structural_analysis/$rep/structural_summary.json"
  if [[ ! -s "$src" ]]; then echo "REFUSE: missing completed summary $src" >&2; write_status runner missing_summary 58; exit 58; fi
  mkdir -p "$out/structural_analysis/$rep"
  cp "$src" "$out/structural_analysis/$rep/structural_summary.json"
done
write_status copy_rep02_rep03 ok 0

# ---- 3. Assemble COMPLETE.json ----
O6U_OUT="$out" "$py" - <<'PY' >> "$out/logs/assemble.log" 2>&1
import os, sys
from pathlib import Path
sys.path.insert(0, "scripts")
from primary_postprocessing_common import validate_primary_manifest
from run_primary_structure_memory_safe import load_realization_summaries, assemble_complete_report
out = Path(os.environ["O6U_OUT"])
m = Path("/root/autodl-tmp/o6u_md_release_3x500ns_v4/config/primary_postprocessing_manifest.approved.json")
validate_primary_manifest(m)
summaries = load_realization_summaries(out)
complete = assemble_complete_report(m, out, summaries)
print("ASSEMBLED status=", complete["status"])
PY
assemble_code=$?
write_status assemble finished "$assemble_code"
if [[ "$assemble_code" -ne 0 ]]; then write_status runner failed_assemble "$assemble_code"; exit "$assemble_code"; fi

# ---- 4. Energy QC ----
set +e
"$py" scripts/gmx_energy_qc.py --manifest "$manifest" --output-root "$out" --mode extract > "$out/logs/energy.log" 2>&1
energy_code=$?
set -e
write_status energy finished "$energy_code"
if [[ "$energy_code" -ne 0 ]]; then write_status runner failed_pre_membrane "$energy_code"; exit "$energy_code"; fi

# ---- 5. Membrane QC ----
set +e
"$py" scripts/analyze_membrane_qc_mdanalysis.py --manifest "$manifest" --output-root "$out" > "$out/logs/membrane.log" 2>&1
membrane_code=$?
set -e
write_status membrane finished "$membrane_code"
if [[ "$membrane_code" -ne 0 ]]; then write_status runner failed_pre_validation "$membrane_code"; exit "$membrane_code"; fi

# ---- 6. Primary validation ----
set +e
"$py" scripts/validate_primary_postprocessing.py --manifest "$manifest" --output-root "$out" > "$out/logs/validate_primary.log" 2>&1
validate_code=$?
set -e
write_status validate_primary finished "$validate_code"
if [[ "$validate_code" -ne 0 ]]; then write_status runner validation_failed "$validate_code"; exit "$validate_code"; fi

write_status runner complete 0

cat <<EOF
REP01_COMPLETION_CHAIN_LAUNCHED
run_id=$run_id
output_root=$out
status_tsv=$status
manifest_sha256=$actual_manifest_sha
rep01_structural_rc=$rep01_code
energy_rc=$energy_code
membrane_rc=$membrane_code
validate_rc=$validate_code
EOF
