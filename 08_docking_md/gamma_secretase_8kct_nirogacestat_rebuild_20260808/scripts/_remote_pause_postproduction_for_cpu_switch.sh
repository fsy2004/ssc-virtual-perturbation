#!/bin/bash
set -euo pipefail

root=/root/autodl-tmp/o6u_md_release_3x500ns_v4
runtime="$root/audit/postproduction_runtime"
install_audit="$root/audit/toolchain_install/20260822"
stamp=$(date '+%Y%m%dT%H%M%S%z')
pause_dir="$root/audit/cpu_instance_switch_pause/$stamp"
mkdir -p "$pause_dir"

declare -A expected
expected[rep01_primary_pbc]='prepare_primary_pbc_trajectories.py --release-root /root/autodl-tmp/o6u_md_release_3x500ns_v4 --replica rep01'
expected[gmx_mmpbsa_install]='install_gmx_mmpbsa_1_6_5_cpu.sh /root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822'
expected[gorder_install]='install_gorder_1_5_0.sh /root/autodl-tmp/tools/gorder_1_5_0_20260822'
expected[fatslim_install]='install_fatslim_0_2_2.sh /root/autodl-tmp/tools/fatslim_0_2_2_20260822'

declare -A pidfiles
pidfiles[rep01_primary_pbc]="$runtime/rep01_primary_pbc.pid"
pidfiles[gmx_mmpbsa_install]="$install_audit/gmx_mmpbsa_install.pid"
pidfiles[gorder_install]="$install_audit/gorder_install.pid"
pidfiles[fatslim_install]="$install_audit/fatslim_install.pid"

descendants() {
  local parent=$1 child
  for child in $(pgrep -P "$parent" 2>/dev/null || true); do
    descendants "$child"
    printf '%s\n' "$child"
  done
}

printf 'pause_reason=cpu_instance_switch\ncreated_at=%s\n' "$(date --iso-8601=seconds)" > "$pause_dir/PAUSE_METADATA.txt"
ps -eo pid,ppid,lstart,etime,stat,%cpu,%mem,rss,args > "$pause_dir/processes_before.txt"
for log in \
  "$runtime/rep01_primary_pbc.log" \
  "$install_audit/gmx_mmpbsa_install.log" \
  "$install_audit/gorder_install.log" \
  "$install_audit/fatslim_install.log"; do
  if [[ -f "$log" ]]; then
    name=$(basename "$log")
    tail -n 200 "$log" > "$pause_dir/${name}.tail.txt" || true
  fi
done

find "$root/analysis/trajectories/8kct_nirogacestat_native/rep01" -maxdepth 1 -type f \
  -printf '%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS%Tz\n' 2>/dev/null \
  | sort > "$pause_dir/rep01_derived_file_inventory.tsv"
for prefix in \
  /root/autodl-tmp/envs/gmx_mmpbsa_1_6_5_cpu_20260822 \
  /root/autodl-tmp/tools/gorder_1_5_0_20260822 \
  /root/autodl-tmp/tools/fatslim_0_2_2_20260822; do
  if [[ -e "$prefix" ]]; then
    du -sh "$prefix" || true
  else
    printf 'MISSING\t%s\n' "$prefix"
  fi
done > "$pause_dir/toolchain_prefix_sizes.txt"

parents=()
all_targets=()
for name in rep01_primary_pbc gmx_mmpbsa_install gorder_install fatslim_install; do
  pidfile=${pidfiles[$name]}
  [[ -f "$pidfile" ]] || continue
  pid=$(tr -d '[:space:]' < "$pidfile")
  [[ "$pid" =~ ^[0-9]+$ ]] || { echo "invalid PID for $name: $pid" >&2; exit 10; }
  if ! kill -0 "$pid" 2>/dev/null; then
    continue
  fi
  cmd=$(ps -o args= -p "$pid")
  [[ "$cmd" == *"${expected[$name]}"* ]] || {
    echo "PID command mismatch for $name: $pid $cmd" >&2
    exit 11
  }
  mapfile -t children < <(descendants "$pid")
  for child in "${children[@]}"; do
    [[ -n "$child" ]] && all_targets+=("$child")
  done
  parents+=("$pid")
done

for pid in "${all_targets[@]}"; do kill -TERM "$pid" 2>/dev/null || true; done
for pid in "${parents[@]}"; do kill -TERM "$pid" 2>/dev/null || true; done

for _ in $(seq 1 20); do
  alive=0
  for pid in "${all_targets[@]}" "${parents[@]}"; do
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then alive=1; fi
  done
  [[ "$alive" -eq 0 ]] && break
  sleep 1
done

forced=()
for pid in "${all_targets[@]}" "${parents[@]}"; do
  if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
    kill -KILL "$pid" 2>/dev/null || true
    forced+=("$pid")
  fi
done
printf '%s\n' "${forced[@]:-}" > "$pause_dir/forced_kill_pids.txt"
sleep 1
ps -eo pid,ppid,lstart,etime,stat,%cpu,%mem,rss,args > "$pause_dir/processes_after.txt"
sync

remaining=$(ps -eo pid,args | grep -E 'prepare_primary_pbc_trajectories|install_gmx_mmpbsa_1_6_5_cpu|install_gorder_1_5_0|install_fatslim_0_2_2|gmx mindist|conda create|rustup-init' | grep -v grep || true)
if [[ -n "$remaining" ]]; then
  printf '%s\n' "$remaining" > "$pause_dir/unexpected_remaining_processes.txt"
  echo 'pause incomplete; matching processes remain' >&2
  cat "$pause_dir/unexpected_remaining_processes.txt" >&2
  exit 12
fi

sha256sum "$pause_dir"/*.txt "$pause_dir"/*.tsv > "$pause_dir/PAUSE_EVIDENCE_SHA256.txt"
echo "PAUSED pause_dir=$pause_dir"
