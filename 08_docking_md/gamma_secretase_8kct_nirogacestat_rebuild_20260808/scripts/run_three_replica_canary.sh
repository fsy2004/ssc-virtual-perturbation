#!/usr/bin/env bash
set -euo pipefail

root="${1:?usage: run_three_replica_canary.sh RELEASE_ROOT [GMX_BIN]}"
gmx_bin="${2:-gmx}"
canary_steps=5000

for rep in rep01 rep02 rep03; do
  work="$root/$rep/work"
  out="$root/canary/$rep"
  test -s "$work/step6.1_equilibration.tpr"
  test ! -e "$out"
  mkdir -p "$out"
  "$gmx_bin" convert-tpr -s "$work/step6.1_equilibration.tpr" \
    -o "$out/canary.tpr" -nsteps "$canary_steps" \
    > "$out/convert_tpr.stdout" 2> "$out/convert_tpr.stderr"
done

pids=()
offsets=(0 8 16)
replicas=(rep01 rep02 rep03)
for index in "${!replicas[@]}"; do
  rep="${replicas[$index]}"
  offset="${offsets[$index]}"
  out="$root/canary/$rep"
  (
    cd "$out"
    "$gmx_bin" mdrun -deffnm canary -ntmpi 1 -ntomp 8 -pin on -pinoffset "$offset" \
      > mdrun.stdout 2> mdrun.stderr
  ) &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
test "$status" -eq 0

# Do not flag configuration fields such as "epsilon-rf = inf". Numerical
# finiteness is validated from the energy records and completed coordinates.
pattern='LINCS WARNING|Too many LINCS warnings|constraint warning|SETTLE.*(error|constraint)|(^|[^A-Za-z])NaN([^A-Za-z]|$)|Fatal error|Segmentation fault'
for rep in rep01 rep02 rep03; do
  out="$root/canary/$rep"
  test -s "$out/canary.gro"
  test -s "$out/canary.edr"
  test -s "$out/canary.log"
  if grep -Ein "$pattern" "$out/canary.log" "$out/mdrun.stderr" > "$out/blocking_scan.txt"; then
    echo "blocking canary warning in $rep" >&2
    exit 3
  fi
  sha256sum "$out/canary.tpr" "$out/canary.gro" "$out/canary.edr" "$out/canary.log" \
    > "$out/SHA256SUMS.txt"
done

touch "$root/canary/COMPLETE"
