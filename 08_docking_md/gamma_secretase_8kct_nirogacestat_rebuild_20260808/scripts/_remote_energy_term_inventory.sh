set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
GMX=/root/GROMACS-2025.2/bin/gmx
OUT="$ROOT/audit/energy_term_inventory"
mkdir -p "$OUT"
"$GMX" --version > "$OUT/gmx_version.txt" 2>&1
printf '0\n' | "$GMX" energy -f "$ROOT/rep01/work/production.edr" -o "$OUT/no_selection.xvg" \
  > "$OUT/gmx_energy_menu.stdout" 2> "$OUT/gmx_energy_menu.stderr"
sha256sum "$ROOT/rep01/work/production.edr" "$OUT/gmx_version.txt" "$OUT/gmx_energy_menu.stdout" "$OUT/gmx_energy_menu.stderr" > "$OUT/SHA256SUMS"
echo '--- stdout tail ---'
tail -n 100 "$OUT/gmx_energy_menu.stdout"
echo '--- stderr tail ---'
tail -n 40 "$OUT/gmx_energy_menu.stderr"
echo '--- hashes ---'
cat "$OUT/SHA256SUMS"
