set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
pgrep -af 'gmx energy.*production.edr|energy_term_inventory' || true
find "$ROOT/audit/energy_term_inventory" -maxdepth 1 -type f -printf '%f %s bytes\n' | sort 2>/dev/null || true
tail -n 120 "$ROOT/audit/energy_term_inventory/gmx_energy_menu.stdout" 2>/dev/null || true
tail -n 30 "$ROOT/audit/energy_term_inventory/gmx_energy_menu.stderr" 2>/dev/null || true
